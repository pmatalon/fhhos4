#pragma once
#include "../../Mesh/Mesh.h"
#include "../../Utils/Utils.h"
#include "../../TestCases/Diffusion/DiffusionTestCase.h"
#include "../../Geometry/CartesianShape.h"
#include "../../Geometry/2D/Triangle.h"
#include "../../Utils/Parallelism.h"
#include "../../Utils/ExportModule.h"
#include "Diff_DGElement.h"
#include "Diff_DGFace.h"
using namespace std;

template <int Dim>
class Diffusion_DG
{
public:
	Mesh<Dim>* _mesh;
	FunctionalBasis<Dim>* Basis;
	SparseMatrix A;
	Vector b;
	Vector SystemSolution;
private:
	DiffusionTestCase<Dim>* _testCase;
	bool _autoPenalization;
	int _penalizationCoefficient;
public:
	Diffusion_DG(Mesh<Dim>* mesh, DiffusionTestCase<Dim>* testCase, string outputDirectory, FunctionalBasis<Dim>* basis, int penalizationCoefficient)
	{ 
		this->_mesh = mesh;
		this->_testCase = testCase;
		this->Basis = basis;
		this->_autoPenalization = penalizationCoefficient == -1;
		if (_autoPenalization)
			this->_penalizationCoefficient = pow(Dim, 2) * pow(Basis->GetDegree() + 1, 2) / this->_mesh->H(); // Ralph-Hartmann
		else
			this->_penalizationCoefficient = penalizationCoefficient;

		//Problem<Dim>::AddFilePrefix("_DG_SIPG_" + basis->Name() + "_pen" + to_string(penalizationCoefficient));
	}

	void PrintDiscretization()
	{
		cout << "Mesh: " << this->_mesh->Description() << endl;
		cout << "Discretization: Discontinuous Galerkin SIPG" << endl;
		cout << "\tPolynomial space: " << (Basis->UsePolynomialSpaceQ() ? "Q" : "P") << endl;
		cout << "\tPolynomial basis: " << Basis->Name() << endl;
		cout << "\tPenalization coefficient: " << _penalizationCoefficient << (_autoPenalization ? " (automatic)" : "") << endl;
		cout << "Local functions: " << Basis->Size() << endl;
		for (BasisFunction<Dim>* phi : Basis->LocalFunctions())
			cout << "\t " << phi->ToString() << endl;
		BigNumber nUnknowns = static_cast<int>(this->_mesh->Elements.size()) * Basis->Size();
		cout << "Unknowns   : " << nUnknowns << endl;
	}

	double L2Error(DomFunction exactSolution)
	{
		double absoluteError = 0;
		double normExactSolution = 0;

		#pragma omp parallel for reduction(+:absoluteError, normExactSolution)
		for (Element<Dim>* element : _mesh->Elements)
		{
			auto approximate = Basis->GetApproximateFunction(SystemSolution, element->Number * Basis->Size());
			absoluteError += element->L2ErrorPow2(approximate, exactSolution);
			normExactSolution += element->Integral([exactSolution](const DomPoint& p) { return pow(exactSolution(p), 2); });
		}

		absoluteError = sqrt(absoluteError);
		normExactSolution = sqrt(normExactSolution);
		return normExactSolution != 0 ? absoluteError / normExactSolution : absoluteError;
	}

	void Assemble(const ActionsArguments& actions, const ExportModule& out) //override
	{
		auto mesh = this->_mesh;
		auto basis = this->Basis;
		auto penalizationCoefficient = this->_penalizationCoefficient;

		if (actions.LogAssembly)
			this->PrintDiscretization();


		BigNumber nUnknowns = static_cast<int>(mesh->Elements.size()) * basis->Size();
		this->b = Vector(nUnknowns);

		if (actions.LogAssembly)
		{
			cout << "--------------------------------------------------------" << endl;
			cout << "Assembly..." << endl;
		}

		CartesianShape<Dim, Dim>::InitReferenceShape()->ComputeAndStoreMassMatrix(basis);
		CartesianShape<Dim, Dim>::InitReferenceShape()->ComputeAndStoreStiffnessMatrix(basis);
		if (Dim == 2)
			Triangle::InitReferenceShape()->ComputeAndStoreMassMatrix((FunctionalBasis<2>*)basis);
		
		//--------------------------------------------//
		// Iteration on the elements: diagonal blocks //
		//--------------------------------------------//

		BigNumber nnzPerElement = basis->Size() * (2 * Dim + 1);
		BigNumber nnzPerElementForExport = actions.Export.AssemblyTermMatrices ? nnzPerElement : 0;
		ThreadLocalCoeffs matrixCoeffs(mesh->Elements.size(), nnzPerElement);
		ThreadLocalCoeffs massMatrixCoeffs(mesh->Elements.size(), nnzPerElementForExport);
		ThreadLocalCoeffs volumicCoeffs(mesh->Elements.size(), nnzPerElementForExport);
		ThreadLocalCoeffs couplingCoeffs(mesh->Elements.size(), nnzPerElementForExport);
		ThreadLocalCoeffs penCoeffs(mesh->Elements.size(), nnzPerElementForExport);

		#pragma omp parallel for
		for (Element<Dim>* e : mesh->Elements)
		{
			Diff_DGElement<Dim>* element = dynamic_cast<Diff_DGElement<Dim>*>(e);
			//cout << "Element " << element->Number << endl;

			for (BasisFunction<Dim>* phi1 : basis->LocalFunctions())
			{
				BigNumber basisFunction1 = element->Number * basis->Size() + phi1->LocalNumber;

				// Current element (block diagonal)
				for (BasisFunction<Dim>* phi2 : basis->LocalFunctions())
				{
					BigNumber basisFunction2 = element->Number * basis->Size() + phi2->LocalNumber;

					//cout << "\t phi" << phi1->LocalNumber << " = " << phi1->ToString() << " phi" << phi2->LocalNumber << " = " << phi2->ToString() << endl;

					double volumicTerm = element->VolumicTerm(phi1, phi2);
					//cout << "\t\t volumic = " << volumicTerm << endl;

					double coupling = 0;
					double penalization = 0;
					for (Face<Dim>* f : element->Faces)
					{
						Diff_DGFace<Dim>* face = dynamic_cast<Diff_DGFace<Dim>*>(f);

						double c = face->CouplingTerm(element, phi1, element, phi2);
						double p = face->PenalizationTerm(element, phi1, element, phi2, penalizationCoefficient);
						coupling += c;
						penalization += p;
						//cout << "\t\t " << face->ToString() << ":\t c=" << c << "\tp=" << p << endl;
					}

					//cout << "\t\t TOTAL = " << volumicTerm + coupling + penalization << endl;

					if (actions.Export.AssemblyTermMatrices)
					{
						volumicCoeffs.Local().Add(basisFunction1, basisFunction2, volumicTerm);
						couplingCoeffs.Local().Add(basisFunction1, basisFunction2, coupling);
						penCoeffs.Local().Add(basisFunction1, basisFunction2, penalization);
					}
					matrixCoeffs.Local().Add(basisFunction1, basisFunction2, volumicTerm + coupling + penalization);
					if (actions.Export.AssemblyTermMatrices)
					{
						double massTerm = element->MassTerm(phi1, phi2);
						massMatrixCoeffs.Local().Add(basisFunction1, basisFunction2, massTerm);
					}
				}

				double rhs = element->SourceTerm(phi1, _testCase->SourceFunction);
				this->b(basisFunction1) = rhs;
			}
		}

		//---------------------------------------------//
		// Iteration on the faces: off-diagonal blocks //
		//---------------------------------------------//

		#pragma omp parallel for
		for (Face<Dim>* f : mesh->Faces)
		{
			Diff_DGFace<Dim>* face = dynamic_cast<Diff_DGFace<Dim>*>(f);
			if (face->IsDomainBoundary)
				continue;

			//cout << "Face " << face->Number << endl;

			for (BasisFunction<Dim>* phi1 : basis->LocalFunctions())
			{
				BigNumber basisFunction1 = face->Element1->Number * basis->Size() + phi1->LocalNumber;
				for (BasisFunction<Dim>* phi2 : basis->LocalFunctions())
				{
					//cout << "\t phi" << phi1->LocalNumber << " = " << phi1->ToString() << " phi" << phi2->LocalNumber << " = " << phi2->ToString() << endl;

					BigNumber basisFunction2 = face->Element2->Number * basis->Size() + phi2->LocalNumber;
					double coupling = face->CouplingTerm(face->Element1, phi1, face->Element2, phi2);
					double penalization = face->PenalizationTerm(face->Element1, phi1, face->Element2, phi2, penalizationCoefficient);

					//cout << "\t\t\t c=" << coupling << "\tp=" << penalization << endl;

					if (actions.Export.AssemblyTermMatrices)
					{
						couplingCoeffs.Local().Add(basisFunction1, basisFunction2, coupling);
						couplingCoeffs.Local().Add(basisFunction2, basisFunction1, coupling);

						penCoeffs.Local().Add(basisFunction1, basisFunction2, penalization);
						penCoeffs.Local().Add(basisFunction2, basisFunction1, penalization);
					}
					matrixCoeffs.Local().Add(basisFunction1, basisFunction2, coupling + penalization);
					matrixCoeffs.Local().Add(basisFunction2, basisFunction1, coupling + penalization);
				}
			}
		}

		//---------------//
		// Matrix export //
		//---------------//

		this->A = SparseMatrix(nUnknowns, nUnknowns);
		matrixCoeffs.Fill(this->A);
		cout << "nnz(A) = " << this->A.nonZeros() << endl;

		if (actions.Export.LinearSystem)
		{
			cout << "Export of the linear system..." << endl;
			out.ExportMatrix(this->A, "A");

			out.ExportVector(this->b, "b");
		}

		if (actions.Export.AssemblyTermMatrices)
		{
			SparseMatrix M(nUnknowns, nUnknowns);
			massMatrixCoeffs.Fill(M);
			out.ExportMatrix(M, "Mass");

			SparseMatrix V(nUnknowns, nUnknowns);
			volumicCoeffs.Fill(V);
			out.ExportMatrix(V, "A_volumic");

			SparseMatrix C(nUnknowns, nUnknowns);
			couplingCoeffs.Fill(C);
			out.ExportMatrix(C, "A_coupling");

			SparseMatrix P(nUnknowns, nUnknowns);
			penCoeffs.Fill(P);
			out.ExportMatrix(P, "A_pen");
		}

	}
};

