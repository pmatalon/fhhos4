#pragma once
#include "../Problem/DiffusionProblem.h"
#include "../Geometry/CartesianShape.h"
#include "../Geometry/2D/Triangle.h"
#include "../Utils/ParallelLoop.h"
#include "Diff_DGElement.h"
#include "Diff_DGFace.h"
using namespace std;

template <int Dim>
class Diffusion_DG : public DiffusionProblem<Dim>
{
public:
	FunctionalBasis<Dim>* Basis;
private:
	bool _autoPenalization;
	int _penalizationCoefficient;
public:
	Diffusion_DG(Mesh<Dim>* mesh, TestCase<Dim>* testCase, string outputDirectory, FunctionalBasis<Dim>* basis, int penalizationCoefficient)
		: DiffusionProblem<Dim>(mesh, testCase, outputDirectory)
	{ 
		this->Basis = basis;
		this->_autoPenalization = penalizationCoefficient == -1;
		if (_autoPenalization)
			this->_penalizationCoefficient = pow(Dim, 2) * pow(Basis->GetDegree() + 1, 2) / this->_mesh->H(); // Ralph-Hartmann
		else
			this->_penalizationCoefficient = penalizationCoefficient;

		this->_fileName = this->_fileName + "_DG_SIPG_" + basis->Name() + "_pen" + to_string(penalizationCoefficient);
	}

	void PrintDiscretization()
	{
		cout << "Mesh: " << this->_mesh->Description() << endl;
		cout << "Discretization: Discontinuous Galerkin SIPG" << endl;
		cout << "\tPolynomial space: " << (Basis->UsePolynomialSpaceQ ? "Q" : "P") << endl;
		cout << "\tPolynomial basis: " << Basis->Name() << endl;
		cout << "\tPenalization coefficient: " << _penalizationCoefficient << (_autoPenalization ? " (automatic)" : "") << endl;
		cout << "Local functions: " << Basis->NumberOfLocalFunctionsInElement(NULL) << endl;
		for (BasisFunction<Dim>* phi : Basis->LocalFunctions)
			cout << "\t " << phi->ToString() << endl;
		BigNumber nUnknowns = static_cast<int>(this->_mesh->Elements.size()) * Basis->NumberOfLocalFunctionsInElement(NULL);
		cout << "Unknowns   : " << nUnknowns << endl;
	}

	double L2Error(DomFunction exactSolution) override
	{
		return Problem<Dim>::L2Error(Basis, this->SystemSolution, exactSolution);
	}

	void Assemble(ActionsArguments actions) override
	{
		auto mesh = this->_mesh;
		auto basis = this->Basis;
		auto penalizationCoefficient = this->_penalizationCoefficient;

		if (actions.LogAssembly)
		{
			this->PrintPhysicalProblem();
			this->PrintDiscretization();
			cout << "Parallelism: " << (BaseParallelLoop::GetDefaultNThreads() == 1 ? "sequential execution" : to_string(BaseParallelLoop::GetDefaultNThreads()) + " threads") << endl;
		}

		string matrixFilePath			= this->GetFilePath("A");
		string matrixVolumicFilePath	= this->GetFilePath("A_volumic");
		string matrixCouplingFilePath	= this->GetFilePath("A_coupling");
		string matrixPenFilePath		= this->GetFilePath("A_pen");
		string massMatrixFilePath		= this->GetFilePath("Mass");
		string rhsFilePath				= this->GetFilePath("b");

		BigNumber nUnknowns = static_cast<int>(mesh->Elements.size()) * basis->NumberOfLocalFunctionsInElement(NULL);
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

		ParallelLoop<Element<Dim>*, EmptyResultChunk> parallelLoop(mesh->Elements);

		vector<NonZeroCoefficients> chunksMatrixCoeffs(parallelLoop.NThreads);
		vector<NonZeroCoefficients> chunksMassMatrixCoeffs(parallelLoop.NThreads);
		vector<NonZeroCoefficients> chunksVolumicCoeffs(parallelLoop.NThreads);
		vector<NonZeroCoefficients> chunksCouplingCoeffs(parallelLoop.NThreads);
		vector<NonZeroCoefficients> chunksPenCoeffs(parallelLoop.NThreads);
		
		for (unsigned int threadNumber = 0; threadNumber < parallelLoop.NThreads; threadNumber++)
		{
			ParallelChunk<EmptyResultChunk>* chunk = parallelLoop.Chunks[threadNumber];

			chunk->ThreadFuture = std::async([this, mesh, basis, &actions, penalizationCoefficient, chunk, &chunksMatrixCoeffs, &chunksMassMatrixCoeffs, &chunksVolumicCoeffs, &chunksCouplingCoeffs, &chunksPenCoeffs]()
			{
				BigNumber nnzApproximate = chunk->Size() * basis->Size() * (2 * Dim + 1);
				NonZeroCoefficients matrixCoeffs(nnzApproximate);
				NonZeroCoefficients massMatrixCoeffs(actions.ExportAssemblyTermMatrices ? nnzApproximate : 0);
				NonZeroCoefficients volumicCoeffs(actions.ExportAssemblyTermMatrices ? nnzApproximate : 0);
				NonZeroCoefficients couplingCoeffs(actions.ExportAssemblyTermMatrices ? nnzApproximate : 0);
				NonZeroCoefficients penCoeffs(actions.ExportAssemblyTermMatrices ? nnzApproximate : 0);

				for (BigNumber iElem = chunk->Start; iElem < chunk->End; iElem++)
				{
					Diff_DGElement<Dim>* element = dynamic_cast<Diff_DGElement<Dim>*>(mesh->Elements[iElem]);
					//cout << "Element " << element->Number << endl;

					for (BasisFunction<Dim>* phi1 : basis->LocalFunctions)
					{
						BigNumber basisFunction1 = basis->GlobalFunctionNumber(element, phi1);

						// Current element (block diagonal)
						for (BasisFunction<Dim>* phi2 : basis->LocalFunctions)
						{
							BigNumber basisFunction2 = basis->GlobalFunctionNumber(element, phi2);

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

							if (actions.ExportAssemblyTermMatrices)
							{
								volumicCoeffs.Add(basisFunction1, basisFunction2, volumicTerm);
								couplingCoeffs.Add(basisFunction1, basisFunction2, coupling);
								penCoeffs.Add(basisFunction1, basisFunction2, penalization);
							}
							matrixCoeffs.Add(basisFunction1, basisFunction2, volumicTerm + coupling + penalization);
							if (actions.ExportAssemblyTermMatrices)
							{
								double massTerm = element->MassTerm(phi1, phi2);
								massMatrixCoeffs.Add(basisFunction1, basisFunction2, massTerm);
							}
						}

						double rhs = element->SourceTerm(phi1, this->_sourceFunction);
						this->b(basisFunction1) = rhs;
					}
				}

				chunksMatrixCoeffs[chunk->ThreadNumber] = matrixCoeffs;
				chunksMassMatrixCoeffs[chunk->ThreadNumber] = massMatrixCoeffs;
				chunksVolumicCoeffs[chunk->ThreadNumber] = volumicCoeffs;
				chunksCouplingCoeffs[chunk->ThreadNumber] = couplingCoeffs;
				chunksPenCoeffs[chunk->ThreadNumber] = penCoeffs;
			}
			);
		}

		BigNumber nnzApproximate = mesh->Elements.size() * basis->Size() * (2 * Dim + 1);
		NonZeroCoefficients matrixCoeffs(nnzApproximate);
		NonZeroCoefficients massMatrixCoeffs(actions.ExportAssemblyTermMatrices ? nnzApproximate : 0);
		NonZeroCoefficients volumicCoeffs(actions.ExportAssemblyTermMatrices ? nnzApproximate : 0);
		NonZeroCoefficients couplingCoeffs(actions.ExportAssemblyTermMatrices ? nnzApproximate : 0);
		NonZeroCoefficients penCoeffs(actions.ExportAssemblyTermMatrices ? nnzApproximate : 0);

		parallelLoop.Wait();

		for (unsigned int threadNumber = 0; threadNumber < parallelLoop.NThreads; threadNumber++)
		{
			matrixCoeffs.Add(chunksMatrixCoeffs[threadNumber]);
			massMatrixCoeffs.Add(chunksMassMatrixCoeffs[threadNumber]);
			volumicCoeffs.Add(chunksVolumicCoeffs[threadNumber]);
			couplingCoeffs.Add(chunksCouplingCoeffs[threadNumber]);
			penCoeffs.Add(chunksPenCoeffs[threadNumber]);
		}

		chunksMatrixCoeffs.clear();
		chunksMassMatrixCoeffs.clear();
		chunksVolumicCoeffs.clear();
		chunksCouplingCoeffs.clear();
		chunksPenCoeffs.clear();

		//---------------------------------------------//
		// Iteration on the faces: off-diagonal blocks //
		//---------------------------------------------//

		ParallelLoop<Face<Dim>*, EmptyResultChunk> parallelLoopFaces(mesh->Faces);

		chunksMatrixCoeffs = vector<NonZeroCoefficients>(parallelLoopFaces.NThreads);
		chunksCouplingCoeffs = vector<NonZeroCoefficients>(parallelLoopFaces.NThreads);
		chunksPenCoeffs = vector<NonZeroCoefficients>(parallelLoopFaces.NThreads);


		for (unsigned int threadNumber = 0; threadNumber < parallelLoopFaces.NThreads; threadNumber++)
		{
			ParallelChunk<EmptyResultChunk>* chunk = parallelLoopFaces.Chunks[threadNumber];

			chunk->ThreadFuture = std::async([this, mesh, basis, &actions, penalizationCoefficient, chunk, &chunksMatrixCoeffs, &chunksCouplingCoeffs, &chunksPenCoeffs]()
			{
				BigNumber nnzApproximate = chunk->Size() * basis->Size() * (2 * Dim + 1);
				NonZeroCoefficients matrixCoeffs(nnzApproximate);
				NonZeroCoefficients couplingCoeffs(actions.ExportAssemblyTermMatrices ? nnzApproximate : 0);
				NonZeroCoefficients penCoeffs(actions.ExportAssemblyTermMatrices ? nnzApproximate : 0);

				for (BigNumber iElem = chunk->Start; iElem < chunk->End; ++iElem)
				{
					Diff_DGFace<Dim>* face = dynamic_cast<Diff_DGFace<Dim>*>(mesh->Faces[iElem]);
					if (face->IsDomainBoundary)
						continue;

					//cout << "Face " << face->Number << endl;

					for (BasisFunction<Dim>* phi1 : basis->LocalFunctions)
					{
						BigNumber basisFunction1 = basis->GlobalFunctionNumber(face->Element1, phi1);
						for (BasisFunction<Dim>* phi2 : basis->LocalFunctions)
						{
							//cout << "\t phi" << phi1->LocalNumber << " = " << phi1->ToString() << " phi" << phi2->LocalNumber << " = " << phi2->ToString() << endl;

							BigNumber basisFunction2 = basis->GlobalFunctionNumber(face->Element2, phi2);
							double coupling = face->CouplingTerm(face->Element1, phi1, face->Element2, phi2);
							double penalization = face->PenalizationTerm(face->Element1, phi1, face->Element2, phi2, penalizationCoefficient);

							//cout << "\t\t\t c=" << coupling << "\tp=" << penalization << endl;

							if (actions.ExportAssemblyTermMatrices)
							{
								couplingCoeffs.Add(basisFunction1, basisFunction2, coupling);
								couplingCoeffs.Add(basisFunction2, basisFunction1, coupling);

								penCoeffs.Add(basisFunction1, basisFunction2, penalization);
								penCoeffs.Add(basisFunction2, basisFunction1, penalization);
							}
							matrixCoeffs.Add(basisFunction1, basisFunction2, coupling + penalization);
							matrixCoeffs.Add(basisFunction2, basisFunction1, coupling + penalization);
						}
					}
				}

				chunksMatrixCoeffs[chunk->ThreadNumber] = matrixCoeffs;
				chunksCouplingCoeffs[chunk->ThreadNumber] = couplingCoeffs;
				chunksPenCoeffs[chunk->ThreadNumber] = penCoeffs;
			}
			);
		}

		parallelLoopFaces.Wait();

		for (unsigned int threadNumber = 0; threadNumber < parallelLoopFaces.NThreads; threadNumber++)
		{
			matrixCoeffs.Add(chunksMatrixCoeffs[threadNumber]);
			couplingCoeffs.Add(chunksCouplingCoeffs[threadNumber]);
			penCoeffs.Add(chunksPenCoeffs[threadNumber]);
		}

		//---------------//
		// Matrix export //
		//---------------//

		this->A = SparseMatrix(nUnknowns, nUnknowns);
		matrixCoeffs.Fill(this->A);
		cout << "nnz(A) = " << this->A.nonZeros() << endl;

		if (actions.ExportLinearSystem)
		{
			cout << "Export..." << endl;
			Eigen::saveMarket(this->A, matrixFilePath);
			cout << "Matrix exported to \t" << matrixFilePath << endl;

			Eigen::saveMarketVector(this->b, rhsFilePath);
			cout << "RHS exported to \t" << rhsFilePath << endl;
		}

		if (actions.ExportAssemblyTermMatrices)
		{
			SparseMatrix M(nUnknowns, nUnknowns);
			massMatrixCoeffs.Fill(M);
			Eigen::saveMarket(M, massMatrixFilePath);
			cout << "Mass matrix exported to \t" << massMatrixFilePath << endl;

			SparseMatrix V(nUnknowns, nUnknowns);
			volumicCoeffs.Fill(V);
			Eigen::saveMarket(V, matrixVolumicFilePath);
			cout << "Volumic part exported to \t" << matrixVolumicFilePath << endl;

			SparseMatrix C(nUnknowns, nUnknowns);
			couplingCoeffs.Fill(C);
			Eigen::saveMarket(C, matrixCouplingFilePath);
			cout << "Coupling part exported to \t" << matrixCouplingFilePath << endl;

			SparseMatrix P(nUnknowns, nUnknowns);
			penCoeffs.Fill(P);
			Eigen::saveMarket(P, matrixPenFilePath);
			cout << "Penalization part exported to \t" << matrixPenFilePath << endl;
		}

	}

	void ExportSolutionToGMSH() override
	{
		this->_mesh->ExportToGMSH(this->Basis, this->SystemSolution, this->GetFilePathPrefix(), "potential");
	}

	void ExportErrorToGMSH(const Vector& coeffs) override
	{
		this->_mesh->ExportToGMSH(this->Basis, coeffs, this->GetFilePathPrefix(), "error");
	}
};

