#pragma once
#include <iostream>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include "../Problem.h"
#include "../Mesh/Mesh.h"
#include "../Mesh/Face.h"
#include "../DG/Poisson_DGTerms.h"
#include "../Utils/NonZeroCoefficients.h"
#include "../Utils/L2.h"
#include "../Utils/Action.h"
using namespace std;

template <short Dim>
class Poisson_DG : public Problem
{
private:
	DiffusionPartition _diffusionPartition;
public:

	Poisson_DG(string solutionName, DiffusionPartition diffusionPartition) 
		: _diffusionPartition(diffusionPartition), Problem(solutionName)
	{	}

	void Assemble(Mesh<Dim>* mesh, FunctionalBasis<Dim>* basis, Poisson_DGTerms<Dim>* dg, int penalizationCoefficient, string outputDirectory, Action action)
	{
		cout << "Problem: Poisson " << Dim << "D";
		if (this->_diffusionPartition.Kappa1 != this->_diffusionPartition.Kappa2)
			cout << ", heterogeneous diffusion coefficient (k1=" << this->_diffusionPartition.Kappa1 << ", k2=" << this->_diffusionPartition.Kappa2 << ")";
		cout << endl;
		cout << "Analytical solution: " ;
		if (this->_solutionName.compare("sine") == 0)
			cout << "sine function";
		else if (this->_solutionName.compare("sine") == 0)
			cout << "polynomial function";
		else if (this->_solutionName.compare("hetero") == 0)
			cout << "heterogeneous-specific piecewise polynomial function";
		else
			cout << "unknown";
		cout << endl;
		cout << "Subdivisions in each cartesian direction: " << mesh->N << endl;
		cout << "Discretization: Discontinuous Galerkin SIPG" << endl;
		cout << "\tPolynomial space: " << (basis->FullTensorization ? "Q" : "P") << endl;
		cout << "\tPolynomial basis: " << (dg->IsGlobalBasis() ? "global" : "") + basis->Name() << endl;

		bool autoPenalization = penalizationCoefficient == -1;
		if (autoPenalization)
			penalizationCoefficient = pow(Dim, 2) * pow(basis->GetDegree() + 1, 2) * mesh->N; // Ralph-Hartmann
		cout << "\tPenalization coefficient: " << penalizationCoefficient << (autoPenalization ? " (automatic)" : "") << endl;
		
		cout << "Local functions: " << basis->NumberOfLocalFunctionsInElement(NULL) << endl;
		for (BasisFunction<Dim>* phi : basis->LocalFunctions)
			cout << "\t " << phi->ToString() << endl;
		BigNumber nUnknowns = static_cast<int>(mesh->Elements.size()) * basis->NumberOfLocalFunctionsInElement(NULL);
		cout << "Unknowns: " << nUnknowns << endl;

		string fileName = "Poisson" + to_string(Dim) + "D" + this->_solutionName + "_n" + to_string(mesh->N) + "_DG_SIPG_" + (dg->IsGlobalBasis() ? "global" : "") + basis->Name() + "_pen" + (autoPenalization ? "-1" : to_string(penalizationCoefficient));
		string matrixFilePath			= outputDirectory + "/" + fileName + "_A.dat";
		string matrixVolumicFilePath	= outputDirectory + "/" + fileName + "_A_volumic.dat";
		string matrixCouplingFilePath	= outputDirectory + "/" + fileName + "_A_coupling.dat";
		string matrixPenFilePath		= outputDirectory + "/" + fileName + "_A_pen.dat";
		string massMatrixFilePath		= outputDirectory + "/" + fileName + "_Mass.dat";
		string rhsFilePath				= outputDirectory + "/" + fileName + "_b.dat";

		this->b = Eigen::VectorXd(nUnknowns);

		BigNumber nnzApproximate = mesh->Elements.size() * basis->NumberOfLocalFunctionsInElement(NULL) * (2 * Dim + 1);
		NonZeroCoefficients matrixCoeffs(nnzApproximate);
		NonZeroCoefficients massMatrixCoeffs((action & Action::ExtractMassMatrix) == Action::ExtractMassMatrix ? nnzApproximate : 0);
		NonZeroCoefficients volumicCoeffs((action & Action::ExtractComponentMatrices) == Action::ExtractComponentMatrices ? nnzApproximate : 0);
		NonZeroCoefficients couplingCoeffs((action & Action::ExtractComponentMatrices) == Action::ExtractComponentMatrices ? nnzApproximate : 0);
		NonZeroCoefficients penCoeffs((action & Action::ExtractComponentMatrices) == Action::ExtractComponentMatrices ? nnzApproximate : 0);

		cout << "Assembly..." << endl;

		//--------------------------------------------//
		// Iteration on the elements: diagonal blocks //
		//--------------------------------------------//

		for (auto element : mesh->Elements)
		{
			//cout << "Element " << element->Number << endl;

			for (BasisFunction<Dim>* phi1 : basis->LocalFunctions)
			{
				BigNumber basisFunction1 = basis->GlobalFunctionNumber(element, phi1);

				// Current element (block diagonal)
				for (BasisFunction<Dim>* phi2 : basis->LocalFunctions)
				{
					BigNumber basisFunction2 = basis->GlobalFunctionNumber(element, phi2);

					//cout << "\t phi" << phi1->LocalNumber << " = " << phi1->ToString() << " phi" << phi2->LocalNumber << " = " << phi2->ToString() << endl;

					double volumicTerm = dg->VolumicTerm(element, phi1, phi2);
					//cout << "\t\t volumic = " << volumicTerm << endl;
					
					double coupling = 0;
					double penalization = 0;
					for (Face<Dim>* face : element->Faces)
					{
						double c = dg->CouplingTerm(face, element, phi1, element, phi2);
						double p = dg->PenalizationTerm(face, element, phi1, element, phi2, penalizationCoefficient);
						coupling += c;
						penalization += p;
						//cout << "\t\t " << face->ToString() << ":\t c=" << c << "\tp=" << p << endl;
					}

					//cout << "\t\t TOTAL = " << volumicTerm + coupling + penalization << endl;
					
					if ((action & Action::ExtractComponentMatrices) == Action::ExtractComponentMatrices)
					{
						volumicCoeffs.Add(basisFunction1, basisFunction2, volumicTerm);
						couplingCoeffs.Add(basisFunction1, basisFunction2, coupling);
						penCoeffs.Add(basisFunction1, basisFunction2, penalization);
					}
					matrixCoeffs.Add(basisFunction1, basisFunction2, volumicTerm + coupling + penalization);
					if ((action & Action::ExtractMassMatrix) == Action::ExtractMassMatrix)
					{
						double massTerm = dg->MassTerm(element, phi1, phi2);
						massMatrixCoeffs.Add(basisFunction1, basisFunction2, massTerm);
					}
				}

				double rhs = dg->RightHandSide(element, phi1);
				this->b(basisFunction1) = rhs;
			}
		}

		//---------------------------------------------//
		// Iteration on the faces: off-diagonal blocks //
		//---------------------------------------------//

		for (auto face : mesh->Faces)
		{
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
					double coupling = dg->CouplingTerm(face, face->Element1, phi1, face->Element2, phi2);
					double penalization = dg->PenalizationTerm(face, face->Element1, phi1, face->Element2, phi2, penalizationCoefficient);

					//cout << "\t\t\t c=" << coupling << "\tp=" << penalization << endl;
					
					if ((action & Action::ExtractComponentMatrices) == Action::ExtractComponentMatrices)
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

		this->A = Eigen::SparseMatrix<double>(nUnknowns, nUnknowns);
		matrixCoeffs.Fill(this->A);
		cout << "nnz(A) = " << this->A.nonZeros() << endl;

		if ((action & Action::ExtractSystem) == Action::ExtractSystem)
		{
			cout << "Export..." << endl;
			Eigen::saveMarket(this->A, matrixFilePath);
			cout << "Matrix exported to \t" << matrixFilePath << endl;

			Eigen::saveMarketVector(this->b, rhsFilePath);
			cout << "RHS exported to \t" << rhsFilePath << endl;
		}

		if ((action & Action::ExtractMassMatrix) == Action::ExtractMassMatrix)
		{
			Eigen::SparseMatrix<double> M(nUnknowns, nUnknowns);
			massMatrixCoeffs.Fill(M);
			Eigen::saveMarket(M, massMatrixFilePath);
			cout << "Mass matrix exported to \t" << massMatrixFilePath << endl;
		}

		if ((action & Action::ExtractComponentMatrices) == Action::ExtractComponentMatrices)
		{
			Eigen::SparseMatrix<double> V(nUnknowns, nUnknowns);
			volumicCoeffs.Fill(V);
			Eigen::saveMarket(V, matrixVolumicFilePath);
			cout << "Volumic part exported to \t" << matrixVolumicFilePath << endl;

			Eigen::SparseMatrix<double> C(nUnknowns, nUnknowns);
			couplingCoeffs.Fill(C);
			Eigen::saveMarket(C, matrixCouplingFilePath);
			cout << "Coupling part exported to \t" << matrixCouplingFilePath << endl;

			Eigen::SparseMatrix<double> P(nUnknowns, nUnknowns);
			penCoeffs.Fill(P);
			Eigen::saveMarket(P, matrixPenFilePath);
			cout << "Penalization part exported to \t" << matrixPenFilePath << endl;
		}

	}
};

