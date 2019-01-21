#pragma once
#include <iostream>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include "Problem.h"
#include "IMesh.h"
#include "FunctionalBasisWithObjects.h"
#include "Face.h"
#include "IPoisson_DGTerms.h"
#include "NonZeroCoefficients.h"
#include "L2.h"
using namespace std;

class Poisson : public Problem
{
private:

public:

	Poisson(string solutionName) : Problem(solutionName)
	{	}

	void DiscretizeDG(IMesh* mesh, FunctionalBasisWithObjects* basis, IPoisson_DGTerms* dg, int penalizationCoefficient, string outputDirectory, bool extractMatrixComponents)
	{
		bool autoPenalization = penalizationCoefficient == -1;
		if (autoPenalization)
			penalizationCoefficient = pow(mesh->Dim, 2) * pow(basis->GetDegree() + 1, 2) * mesh->N; // Ralph-Hartmann

		cout << "Discretization: Discontinuous Galerkin SIPG" << endl;
		cout << "\tPenalization coefficient: " << penalizationCoefficient << endl;
		cout << "\tBasis of polynomials: " << (dg->IsGlobalBasis() ? "global" : "") + basis->Name() << endl;
		
		cout << "Local functions: " << basis->NumberOfLocalFunctionsInElement(NULL) << endl;
		for (BasisFunction* phi : basis->LocalFunctions)
			cout << "\t " << phi->ToString() << endl;
		BigNumber nUnknowns = static_cast<int>(mesh->Elements.size()) * basis->NumberOfLocalFunctionsInElement(NULL);
		cout << "Unknowns: " << nUnknowns << endl;

		string fileName = "Poisson" + to_string(mesh->Dim) + "D" + this->_solutionName + "_n" + to_string(mesh->N) + "_DG_SIPG_" + (dg->IsGlobalBasis() ? "global" : "") + basis->Name() + "_pen" + (autoPenalization ? "-1" : to_string(penalizationCoefficient));
		string matrixFilePath			= outputDirectory + "/" + fileName + "_A.dat";
		string matrixVolumicFilePath	= outputDirectory + "/" + fileName + "_A_volumic.dat";
		string matrixCouplingFilePath	= outputDirectory + "/" + fileName + "_A_coupling.dat";
		string matrixPenFilePath		= outputDirectory + "/" + fileName + "_A_pen.dat";
		string massMatrixFilePath		= outputDirectory + "/" + fileName + "_Mass.dat";
		string rhsFilePath				= outputDirectory + "/" + fileName + "_b.dat";

		this->b = Eigen::VectorXd(nUnknowns);

		BigNumber nnzApproximate = extractMatrixComponents ? mesh->Elements.size() * basis->NumberOfLocalFunctionsInElement(NULL) * (2 * mesh->Dim + 1) : 0;
		NonZeroCoefficients matrixCoeffs(nnzApproximate);
		NonZeroCoefficients massMatrixCoeffs(nnzApproximate);
		NonZeroCoefficients volumicCoeffs(nnzApproximate);
		NonZeroCoefficients couplingCoeffs(nnzApproximate);
		NonZeroCoefficients penCoeffs(nnzApproximate);

		//--------------------------------------------//
		// Iteration on the elements: diagonal blocks //
		//--------------------------------------------//

		for (auto element : mesh->Elements)
		{
			//cout << "Element " << element->Number << endl;

			for (BasisFunction* phi1 : basis->LocalFunctions)
			{
				BigNumber basisFunction1 = basis->GlobalFunctionNumber(element, phi1);

				// Current element (block diagonal)
				for (BasisFunction* phi2 : basis->LocalFunctions)
				{
					BigNumber basisFunction2 = basis->GlobalFunctionNumber(element, phi2);

					//cout << "\t phi" << phi1->LocalNumber << " = " << phi1->ToString() << " phi" << phi2->LocalNumber << " = " << phi2->ToString() << endl;

					double volumicTerm = dg->VolumicTerm(element, phi1, phi2);
					//cout << "\t\t volumic = " << volumicTerm << endl;
					double massTerm = dg->MassTerm(element, phi1, phi2);
					
					double coupling = 0;
					double penalization = 0;
					for (Face* face : element->Faces)
					{
						double c = dg->CouplingTerm(face, element, phi1, element, phi2);
						double p = dg->PenalizationTerm(face, element, phi1, element, phi2, penalizationCoefficient);
						coupling += c;
						penalization += p;
						//cout << "\t\t " << face->ToString() << ":\t c=" << c << "\tp=" << p << endl;
					}

					//cout << "\t\t TOTAL = " << volumicTerm + coupling + penalization << endl;
					
					if (extractMatrixComponents)
					{
						volumicCoeffs.Add(basisFunction1, basisFunction2, volumicTerm);
						couplingCoeffs.Add(basisFunction1, basisFunction2, coupling);
						penCoeffs.Add(basisFunction1, basisFunction2, penalization);
					}
					matrixCoeffs.Add(basisFunction1, basisFunction2, volumicTerm + coupling + penalization);
					massMatrixCoeffs.Add(basisFunction1, basisFunction2, massTerm);
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

			for (BasisFunction* phi1 : basis->LocalFunctions)
			{
				BigNumber basisFunction1 = basis->GlobalFunctionNumber(face->Element1, phi1);
				for (BasisFunction* phi2 : basis->LocalFunctions)
				{
					BigNumber basisFunction2 = basis->GlobalFunctionNumber(face->Element2, phi2);
					double coupling = dg->CouplingTerm(face, face->Element1, phi1, face->Element2, phi2);
					double penalization = dg->PenalizationTerm(face, face->Element1, phi1, face->Element2, phi2, penalizationCoefficient);
					
					if (extractMatrixComponents)
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
		Eigen::saveMarket(this->A, matrixFilePath);

		Eigen::SparseMatrix<double> M(nUnknowns, nUnknowns);
		massMatrixCoeffs.Fill(M);
		Eigen::saveMarket(M, massMatrixFilePath);

		Eigen::saveMarketVector(this->b, rhsFilePath);

		if (extractMatrixComponents)
		{
			Eigen::SparseMatrix<double> V(nUnknowns, nUnknowns);
			volumicCoeffs.Fill(V);
			Eigen::saveMarket(V, matrixVolumicFilePath);

			Eigen::SparseMatrix<double> C(nUnknowns, nUnknowns);
			couplingCoeffs.Fill(C);
			Eigen::saveMarket(C, matrixCouplingFilePath);

			Eigen::SparseMatrix<double> P(nUnknowns, nUnknowns);
			penCoeffs.Fill(P);
			Eigen::saveMarket(P, matrixPenFilePath);
		}

		cout << "Matrix exported to \t" << matrixFilePath << endl;
		cout << "RHS exported to \t" << rhsFilePath << endl;
	}
};
