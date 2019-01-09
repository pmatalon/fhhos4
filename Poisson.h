#pragma once
#include <iostream>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include "Problem.h"
#include "IMesh.h"
#include "FunctionalBasisWithObjects.h"
#include "ElementInterface.h"
#include "IPoisson_DGTerms.h"
#include "NonZeroCoefficients.h"
#include "L2.h"
using namespace std;

template <class IBasisFunction>
class Poisson : public Problem
{
private:

public:

	Poisson(string solutionName) : Problem(solutionName)
	{	}

	void DiscretizeDG(IMesh* mesh, FunctionalBasisWithObjects<IBasisFunction>* basis, IPoisson_DGTerms<IBasisFunction>* dg, int penalizationCoefficient, string outputDirectory, bool extractMatrixComponents)
	{
		bool autoPenalization = penalizationCoefficient == -1;
		if (autoPenalization)
			penalizationCoefficient = pow(mesh->Dim, 2) * pow(basis->GetDegree() + 1, 2) * mesh->N; // Ralph-Hartmann

		cout << "Discretization: Discontinuous Galerkin SIPG" << endl;
		cout << "\tPenalization coefficient: " << penalizationCoefficient << endl;
		cout << "\tBasis of polynomials: " << (dg->IsGlobalBasis() ? "global" : "") + basis->Name() << endl;
		
		cout << "Local functions: " << basis->NumberOfLocalFunctionsInElement(0) << endl;
		for (int localFunctionNumber = 0; localFunctionNumber < basis->NumberOfLocalFunctionsInElement(NULL); localFunctionNumber++)
		{
			IBasisFunction* localFunction = basis->GetLocalBasisFunction(NULL, localFunctionNumber);
			cout << "\t " << localFunction->ToString() << endl;
		}
		BigNumber nUnknowns = static_cast<int>(mesh->Elements.size()) * basis->NumberOfLocalFunctionsInElement(NULL);
		cout << "Unknowns: " << nUnknowns << endl;

		string fileName = "Poisson" + to_string(mesh->Dim) + "D" + this->_solutionName + "_n" + to_string(mesh->N) + "_DG_SIPG_" + (dg->IsGlobalBasis() ? "global" : "") + basis->Name() + "_pen" + (autoPenalization ? "-1" : to_string(penalizationCoefficient));
		string matrixFilePath			= outputDirectory + "/" + fileName + "_A.dat";
		string matrixVolumicFilePath	= outputDirectory + "/" + fileName + "_A_volumic.dat";
		string matrixCouplingFilePath	= outputDirectory + "/" + fileName + "_A_coupling.dat";
		string matrixPenFilePath		= outputDirectory + "/" + fileName + "_A_pen.dat";
		string rhsFilePath				= outputDirectory + "/" + fileName + "_b.dat";

		this->b = Eigen::VectorXd(nUnknowns);

		BigNumber nnzApproximate = extractMatrixComponents ? mesh->Elements.size() * basis->NumberOfLocalFunctionsInElement(NULL) * (2 * mesh->Dim + 1) : 0;
		NonZeroCoefficients matrixCoeffs(nnzApproximate);
		NonZeroCoefficients volumicCoeffs(nnzApproximate);
		NonZeroCoefficients couplingCoeffs(nnzApproximate);
		NonZeroCoefficients penCoeffs(nnzApproximate);

		//--------------------------------------------//
		// Iteration on the elements: diagonal blocks //
		//--------------------------------------------//

		for (Element* element : mesh->Elements)
		{
			//cout << "Element " << element->Number << endl;
			vector<ElementInterface*> elementInterfaces = element->Interfaces;

			for (int localFunctionNumber1 = 0; localFunctionNumber1 < basis->NumberOfLocalFunctionsInElement(element); localFunctionNumber1++)
			{
				IBasisFunction* localFunction1 = basis->GetLocalBasisFunction(element, localFunctionNumber1);
				BigNumber basisFunction1 = basis->GlobalFunctionNumber(element, localFunctionNumber1);

				// Current element (block diagonal)
				for (int localFunctionNumber2 = 0; localFunctionNumber2 < basis->NumberOfLocalFunctionsInElement(element); localFunctionNumber2++)
				{
					IBasisFunction* localFunction2 = basis->GetLocalBasisFunction(element, localFunctionNumber2);
					BigNumber basisFunction2 = basis->GlobalFunctionNumber(element, localFunctionNumber2);

					double volumicTerm = dg->VolumicTerm(element, localFunction1, localFunction2);
					
					double coupling = 0;
					double penalization = 0;
					for (ElementInterface* elemInterface : elementInterfaces)
					{
						coupling += dg->CouplingTerm(elemInterface, element, localFunction1, element, localFunction2);
						penalization += dg->PenalizationTerm(elemInterface, element, localFunction1, element, localFunction2, penalizationCoefficient);
						//cout << "\t\t " << elemInterface->ToString() << ":\t c=" << coupling << "\tp=" << penalization << endl;
					}

					//cout << "\t\t TOTAL = " << volumicTerm + coupling + penalization << endl;
					
					if (extractMatrixComponents)
					{
						volumicCoeffs.Add(basisFunction1, basisFunction2, volumicTerm);
						couplingCoeffs.Add(basisFunction1, basisFunction2, coupling);
						penCoeffs.Add(basisFunction1, basisFunction2, penalization);
					}
					matrixCoeffs.Add(basisFunction1, basisFunction2, volumicTerm + coupling + penalization);
				}

				double rhs = dg->RightHandSide(element, localFunction1);
				this->b(basisFunction1) = rhs;
			}
		}

		//--------------------------------------------------//
		// Iteration on the interfaces: off-diagonal blocks //
		//--------------------------------------------------//

		for (ElementInterface* interface : mesh->Interfaces)
		{
			if (interface->IsDomainBoundary)
				continue;

			for (int localFunctionNumber1 = 0; localFunctionNumber1 < basis->NumberOfLocalFunctionsInElement(interface->Element1); localFunctionNumber1++)
			{
				IBasisFunction* localFunction1 = basis->GetLocalBasisFunction(interface->Element1, localFunctionNumber1);
				BigNumber basisFunction1 = basis->GlobalFunctionNumber(interface->Element1, localFunctionNumber1);
				for (int localFunctionNumber2 = 0; localFunctionNumber2 < basis->NumberOfLocalFunctionsInElement(interface->Element2); localFunctionNumber2++)
				{
					IBasisFunction* localFunction2 = basis->GetLocalBasisFunction(interface->Element2, localFunctionNumber2);
					BigNumber basisFunction2 = basis->GlobalFunctionNumber(interface->Element2, localFunctionNumber2);
					double coupling = dg->CouplingTerm(interface, interface->Element1, localFunction1, interface->Element2, localFunction2);
					double penalization = dg->PenalizationTerm(interface, interface->Element1, localFunction1, interface->Element2, localFunction2, penalizationCoefficient);
					
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

