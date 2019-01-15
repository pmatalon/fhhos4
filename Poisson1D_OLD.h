#pragma once
#include <functional>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include "Problem.h"
#include "FunctionalBasisWithNumbers.h"
#include "IPoisson1D_DGTerms.h"
#include "CartesianGrid1DOLD.h"
#include "Element.h"
#include <math.h>
#include "NonZeroCoefficients.h"
using namespace std;

class Poisson1D : public Problem
{
private:
	std::function<double(double)> _sourceFunction;
public:
	Poisson1D(string solutionName, function<double(double)> sourceFunction) : Problem(solutionName)
	{
		this->_sourceFunction = sourceFunction;
	}

	void DiscretizeDG(CartesianGrid1DOLD* mesh, FunctionalBasisWithNumbers* basis, IPoisson1D_DGTerms* dg, int penalizationCoefficient, string outputDirectory, bool extractMatrixComponents)
	{
		bool autoPenalization = penalizationCoefficient == -1;
		if (autoPenalization)
			penalizationCoefficient = pow(basis->GetDegree() + 1, 2) * mesh->N; // Ralph-Hartmann

		cout << "Discretization: Discontinuous Galerkin SIPG" << endl;
		cout << "\tPenalization coefficient: " << penalizationCoefficient << endl;
		cout << "\tBasis of polynomials: " << (dg->IsGlobalBasis() ? "global" : "") + basis->Name() << endl;

		cout << "Local functions: " << basis->NumberOfLocalFunctionsInElement(0) << endl;
		for (int localFunctionNumber = 0; localFunctionNumber < basis->NumberOfLocalFunctionsInElement(0); localFunctionNumber++)
		{
			IBasisFunction1D* localFunction = basis->GetLocalBasisFunction(0, localFunctionNumber);
			cout << "\t " << localFunction->ToString() << endl;
		}
		int nUnknowns = mesh->NElements() * basis->NumberOfLocalFunctionsInElement(0);
		cout << "Unknowns: " << nUnknowns << endl;

		string fileName = "Poisson1D" + this->_solutionName + "_n" + to_string(mesh->NElements()) + "_DG_SIPG_" + (dg->IsGlobalBasis() ? "global" : "") + basis->Name() + "_pen" + (autoPenalization ? "-1" : to_string(penalizationCoefficient));
		string matrixFilePath			= outputDirectory + "/" + fileName + "_A.dat";
		string matrixVolumicFilePath	= outputDirectory + "/" + fileName + "_A_volumic.dat";
		string matrixCouplingFilePath	= outputDirectory + "/" + fileName + "_A_coupling.dat";
		string matrixPenFilePath		= outputDirectory + "/" + fileName + "_A_pen.dat";
		string massMatrixFilePath		= outputDirectory + "/" + fileName + "_Mass.dat";
		string rhsFilePath				= outputDirectory + "/" + fileName + "_b.dat";

		BigNumber nnzApproximate = extractMatrixComponents ? mesh->NElements() * basis->NumberOfLocalFunctionsInElement(0) * 3 : 0;
		NonZeroCoefficients matrixCoeffs(nnzApproximate);
		NonZeroCoefficients massMatrixCoeffs(nnzApproximate);
		NonZeroCoefficients volumicCoeffs(nnzApproximate);
		NonZeroCoefficients couplingCoeffs(nnzApproximate);
		NonZeroCoefficients penCoeffs(nnzApproximate);

		this->b = Eigen::VectorXd(nUnknowns);

		for (int element = 0; element < mesh->NElements(); element++) // interval [element/n, (element+1)/n]
		{
			for (int localFunctionNumber1 = 0; localFunctionNumber1 < basis->NumberOfLocalFunctionsInElement(element); localFunctionNumber1++)
			{
				IBasisFunction1D* localFunction1 = basis->GetLocalBasisFunction(element, localFunctionNumber1);
				BigNumber basisFunction1 = basis->GlobalFunctionNumber(element, localFunctionNumber1);

				if (!mesh->IsFirstElement(element))
				{
					// Left neighbour
					int neighbour = element - 1;
					int interface = mesh->GetInterface(element, neighbour);

					//cout << "Interface " << interface << " between element " << element << " and " << neighbour << endl;

					for (int localFunctionNumber2 = 0; localFunctionNumber2 < basis->NumberOfLocalFunctionsInElement(neighbour); localFunctionNumber2++)
					{
						IBasisFunction1D* localFunction2 = basis->GetLocalBasisFunction(element, localFunctionNumber2);
						BigNumber basisFunction2 = basis->GlobalFunctionNumber(neighbour, localFunctionNumber2);
						
						double couplingTerm = dg->CouplingTerm(interface, element, localFunction1, neighbour, localFunction2);
						double penalization = dg->PenalizationTerm(interface, element, localFunction1, neighbour, localFunction2, penalizationCoefficient);

						//cout << "\t func1 = " << localFunction1->ToString() << " func2 = " << localFunction2->ToString() << endl;
						//cout << "\t\t Interface " << interface << ":\t c=" << couplingTerm << "\tp=" << penalization << endl;

						if (extractMatrixComponents)
						{
							couplingCoeffs.Add(basisFunction1, basisFunction2, couplingTerm);
							penCoeffs.Add(basisFunction1, basisFunction2, penalization);
						}
						matrixCoeffs.Add(basisFunction1, basisFunction2, couplingTerm + penalization);
					}
				}

				// Current element (block diagonal)
				//cout << "Element " << element << endl;
				for (int localFunctionNumber2 = 0; localFunctionNumber2 < basis->NumberOfLocalFunctionsInElement(element); localFunctionNumber2++)
				{
					IBasisFunction1D* localFunction2 = basis->GetLocalBasisFunction(element, localFunctionNumber2);
					int basisFunction2 = basis->GlobalFunctionNumber(element, localFunctionNumber2);

					double volumicTerm = dg->VolumicTerm(element, localFunction1, localFunction2);

					double couplingTermLeft = dg->CouplingTerm(mesh->LeftFace(element), element, localFunction1, element, localFunction2);
					double penalizationLeft = dg->PenalizationTerm(mesh->LeftFace(element), element, localFunction1, element, localFunction2, penalizationCoefficient);

					double couplingTermRight = dg->CouplingTerm(mesh->RightFace(element), element, localFunction1, element, localFunction2);
					double penalizationRight = dg->PenalizationTerm(mesh->RightFace(element), element, localFunction1, element, localFunction2, penalizationCoefficient);

					//cout << "\t func1 = " << localFunction1->ToString() << " func2 = " << localFunction2->ToString() << endl;
					//cout << "\t\t Volumic = " << volumicTerm << endl;
					//cout << "\t\t Interface left:\t c=" << couplingTermLeft << "\tp=" << penalizationLeft << endl;
					//cout << "\t\t Interface right:\t c=" << couplingTermRight << "\tp=" << penalizationRight << endl;

					if (extractMatrixComponents)
					{
						volumicCoeffs.Add(basisFunction1, basisFunction2, volumicTerm);
						couplingCoeffs.Add(basisFunction1, basisFunction2, couplingTermLeft + couplingTermRight);
						penCoeffs.Add(basisFunction1, basisFunction2, penalizationRight + penalizationLeft);

						double massTerm = dg->MassTerm(element, localFunction1, localFunction2);
						massMatrixCoeffs.Add(basisFunction1, basisFunction2, massTerm);
					}
					matrixCoeffs.Add(basisFunction1, basisFunction2, volumicTerm + couplingTermLeft + couplingTermRight + penalizationRight + penalizationLeft);

					//cout << "v = " << volumicTerm << " couplingTermLeft=" << couplingTermLeft << " couplingTermRight=" << couplingTermRight << " penalizationRight=" << penalizationRight << " penalizationRight=" << penalizationRight;
				}

				if (!mesh->IsLastElement(element))
				{
					// Right neighbour
					int neighbour = element + 1;
					int interface = mesh->GetInterface(element, neighbour);

					//cout << "Interface " << interface << " between element " << element << " and " << neighbour << endl;

					for (int localFunctionNumber2 = 0; localFunctionNumber2 < basis->NumberOfLocalFunctionsInElement(neighbour); localFunctionNumber2++)
					{
						IBasisFunction1D* localFunction2 = basis->GetLocalBasisFunction(element, localFunctionNumber2);
						int basisFunction2 = basis->GlobalFunctionNumber(neighbour, localFunctionNumber2);

						double couplingTerm = dg->CouplingTerm(interface, element, localFunction1, neighbour, localFunction2);
						double penalization = dg->PenalizationTerm(interface, element, localFunction1, neighbour, localFunction2, penalizationCoefficient);

						//cout << "\t func1 = " << localFunction1->ToString() << " func2 = " << localFunction2->ToString() << endl;
						//cout << "\t\t Interface " << interface << ":\t c=" << couplingTerm << "\tp=" << penalization << endl;

						if (extractMatrixComponents)
						{
							couplingCoeffs.Add(basisFunction1, basisFunction2, couplingTerm);
							penCoeffs.Add(basisFunction1, basisFunction2, penalization);
						}
						matrixCoeffs.Add(basisFunction1, basisFunction2, couplingTerm + penalization);
					}
				}

				double rhs = dg->RightHandSide(element, localFunction1);
				this->b(basisFunction1) = rhs;
			}
		}

		this->A = Eigen::SparseMatrix<double>(nUnknowns, nUnknowns);
		matrixCoeffs.Fill(this->A);
		Eigen::saveMarket(this->A, matrixFilePath);
		Eigen::saveMarketVector(this->b, rhsFilePath);

		if (extractMatrixComponents)
		{
			Eigen::SparseMatrix<double> M(nUnknowns, nUnknowns);
			massMatrixCoeffs.Fill(M);
			Eigen::saveMarket(M, massMatrixFilePath);

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

		cout << "Matrix exported to " << matrixFilePath << endl;
		cout << "RHS exported to " << rhsFilePath << endl;
	}

	~Poisson1D()
	{
	}
};

