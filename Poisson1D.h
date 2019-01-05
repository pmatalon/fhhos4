#pragma once
#include <functional>
#include "FunctionalBasisWithNumbers.h"
#include "IPoisson1D_DGTerms.h"
#include "CartesianGrid1D.h"
#include "IBasisFunction1D.h"
#include "Element.h"
#include "FileMatrix.h"
#include "FileVector.h"
#include <math.h>
using namespace std;

class Poisson1D
{
private:
	string _solution;
	std::function<double(double)> _sourceFunction;
public:
	Poisson1D(string solution, function<double(double)> sourceFunction)
	{
		cout << "----------------------------------------" << endl;
		cout << "----------------------------------------" << endl;
		this->_solution = solution;
		this->_sourceFunction = sourceFunction;
	}

	void DiscretizeDG(CartesianGrid1D* grid, FunctionalBasisWithNumbers* basis, IPoisson1D_DGTerms* dg, int penalizationCoefficient, string outputDirectory, bool extractMatrixComponents)
	{
		bool autoPenalization = penalizationCoefficient == -1;
		if (autoPenalization)
			penalizationCoefficient = pow(basis->GetDegree() + 1, 2) * grid->NElements(); // Ralph-Hartmann

		cout << "Discretization: Discontinuous Galerkin SIPG" << endl;
		cout << "\tPenalization coefficient: " << penalizationCoefficient << endl;
		cout << "\tBasis of polynomials: " << (dg->IsGlobalBasis() ? "global" : "") + basis->Name() << endl;

		cout << "Local functions: " << basis->NumberOfLocalFunctionsInElement(0) << endl;
		for (int localFunctionNumber = 0; localFunctionNumber < basis->NumberOfLocalFunctionsInElement(0); localFunctionNumber++)
		{
			IBasisFunction1D* localFunction = basis->GetLocalBasisFunction(0, localFunctionNumber);
			cout << "\t " << localFunction->ToString() << endl;
		}
		int nUnknowns = grid->NElements() * basis->NumberOfLocalFunctionsInElement(0);
		//cout << "Unknowns: " << nUnknowns << endl;

		string fileName = "Poisson1D" + this->_solution + "_n" + to_string(grid->NElements()) + "_DG_SIPG_" + (dg->IsGlobalBasis() ? "global" : "") + basis->Name() + "_pen" + (autoPenalization ? "-1" : to_string(penalizationCoefficient));
		string matrixFilePath = outputDirectory + "/" + fileName + "_A.dat";
		string matrixVolumicFilePath = outputDirectory + "/" + fileName + "_A_volumic.dat";
		string matrixCouplingFilePath = outputDirectory + "/" + fileName + "_A_coupling.dat";
		string matrixPenFilePath = outputDirectory + "/" + fileName + "_A_pen.dat";
		string massMatrixFilePath = outputDirectory + "/" + fileName + "_Mass.dat";
		FileMatrix* fileMatrix = new FileMatrix(nUnknowns, nUnknowns, matrixFilePath);
		FileMatrix* fileMatrixVolumic = NULL;
		FileMatrix* fileMatrixCoupling = NULL;
		FileMatrix* fileMatrixPen = NULL;
		FileMatrix* fileMassMatrix = NULL;
		if (extractMatrixComponents)
		{
			fileMatrixVolumic = new FileMatrix(nUnknowns, nUnknowns, matrixVolumicFilePath);
			fileMatrixCoupling = new FileMatrix(nUnknowns, nUnknowns, matrixCouplingFilePath);
			fileMatrixPen = new FileMatrix(nUnknowns, nUnknowns, matrixPenFilePath);
			fileMassMatrix = new FileMatrix(nUnknowns, nUnknowns, massMatrixFilePath);
		}

		string rhsFilePath = outputDirectory + "/" + fileName + "_b.dat";
		FileVector* fileRHS = new FileVector(rhsFilePath);

		for (int element = 0; element < grid->NElements(); element++) // interval [element/n, (element+1)/n]
		{
			for (int localFunctionNumber1 = 0; localFunctionNumber1 < basis->NumberOfLocalFunctionsInElement(element); localFunctionNumber1++)
			{
				IBasisFunction1D* localFunction1 = basis->GetLocalBasisFunction(element, localFunctionNumber1);
				BigNumber basisFunction1 = basis->GlobalFunctionNumber(element, localFunctionNumber1);


				if (!grid->IsFirstElement(element))
				{
					// Left neighbour
					int neighbour = element - 1;
					int interface = grid->GetInterface(element, neighbour);

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
							fileMatrixCoupling->Add(basisFunction1, basisFunction2, couplingTerm);
							fileMatrixPen->Add(basisFunction1, basisFunction2, penalization);
						}
						fileMatrix->Add(basisFunction1, basisFunction2, couplingTerm + penalization);
					}
				}

				// Current element (block diagonal)
				//cout << "Element " << element << endl;
				for (int localFunctionNumber2 = 0; localFunctionNumber2 < basis->NumberOfLocalFunctionsInElement(element); localFunctionNumber2++)
				{
					IBasisFunction1D* localFunction2 = basis->GetLocalBasisFunction(element, localFunctionNumber2);
					int basisFunction2 = basis->GlobalFunctionNumber(element, localFunctionNumber2);

					double volumicTerm = dg->VolumicTerm(element, localFunction1, localFunction2);

					double couplingTermLeft = dg->CouplingTerm(grid->LeftInterface(element), element, localFunction1, element, localFunction2);
					double penalizationLeft = dg->PenalizationTerm(grid->LeftInterface(element), element, localFunction1, element, localFunction2, penalizationCoefficient);

					double couplingTermRight = dg->CouplingTerm(grid->RightInterface(element), element, localFunction1, element, localFunction2);
					double penalizationRight = dg->PenalizationTerm(grid->RightInterface(element), element, localFunction1, element, localFunction2, penalizationCoefficient);

					//cout << "\t func1 = " << localFunction1->ToString() << " func2 = " << localFunction2->ToString() << endl;
					//cout << "\t\t Volumic = " << volumicTerm << endl;
					//cout << "\t\t Interface left:\t c=" << couplingTermLeft << "\tp=" << penalizationLeft << endl;
					//cout << "\t\t Interface right:\t c=" << couplingTermRight << "\tp=" << penalizationRight << endl;

					if (extractMatrixComponents)
					{
						fileMatrixVolumic->Add(basisFunction1, basisFunction2, volumicTerm);
						fileMatrixCoupling->Add(basisFunction1, basisFunction2, couplingTermLeft + couplingTermRight);
						fileMatrixPen->Add(basisFunction1, basisFunction2, penalizationRight + penalizationLeft);

						double massTerm = dg->MassTerm(element, localFunction1, localFunction2);
						fileMassMatrix->Add(basisFunction1, basisFunction2, massTerm);
					}
					fileMatrix->Add(basisFunction1, basisFunction2, volumicTerm + couplingTermLeft + couplingTermRight + penalizationRight + penalizationLeft);


					//cout << "v = " << volumicTerm << " couplingTermLeft=" << couplingTermLeft << " couplingTermRight=" << couplingTermRight << " penalizationRight=" << penalizationRight << " penalizationRight=" << penalizationRight;
				}

				if (!grid->IsLastElement(element))
				{
					// Right neighbour
					int neighbour = element + 1;
					int interface = grid->GetInterface(element, neighbour);

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
							fileMatrixCoupling->Add(basisFunction1, basisFunction2, couplingTerm);
							fileMatrixPen->Add(basisFunction1, basisFunction2, penalization);
						}
						fileMatrix->Add(basisFunction1, basisFunction2, couplingTerm + penalization);
					}
				}

				double rhs = dg->RightHandSide(element, localFunction1);
				fileRHS->Add(rhs);
			}
		}

		if (extractMatrixComponents)
		{
			delete fileMatrixVolumic;
			delete fileMatrixCoupling;
			delete fileMatrixPen;
			delete fileMassMatrix;
		}
		delete fileMatrix;
		delete fileRHS;

		cout << "Matrix exported to " << matrixFilePath << endl;
		cout << "RHS exported to " << rhsFilePath << endl;
	}

	~Poisson1D()
	{
	}
};

