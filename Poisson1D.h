#pragma once
#include <functional>
#include "FunctionalBasisWithNumbers.h"
#include "CartesianGrid1D.h"
#include "IBasisFunction1D.h"
#include "Element.h"
#include "FileMatrix.h"
#include "FileVector.h"
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

	void DiscretizeDG(CartesianGrid1D* grid, FunctionalBasisWithNumbers* basis, int penalizationCoefficient, string outputDirectory)
	{
		cout << "Discretization: Discontinuous Galerkin SIPG" << endl;
		cout << "\tPenalization coefficient: " << penalizationCoefficient << endl;
		cout << "\tBasis of polynomials: " << basis->Name() << endl;

		int nUnknowns = grid->NElements() * basis->NumberOfLocalFunctionsInElement(0);
		//cout << "Unknowns: " << nUnknowns << endl;

		//string terms = "volumic";
		//string terms = "coupling";
		//string terms = "penalization";
		string terms = "";

		string fileName = "Poisson1D" + this->_solution + "_n" + to_string(grid->NElements()) + "_DG_SIPG_" + basis->Name() + "_pen" + to_string(penalizationCoefficient);
		string matrixFilePath = outputDirectory + "/" + fileName + "_A" + terms + ".dat";
		FileMatrix* fileMatrix = new FileMatrix(nUnknowns, nUnknowns, matrixFilePath);

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
					for (int localFunctionNumber2 = 0; localFunctionNumber2 < basis->NumberOfLocalFunctionsInElement(neighbour); localFunctionNumber2++)
					{
						IBasisFunction1D* localFunction2 = basis->GetLocalBasisFunction(element, localFunctionNumber2);

						double couplingTerm = basis->CouplingTerm(interface, element, localFunction1, neighbour, localFunction2);
						double penalization = basis->PenalizationTerm(interface, element, localFunction1, neighbour, localFunction2);

						BigNumber basisFunction2 = basis->GlobalFunctionNumber(neighbour, localFunctionNumber2);

						if (terms.compare("volumic") == 0) {}
						else if (terms.compare("coupling") == 0)
							fileMatrix->Add(basisFunction1, basisFunction2, couplingTerm);
						else if (terms.compare("penalization") == 0)
							fileMatrix->Add(basisFunction1, basisFunction2, penalization);
						else
							fileMatrix->Add(basisFunction1, basisFunction2, couplingTerm + penalization);

					}
				}

				// Current element (block diagonal)
				for (int localFunctionNumber2 = 0; localFunctionNumber2 < basis->NumberOfLocalFunctionsInElement(element); localFunctionNumber2++)
				{
					IBasisFunction1D* localFunction2 = basis->GetLocalBasisFunction(element, localFunctionNumber2);

					double volumicTerm = basis->VolumicTerm(element, localFunction1, localFunction2);

					double couplingTermLeft = basis->CouplingTerm(grid->LeftInterface(element), element, localFunction1, element, localFunction2);
					double penalizationLeft = basis->PenalizationTerm(grid->LeftInterface(element), element, localFunction1, element, localFunction2);

					double couplingTermRight = basis->CouplingTerm(grid->RightInterface(element), element, localFunction1, element, localFunction2);
					double penalizationRight = basis->PenalizationTerm(grid->RightInterface(element), element, localFunction1, element, localFunction2);

					int basisFunction2 = basis->GlobalFunctionNumber(element, localFunctionNumber2);

					if (terms.compare("volumic") == 0)
						fileMatrix->Add(basisFunction1, basisFunction2, volumicTerm);
					else if (terms.compare("coupling") == 0)
						fileMatrix->Add(basisFunction1, basisFunction2, couplingTermLeft + couplingTermRight);
					else if (terms.compare("penalization") == 0)
						fileMatrix->Add(basisFunction1, basisFunction2, penalizationRight + penalizationLeft);
					else
						fileMatrix->Add(basisFunction1, basisFunction2, volumicTerm + couplingTermLeft + couplingTermRight + penalizationRight + penalizationLeft);


					cout << "v = " << volumicTerm << " couplingTermLeft=" << couplingTermLeft << " couplingTermRight=" << couplingTermRight << " penalizationRight=" << penalizationRight << " penalizationRight=" << penalizationRight;
				}

				if (!grid->IsLastElement(element))
				{
					// Right neighbour
					int neighbour = element + 1;
					int interface = grid->GetInterface(element, neighbour);
					for (int localFunctionNumber2 = 0; localFunctionNumber2 < basis->NumberOfLocalFunctionsInElement(neighbour); localFunctionNumber2++)
					{
						IBasisFunction1D* localFunction2 = basis->GetLocalBasisFunction(element, localFunctionNumber2);

						double couplingTerm = basis->CouplingTerm(interface, element, localFunction1, neighbour, localFunction2);
						double penalization = basis->PenalizationTerm(interface, element, localFunction1, neighbour, localFunction2);

						int basisFunction2 = basis->GlobalFunctionNumber(neighbour, localFunctionNumber2);

						if (terms.compare("volumic") == 0) {}
						else if (terms.compare("coupling") == 0)
							fileMatrix->Add(basisFunction1, basisFunction2, couplingTerm);
						else if (terms.compare("penalization") == 0)
							fileMatrix->Add(basisFunction1, basisFunction2, penalization);
						else
							fileMatrix->Add(basisFunction1, basisFunction2, couplingTerm + penalization);
					}
				}

				double rhs = basis->RightHandSide(element, localFunction1);
				fileRHS->Add(rhs);
			}
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

