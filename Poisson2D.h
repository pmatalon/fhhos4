#include <iostream>
#include <functional>
#include "CartesianGrid2D.h"
#include "FunctionalBasisWithObjects.h"
#include "ElementInterface.h"
#include "FileMatrix.h"
#include "FileVector.h"
using namespace std;

#pragma once
class Poisson2D
{
private:
	string _solution;
	function<double(double, double)> _sourceFunction;
public:
	Poisson2D(string solution, function<double(double, double)> sourceFunction)
	{
		this->_solution = solution;
		this->_sourceFunction = sourceFunction;
	}

	void DiscretizeDG(CartesianGrid2D* grid, FunctionalBasisWithObjects* basis, int penalizationCoefficient, string outputDirectory)
	{
		cout << "Discretization: Discontinuous Galerkin SIPG" << endl;
		cout << "\tPenalization coefficient: " << penalizationCoefficient << endl;
		cout << "\tBasis of polynomials: " << basis->Name() << endl;
		
		BigNumber nUnknowns = static_cast<int>(grid->Elements.size()) * basis->NumberOfLocalFunctionsInElement(0);
		cout << "Unknowns: " << nUnknowns << endl;

		//string terms = "volumic";
		//string terms = "coupling";
		//string terms = "penalization";
		string terms = "";

		string fileName = "Poisson2D" + this->_solution + "_n" + to_string(grid->N) + "_DG_SIPG_" + basis->Name() + "_pen" + to_string(penalizationCoefficient);
		string matrixFilePath = outputDirectory + "/" + fileName + "_A" + terms + ".dat";
		FileMatrix* fileMatrix = new FileMatrix(nUnknowns, nUnknowns, matrixFilePath);

		string rhsFilePath = outputDirectory + "/" + fileName + "_b.dat";
		FileVector* fileRHS = new FileVector(rhsFilePath);

		//--------------------------------------------//
		// Iteration on the elements: diagonal blocks //
		//--------------------------------------------//

		for (BigNumber k = 0; k < grid->Elements.size(); k++)
		{
			Element* element = grid->Elements[k];
			vector<ElementInterface*> elementInterfaces = element->Interfaces;

			for (int localFunctionNumber1 = 0; localFunctionNumber1 < basis->NumberOfLocalFunctionsInElement(element); localFunctionNumber1++)
			{
				IBasisFunction2D* localFunction1 = basis->GetLocalBasisFunction(element, localFunctionNumber1);
				BigNumber basisFunction1 = basis->GlobalFunctionNumber(element, localFunctionNumber1);

				// Current element (block diagonal)
				for (int localFunctionNumber2 = 0; localFunctionNumber2 < basis->NumberOfLocalFunctionsInElement(element); localFunctionNumber2++)
				{
					IBasisFunction2D* localFunction2 = basis->GetLocalBasisFunction(element, localFunctionNumber2);
					BigNumber basisFunction2 = basis->GlobalFunctionNumber(element, localFunctionNumber2);

					double volumicTerm = basis->VolumicTerm(element, localFunction1, localFunction2);

					double coupling = 0;
					double penalization = 0;
					for (vector<ElementInterface*>::iterator it = elementInterfaces.begin(); it != elementInterfaces.end(); ++it)
					{
						ElementInterface* elemInterface = *it;
						coupling += basis->CouplingTerm(elemInterface, element, localFunction1, element, localFunction2);
						penalization += basis->PenalizationTerm(elemInterface, element, localFunction1, element, localFunction2);
					}

					if (terms.compare("volumic") == 0)
						fileMatrix->Add(basisFunction1, basisFunction2, volumicTerm);
					else if (terms.compare("coupling") == 0)
						fileMatrix->Add(basisFunction1, basisFunction2, coupling);
					else if (terms.compare("penalization") == 0)
						fileMatrix->Add(basisFunction1, basisFunction2, penalization);
					else
						fileMatrix->Add(basisFunction1, basisFunction2, volumicTerm + coupling + penalization);
				}

				double rhs = basis->RightHandSide(element, localFunction1);
				fileRHS->Add(rhs);
			}
		}

		//--------------------------------------------------//
		// Iteration on the interfaces: off-diagonal blocks //
		//--------------------------------------------------//

		for (BigNumber k = 0; k < grid->Interfaces.size(); k++)
		{
			auto interface = grid->Interfaces[k];
			if (interface->IsDomainBoundary)
				continue;

			for (int localFunctionNumber1 = 0; localFunctionNumber1 < basis->NumberOfLocalFunctionsInElement(interface->Element1); localFunctionNumber1++)
			{
				IBasisFunction2D* localFunction1 = basis->GetLocalBasisFunction(interface->Element1, localFunctionNumber1);
				BigNumber basisFunction1 = basis->GlobalFunctionNumber(interface->Element1, localFunctionNumber1);
				for (int localFunctionNumber2 = 0; localFunctionNumber2 < basis->NumberOfLocalFunctionsInElement(interface->Element2); localFunctionNumber2++)
				{
					IBasisFunction2D* localFunction2 = basis->GetLocalBasisFunction(interface->Element2, localFunctionNumber2);
					BigNumber basisFunction2 = basis->GlobalFunctionNumber(interface->Element2, localFunctionNumber2);
					double coupling = basis->CouplingTerm(interface, interface->Element1, localFunction1, interface->Element2, localFunction2);
					double penalization = basis->PenalizationTerm(interface, interface->Element1, localFunction1, interface->Element2, localFunction2);
					
					if (terms.compare("volumic") == 0)
					{ }
					else if (terms.compare("coupling") == 0)
					{
						fileMatrix->Add(basisFunction1, basisFunction2, coupling);
						fileMatrix->Add(basisFunction2, basisFunction1, coupling);
					}
					else if (terms.compare("penalization") == 0)
					{
						fileMatrix->Add(basisFunction1, basisFunction2, penalization);
						fileMatrix->Add(basisFunction2, basisFunction1, penalization);
					}
					else
					{
						fileMatrix->Add(basisFunction1, basisFunction2, coupling + penalization);
						fileMatrix->Add(basisFunction2, basisFunction1, coupling + penalization);
					}
				}
			}
		}

		delete fileMatrix;
		delete fileRHS;

		cout << "Matrix exported to \t" << matrixFilePath << endl;
		cout << "RHS exported to \t" << rhsFilePath << endl;
	}

	~Poisson2D()
	{
	}
};

