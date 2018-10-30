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
	function<double(double, double)> _sourceFunction;
public:
	Poisson2D(function<double(double, double)> sourceFunction)
	{
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

		string fileName = "Poisson2D_n" + to_string(grid->N) + "_DG_SIPG_" + basis->Name() + "_pen" + to_string(penalizationCoefficient);
		string matrixFilePath = outputDirectory + "/" + fileName + "_A" + terms + ".dat";
		FileMatrix* fileMatrix = new FileMatrix(nUnknowns, nUnknowns, matrixFilePath);

		string rhsFilePath = outputDirectory + "/" + fileName + "_b.dat";
		FileVector* fileRHS = new FileVector(rhsFilePath);

		/*function<double(double)> testFunc1D = [](double x) {
			return x;
		};

		function<double(double, double)> testFunc2D = [](double x, double y) {
			return x * y;
		};

		double intTest1D = Utils::Integral(testFunc1D, 1, 4);
		cout << "Integral test 1D = " << intTest1D << endl;
		double intTest2D = Utils::Integral(testFunc2D, 1, 2, 1, 1);
		cout << "Integral test 2D = " << intTest2D << endl;*/

		//--------------------------------------------//
		// Iteration on the elements: diagonal blocks //
		//--------------------------------------------//

		for (BigNumber k = 0; k < grid->Elements.size(); k++)
		{
			Element* element = grid->Elements[k];
			vector<ElementInterface*> elementInterfaces = element->Interfaces;

			for (int localFunction1 = 0; localFunction1 < basis->NumberOfLocalFunctionsInElement(element); localFunction1++)
			{
				BigNumber basisFunction1 = basis->GlobalFunctionNumber(element, localFunction1);

				// Current element (block diagonal)
				for (int localFunction2 = 0; localFunction2 < basis->NumberOfLocalFunctionsInElement(element); localFunction2++)
				{
					BigNumber basisFunction2 = basis->GlobalFunctionNumber(element, localFunction2);

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

			for (int localFunction1 = 0; localFunction1 < basis->NumberOfLocalFunctionsInElement(interface->Element1); localFunction1++)
			{
				BigNumber basisFunction1 = basis->GlobalFunctionNumber(interface->Element1, localFunction1);
				for (int localFunction2 = 0; localFunction2 < basis->NumberOfLocalFunctionsInElement(interface->Element2); localFunction2++)
				{
					BigNumber basisFunction2 = basis->GlobalFunctionNumber(interface->Element2, localFunction2);
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

