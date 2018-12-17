#include <iostream>
#include <functional>
//#include "CartesianGrid2D.h"
#include "IMesh.h"
#include "FunctionalBasisWithObjects.h"
#include "ElementInterface.h"
#include "FileMatrix.h"
#include "FileVector.h"
#include "MonomialBasis2D.h"
using namespace std;

#pragma once
template <class IBasisFunction>
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

	void DiscretizeDG(IMesh* grid, FunctionalBasisWithObjects<IBasisFunction>* basis, int penalizationCoefficient, string outputDirectory, bool extractMatrixComponents)
	{
		bool autoPenalization = penalizationCoefficient == -1;
		if (autoPenalization)
			penalizationCoefficient = 4 * pow(basis->GetDegree() + 1, 2) * grid->N; // Ralph-Hartmann

		cout << "Discretization: Discontinuous Galerkin SIPG" << endl;
		cout << "\tPenalization coefficient: " << penalizationCoefficient << endl;
		cout << "\tBasis of polynomials: " << basis->Name() << endl;
		
		cout << "Local functions: " << basis->NumberOfLocalFunctionsInElement(0) << endl;
		for (int localFunctionNumber = 0; localFunctionNumber < basis->NumberOfLocalFunctionsInElement(NULL); localFunctionNumber++)
		{
			IBasisFunction* localFunction = basis->GetLocalBasisFunction(NULL, localFunctionNumber);
			cout << "\t " << localFunction->ToString() << endl;
		}
		BigNumber nUnknowns = static_cast<int>(grid->Elements.size()) * basis->NumberOfLocalFunctionsInElement(0);
		cout << "Unknowns: " << nUnknowns << endl;

		string fileName = "Poisson" + to_string(grid->Dim) + "D" + this->_solution + "_n" + to_string(grid->N) + "_DG_SIPG_" + basis->Name() + "_pen" + (autoPenalization ? "-1" : to_string(penalizationCoefficient));
		string matrixFilePath			= outputDirectory + "/" + fileName + "_A.dat";
		string matrixVolumicFilePath	= outputDirectory + "/" + fileName + "_A_volumic.dat";
		string matrixCouplingFilePath	= outputDirectory + "/" + fileName + "_A_coupling.dat";
		string matrixPenFilePath		= outputDirectory + "/" + fileName + "_A_pen.dat";
		FileMatrix* fileMatrix			= new FileMatrix(nUnknowns, nUnknowns, matrixFilePath);
		FileMatrix* fileMatrixVolumic = NULL;
		FileMatrix* fileMatrixCoupling = NULL;
		FileMatrix* fileMatrixPen = NULL;
		if (extractMatrixComponents)
		{
			fileMatrixVolumic	= new FileMatrix(nUnknowns, nUnknowns, matrixVolumicFilePath);
			fileMatrixCoupling	= new FileMatrix(nUnknowns, nUnknowns, matrixCouplingFilePath);
			fileMatrixPen		= new FileMatrix(nUnknowns, nUnknowns, matrixPenFilePath);
		}

		string rhsFilePath = outputDirectory + "/" + fileName + "_b.dat";
		FileVector* fileRHS = new FileVector(rhsFilePath);

		//--------------------------------------------//
		// Iteration on the elements: diagonal blocks //
		//--------------------------------------------//

		for (BigNumber k = 0; k < grid->Elements.size(); k++)
		{
			Element* element = grid->Elements[k];
			cout << "Element " << k << endl;
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

					double volumicTerm = basis->VolumicTerm(element, localFunction1, localFunction2);

					//cout << "\t func1 = " << dynamic_cast<Monomial2D*>(localFunction1)->ToString() << " func2 = " << dynamic_cast<Monomial2D*>(localFunction2)->ToString() << endl;
					cout << "\t\t Volumic = " << volumicTerm << endl;

					double coupling = 0;
					double penalization = 0;
					for (ElementInterface* elemInterface : elementInterfaces)
					{
						coupling += basis->CouplingTerm(elemInterface, element, localFunction1, element, localFunction2);
						penalization += basis->PenalizationTerm(elemInterface, element, localFunction1, element, localFunction2);
						cout << "\t\t " << elemInterface->ToString() << ":\t c=" << coupling << "\tp=" << penalization << endl;
					}

					cout << "\t\t TOTAL = " << volumicTerm + coupling + penalization << endl;
					
					if (extractMatrixComponents)
					{
						fileMatrixVolumic->Add(basisFunction1, basisFunction2, volumicTerm);
						fileMatrixCoupling->Add(basisFunction1, basisFunction2, coupling);
						fileMatrixPen->Add(basisFunction1, basisFunction2, penalization);
					}
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

			cout << interface->ToString() << endl;

			for (int localFunctionNumber1 = 0; localFunctionNumber1 < basis->NumberOfLocalFunctionsInElement(interface->Element1); localFunctionNumber1++)
			{
				IBasisFunction* localFunction1 = basis->GetLocalBasisFunction(interface->Element1, localFunctionNumber1);
				BigNumber basisFunction1 = basis->GlobalFunctionNumber(interface->Element1, localFunctionNumber1);
				for (int localFunctionNumber2 = 0; localFunctionNumber2 < basis->NumberOfLocalFunctionsInElement(interface->Element2); localFunctionNumber2++)
				{
					IBasisFunction* localFunction2 = basis->GetLocalBasisFunction(interface->Element2, localFunctionNumber2);
					BigNumber basisFunction2 = basis->GlobalFunctionNumber(interface->Element2, localFunctionNumber2);
					double coupling = basis->CouplingTerm(interface, interface->Element1, localFunction1, interface->Element2, localFunction2);
					double penalization = basis->PenalizationTerm(interface, interface->Element1, localFunction1, interface->Element2, localFunction2);
					
					if (extractMatrixComponents)
					{
						fileMatrixCoupling->Add(basisFunction1, basisFunction2, coupling);
						fileMatrixCoupling->Add(basisFunction2, basisFunction1, coupling);

						fileMatrixPen->Add(basisFunction1, basisFunction2, penalization);
						fileMatrixPen->Add(basisFunction2, basisFunction1, penalization);
					}
					fileMatrix->Add(basisFunction1, basisFunction2, coupling + penalization);
					fileMatrix->Add(basisFunction2, basisFunction1, coupling + penalization);
				}
			}
		}

		if (extractMatrixComponents)
		{
			delete fileMatrixVolumic;
			delete fileMatrixCoupling;
			delete fileMatrixPen;
		}
		delete fileMatrix;
		delete fileRHS;

		cout << "Matrix exported to \t" << matrixFilePath << endl;
		cout << "RHS exported to \t" << rhsFilePath << endl;
	}
};

