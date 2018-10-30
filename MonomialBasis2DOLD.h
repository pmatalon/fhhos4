#pragma once
#include <functional>
#include <math.h>
#include "CartesianGrid2D.h"
#include "Element.h"
#include "Square.h"
#include "ElementInterface.h"
#include "FunctionalBasisWithObjects.h"
#include "Utils.h"

class MonomialBasis2DOLD : public FunctionalBasisWithObjects
{
private:
	int _maxPolynomialDegree;
	int _penalizationCoefficient;
	std::function<double(double, double)> _sourceFunction;
	CartesianGrid2D* _grid;
	int _numberOfLocalFunctions = 0;
	int** _localFunctions;

public:
	MonomialBasis2DOLD(int maxPolynomialDegree, CartesianGrid2D* grid, int penalizationCoefficient, std::function<double(double, double)> sourceFunction)
	{
		this->_maxPolynomialDegree = maxPolynomialDegree;
		this->_grid = grid;
		this->_penalizationCoefficient = penalizationCoefficient;
		this->_sourceFunction = sourceFunction;

		for (int i = 0; i <= this->_maxPolynomialDegree; i++)
			for (int j = 0; i + j <= this->_maxPolynomialDegree; j++)
				this->_numberOfLocalFunctions++;

		this->_localFunctions = new int*[this->_numberOfLocalFunctions];
		int functionNumber = 0;
		for (int i = 0; i <= this->_maxPolynomialDegree; i++)
		{
			for (int j = 0; i + j <= this->_maxPolynomialDegree; j++)
			{
				int* monomialDegrees = new int[2];
				monomialDegrees[0] = i;
				monomialDegrees[1] = j;
				this->_localFunctions[functionNumber++] = monomialDegrees;
			}
		}
	}

	std::string Name()
	{
		return "oldmonomials_p" + std::to_string(this->_maxPolynomialDegree);
	}

	~MonomialBasis2DOLD() 
	{
		delete[] _localFunctions;
	}

	int NumberOfLocalFunctionsInElement(Element* element)
	{
		return this->_numberOfLocalFunctions;
	}

	BigNumber GlobalFunctionNumber(Element* element, int localFunctionNumber)
	{
		BigNumber numberOfFunctions = NumberOfLocalFunctionsInElement(0);
		return element->Number * numberOfFunctions + localFunctionNumber + 1; // +1 so that the numbers start at 1
	}

	double VolumicTerm(Element* element, int localFunctionNumber1, int localFunctionNumber2)
	{
		Square* square = static_cast<Square*>(element);
		return this->VolumicTerm(square, localFunctionNumber1, localFunctionNumber2);
	}

	double VolumicTerm(Square* element, int localFunctionNumber1, int localFunctionNumber2)
	{
		int* monomialDegrees1 = this->_localFunctions[localFunctionNumber1];
		int* monomialDegrees2 = this->_localFunctions[localFunctionNumber2];
		int i1 = monomialDegrees1[0]; // monomial degree 1 for X
		int j1 = monomialDegrees1[1]; // monomial degree 1 for Y
		int i2 = monomialDegrees2[0]; // monomial degree 2 for X
		int j2 = monomialDegrees2[1]; // monomial degree 2 for Y

		double integralDerivX = 0;
		if (i1 == 0 || i2 == 0)
			integralDerivX = 0;
		else
		{
			double integralOverX = (double)(i1 * i2) / (double)(i1 + i2 - 1) * (pow(element->X + element->Width, i1 + i2 - 1) - pow(element->X, i1 + i2 - 1));
			double integralOverY = 1/(double)(j1 + j2 + 1) * (pow(element->Y + element->Width, j1 + j2 + 1) - pow(element->Y, j1 + j2 + 1));
			integralDerivX = integralOverX * integralOverY;
		}

		double integralDerivY = 0;
		if (j1 == 0 || j2 == 0)
			integralDerivY = 0;
		else
		{
			double integralOverY = (double)(j1 * j2) / (double)(j1 + j2 - 1) * (pow(element->Y + element->Width, j1 + j2 - 1) - pow(element->Y, j1 + j2 - 1));
			double integralOverX = 1 / (double)(i1 + i2 + 1) * (pow(element->X + element->Width, i1 + i2 + 1) - pow(element->X, i1 + i2 + 1));
			integralDerivY = integralOverY * integralOverX;
		}

		return integralDerivX + integralDerivY;
	}

	double CouplingTerm(ElementInterface* interface, Element* element1, int localFunctionNumber1, Element* element2, int localFunctionNumber2)
	{
		return 0;
	}
	
	double PenalizationTerm(ElementInterface* interface, Element* element1, int localFunctionNumber1, Element* element2, int localFunctionNumber2)
	{
		return 0;
	}

	double RightHandSide(Element* element, int localFunctionNumber)
	{
		Square* square = static_cast<Square*>(element);
		return this->RightHandSide(square, localFunctionNumber);
	}

	double RightHandSide(Square* element, int localFunctionNumber)
	{
		int* monomialDegrees = this->_localFunctions[localFunctionNumber];
		int degreeX = localFunctionNumber = monomialDegrees[0];
		int degreeY = localFunctionNumber = monomialDegrees[1];
		std::function<double(double, double)> sourceTimesBasisFunction = [this, degreeX, degreeY](double x, double y) {
			return this->_sourceFunction(x, y) * pow(x, degreeX) * pow(y, degreeY);
		};
		return Utils::Integral(sourceTimesBasisFunction, element->X, element->X + element->Width, element->Y, element->Y + element->Width);
	}

};