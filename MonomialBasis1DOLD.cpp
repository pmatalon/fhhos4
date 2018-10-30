#include "MonomialBasis1DOLD.h"
#include <math.h>
#include "gauss_legendre.h"
#include "CartesianGrid1D.h"

MonomialBasis1DOLD::MonomialBasis1DOLD(int maxPolynomialDegree, CartesianGrid1D* grid, int penalizationCoefficient, std::function<double(double)> sourceFunction)
{
	this->_maxPolynomialDegree = maxPolynomialDegree;
	this->_grid = grid;
	this->_penalizationCoefficient = penalizationCoefficient;
	this->_sourceFunction = sourceFunction;

	// Interval [0,       1/n]: B(   0)   = 1, B(   1) = X,   B(   2) = X^2,  ... , B(   p) = X^p
	// Interval [1/n,     2/n]: B( p+1)   = 1, B( p+2) = X,   B( p+3) = X^2,  ... , B(2p+1) = X^p
	// Interval [2/n,     3/n]: B(2(p+1)) = 1, B(2(p+1)+1) = X, B(2(p+1)+1) = X^2,  ... , B(2(p+1)+p) = X^p
	// ...
	// Interval [k/n, (k+1)/n]: B(k(p+1)) = 1, B(k(p+1)+1) = X, B(k(p+1)+2) = X^2,  ... , B(k(p+1)+p) = X^p
	// ...
	// Interval [(n-1)/n,   1]: 
}

std::string MonomialBasis1DOLD::Name()
{
	return "oldmonomials_p" + std::to_string(this->_maxPolynomialDegree);
}

int MonomialBasis1DOLD::NumberOfLocalFunctionsInElement(BigNumber element)
{
	return this->_maxPolynomialDegree + 1;
}

BigNumber MonomialBasis1DOLD::GlobalFunctionNumber(BigNumber element, int localFunctionNumber)
{
	return element * (this->_maxPolynomialDegree + 1) + localFunctionNumber + 1; // +1 so that the numbers start at 1
}

double MonomialBasis1DOLD::VolumicTerm(BigNumber element, int localFunctionNumber1, int localFunctionNumber2)
{
	int i = localFunctionNumber1; // monomial degree 1
	int j = localFunctionNumber2; // monomial degree 2

	if (i == 0 || j == 0)
		return 0;
	return (double)(i * j) / (double)(i + j - 1) * (pow(this->_grid->XRight(element), i + j - 1) - pow(this->_grid->XLeft(element), i + j - 1));
}

double MonomialBasis1DOLD::CouplingTerm(BigNumber interface, BigNumber element1, int localFunctionNumber1, BigNumber element2, int localFunctionNumber2)
{
	int i = localFunctionNumber1; // monomial degree 1
	int j = localFunctionNumber2; // monomial degree 2

	// Boundary points
	if (this->_grid->IsBoundaryLeft(interface))
	{
		if (this->_grid->IsFirstElement(element1) && this->_grid->IsFirstElement(element2))
		{
			if ((i == 0 && j == 1) || (i == 1 && j == 0))
				return 1;
			return 0;
		}
		return 0;
	}
	else if (this->_grid->IsBoundaryRight(interface))
	{
		if (this->_grid->IsLastElement(element1) && this->_grid->IsLastElement(element2))
			return -i-j;
		return 0;
	}

	double x = this->_grid->X(interface);

	// Interior points
	if (element1 == element2)
	{
		int factor = this->_grid->IsLeftInterface(element1, interface) ? 1 : -1;

		if (i == 0 && j == 0)
			return 0;
		else if (i == 0 && j != 0)
			return factor * 0.5 * j * pow(x, j - 1);
		else if (i != 0 && j == 0)
			return factor * 0.5 * i * pow(x, i - 1);
		else
			return factor * (double)(i + j) / 2 * pow(x, i + j - 1);
	}
	else if (element1 == element2 - 1)
	{
		if (i == 0 && j == 0)
			return 0;
		else if (i == 0 && j != 0)
			return -0.5 * j * pow(x, j - 1);
		else if (i != 0 && j == 0)
			return 0.5 * i * pow(x, i - 1);
		else
			return (double)(i - j) / 2 * pow(x, i + j - 1);
	}
	else if (element1 == element2 + 1)
	{
		return this->CouplingTerm(interface, element2, localFunctionNumber2, element1, localFunctionNumber1);
	}
	return 0;
}

double MonomialBasis1DOLD::PenalizationTerm(BigNumber point, BigNumber element1, int localFunctionNumber1, BigNumber element2, int localFunctionNumber2)
{
	int i = localFunctionNumber1; // monomial degree 1
	int j = localFunctionNumber2; // monomial degree 2

	if (element1 == element2)
	{
		return this->_penalizationCoefficient * pow(this->_grid->X(point), i + j);
	}
	else if (element2 == element1 - 1)
	{
		return -this->_penalizationCoefficient * pow(this->_grid->X(point), i + j);
	}
	else if (element2 == element1 + 1)
	{
		return this->PenalizationTerm(point, element2, localFunctionNumber2, element1, localFunctionNumber1);
	}
	else
		return 0;
}

double MonomialBasis1DOLD::RightHandSide(BigNumber element, int localFunctionNumber)
{
	int degree = localFunctionNumber;
	std::function<double (double)> sourceTimesBasisFunction = [this, degree](double x) {
		return this->_sourceFunction(x) * pow(x, degree);
	};
	return integral(this->_grid->XLeft(element), this->_grid->XRight(element), sourceTimesBasisFunction);
}

MonomialBasis1DOLD::~MonomialBasis1DOLD()
{
}
