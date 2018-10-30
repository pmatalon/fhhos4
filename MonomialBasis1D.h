#pragma once
#include "FunctionalBasis1D.h"
#include "Monomial1D.h"

using namespace std;

class MonomialBasis1D : public FunctionalBasis1D
{
private:
	int _maxPolynomialDegree;

public:
	MonomialBasis1D(int maxPolynomialDegree, CartesianGrid1D* grid, int penalizationCoefficient, function<double(double)> sourceFunction)
		:FunctionalBasis1D(grid, penalizationCoefficient, sourceFunction)
	{
		this->_maxPolynomialDegree = maxPolynomialDegree;

		for (int i = 0; i <= maxPolynomialDegree; i++)
			this->_localFunctions[i] = new Monomial1D(i);
	}

	string Name()
	{
		return "monomials_p" + std::to_string(this->_maxPolynomialDegree);
	}

	double VolumicTerm(BigNumber element, int localFunctionNumber1, int localFunctionNumber2)
	{
		Monomial1D* func1 = (Monomial1D*)this->_localFunctions[localFunctionNumber1];
		Monomial1D* func2 = (Monomial1D*)this->_localFunctions[localFunctionNumber2];

		int i = func1->Degree;
		int j = func2->Degree;

		if (i == 0 || j == 0)
			return 0;
		return (double)(i * j) / (double)(i + j - 1) * (pow(this->_grid->XRight(element), i + j - 1) - pow(this->_grid->XLeft(element), i + j - 1));
	}
};