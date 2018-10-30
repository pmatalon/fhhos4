#pragma once
#include <functional>
#include "FunctionalBasisWithNumbers.h"
#include "CartesianGrid1D.h"
#include "Element.h"

class MonomialBasis1DOLD : public FunctionalBasisWithNumbers
{
private:
	int _maxPolynomialDegree;
	int _penalizationCoefficient;
	std::function<double(double)> _sourceFunction;
	CartesianGrid1D* _grid;
public:
	MonomialBasis1DOLD(int maxPolynomialDegree, CartesianGrid1D* grid, int penalizationCoefficient, std::function<double(double)> sourceFunction);
	std::string Name();
	int NumberOfLocalFunctionsInElement(BigNumber element);
	BigNumber GlobalFunctionNumber(BigNumber element, int localFunctionNumber);
	double VolumicTerm(BigNumber element, int localFunctionNumber1, int localFunctionNumber2);
	double CouplingTerm(BigNumber interface, BigNumber element1, int localFunctionNumber1, BigNumber element2, int localFunctionNumber2);
	double PenalizationTerm(BigNumber interface, BigNumber element1, int localFunctionNumber1, BigNumber element2, int localFunctionNumber2);
	double RightHandSide(BigNumber element, int localFunctionNumber);
	~MonomialBasis1DOLD();
};

