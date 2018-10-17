#pragma once
#include <functional>
#include "CartesianGrid1D.h"

class ModalBasis
{
private:
	int _maxPolynomialDegree;
	int _penalizationCoefficient;
	std::function<double(double)> _sourceFunction;
	CartesianGrid1D* _grid;
public:
	ModalBasis(int degree, CartesianGrid1D* grid, int penalizationCoefficient, std::function<double(double)> sourceFunction);
	int NumberOfLocalFunctionsInElement(int element);
	int GlobalFunctionNumber(int element, int localFunctionNumber);
	double VolumicTerm(int element, int localFunctionNumber1, int localFunctionNumber2);
	double CouplingTerm(int interface, int element1, int localFunctionNumber1, int element2, int localFunctionNumber2);
	double PenalizationTerm(int interface, int element1, int localFunctionNumber1, int element2, int localFunctionNumber2);
	double RightHandSide(int element, int localFunctionNumber);
	~ModalBasis();
};

