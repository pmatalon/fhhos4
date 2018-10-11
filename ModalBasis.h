#pragma once
#include <functional>

class ModalBasis
{
private:
	int _maxPolynomialDegree;
	int _penalizationCoefficient;
	std::function<double(double)> _sourceFunction;
public:
	ModalBasis(int degree, int nIntervals, int penalizationCoefficient, std::function<double(double)> sourceFunction);
	~ModalBasis();
	double VolumicTerm(int basisFunctionNumber1, int basisFunctionNumber2, double xLeft, double xRight);
	double CouplingTerm(double x, int basisFunctionNumber1, int basisFunctionNumber2);
	double PenalizationTerm(double x, int basisFunctionNumber1, int basisFunctionNumber2);
	double RightHandSide(double xLeft, double xRight, int basisFunctionNumber);
};

