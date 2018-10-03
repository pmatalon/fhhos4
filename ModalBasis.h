#pragma once
class ModalBasis
{
private:
	int _degree;
public:
	ModalBasis(int degree, int nIntervals);
	~ModalBasis();
	double IntegralOfGradients(int basisFunctionNumber1, int basisFunctionNumber2, double xLeft, double xRight);
};

