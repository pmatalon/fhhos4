#include "ModalBasis.h"
#include <math.h>

ModalBasis::ModalBasis(int degree, int nIntervals)
{
	this->_degree = degree;

	// Interval [0,       1/n]: B(   0) = 1, B(   1) = X,   B(   2) = X^2,  ... , B(   p) = X^p
	// Interval [1/n,     2/n]: B( p+1) = 1, B( p+2) = X,   B( p+3) = X^2,  ... , B(2p+1) = X^p
	// Interval [2/n,     3/n]: B(2p+2) = 1, B(2p+3) = X,   B(2p+4) = X^2,  ... , B(3p+2) = X^p
	// ...
	// Interval [k/n, (k+1)/n]: B(k(p+1))=1, B(k(p+1)+1)=X, B(k(p+1)+2)=X^2, ..., B((p+1)+p)=X^p
	// ...
	// Interval [(n-1)/n,   1]: 
}

inline double ModalBasis::IntegralOfGradients(int basisFunctionNumber1, int basisFunctionNumber2, double xLeft, double xRight)
{
	if (abs(basisFunctionNumber2 - basisFunctionNumber1) > this->_degree)
		return 0;
	int polyDegree1 = 0;
	int p = basisFunctionNumber1 + basisFunctionNumber2;
	// Primitive of X^p = X^(p+1)/(p+1)
	double integral = (pow(xRight, p + 1) - pow(xLeft, p + 1)) / (p + 1);
	return integral;
}


ModalBasis::~ModalBasis()
{
}
