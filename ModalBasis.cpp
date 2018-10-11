#include "ModalBasis.h"
#include <math.h>
#include "gauss_legendre.h"

ModalBasis::ModalBasis(int maxPolynomialDegree, int nIntervals, int penalizationCoefficient, std::function<double(double)> sourceFunction)
{
	this->_maxPolynomialDegree = maxPolynomialDegree;
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

double ModalBasis::VolumicTerm(int basisFunctionNumber1, int basisFunctionNumber2, double xLeft, double xRight)
{
	if (abs(basisFunctionNumber2 - basisFunctionNumber1) > this->_maxPolynomialDegree)
		return 0;
	int polyDegree1 = basisFunctionNumber1 % (this->_maxPolynomialDegree + 1);
	int polyDegree2 = basisFunctionNumber2 % (this->_maxPolynomialDegree + 1);
	if (polyDegree1 == 0 || polyDegree2 == 0)
		return 0;
	double integral = (double)(polyDegree1 * polyDegree2) / (double)(polyDegree1 + polyDegree2 - 1) * (pow(xRight, polyDegree1 + polyDegree2 - 1) - pow(xLeft, polyDegree1 + polyDegree2 - 1));
	return integral;
}

double ModalBasis::CouplingTerm(double x, int basisFunctionNumber1, int basisFunctionNumber2)
{
	//if (abs(basisFunctionNumber2 - basisFunctionNumber1) > 2*(this->_maxPolynomialDegree + 1))
	//	return 0;
	int polyDegree1 = basisFunctionNumber1 % (this->_maxPolynomialDegree + 1);
	int polyDegree2 = basisFunctionNumber2 % (this->_maxPolynomialDegree + 1);

	if (polyDegree1 == 0 && polyDegree2 == 0)
		return -this->_penalizationCoefficient;
	else if (polyDegree1 == 0 && polyDegree2 != 0)
		return -0.5* polyDegree2*pow(x, polyDegree2 - 1) - this->_penalizationCoefficient*pow(x, polyDegree2);
	else if (polyDegree1 != 0 && polyDegree2 == 0)
		return 0.5*polyDegree1*pow(x, polyDegree1 - 1) - this->_penalizationCoefficient*pow(x, polyDegree1);
	else
		return 0.5*(polyDegree1 - polyDegree2)*pow(x, polyDegree1 + polyDegree2 - 1) - this->_penalizationCoefficient*pow(x, polyDegree1 + polyDegree2);
}

double ModalBasis::PenalizationTerm(double x, int basisFunctionNumber1, int basisFunctionNumber2)
{
	int polyDegree1 = basisFunctionNumber1 % (this->_maxPolynomialDegree + 1);
	int polyDegree2 = basisFunctionNumber2 % (this->_maxPolynomialDegree + 1);

	if (polyDegree1 == 0 && polyDegree2 == 0)
		return -this->_penalizationCoefficient;
	else if (polyDegree1 == 0 && polyDegree2 != 0)
		return -0.5* polyDegree2*pow(x, polyDegree2 - 1) - this->_penalizationCoefficient*pow(x, polyDegree2);
	else if (polyDegree1 != 0 && polyDegree2 == 0)
		return 0.5*polyDegree1*pow(x, polyDegree1 - 1) - this->_penalizationCoefficient*pow(x, polyDegree1);
	else
		return 0.5*(polyDegree1 - polyDegree2)*pow(x, polyDegree1 + polyDegree2 - 1) - this->_penalizationCoefficient*pow(x, polyDegree1 + polyDegree2);
}

double ModalBasis::RightHandSide(double xLeft, double xRight, int basisFunctionNumber)
{
	int degree = basisFunctionNumber % (this->_maxPolynomialDegree + 1);
	std::function<double (double)> sourceTimesBasisFunction = [this, degree](double x) {
		return this->_sourceFunction(x) * pow(x, degree);
	};
	return integral(xLeft, xRight, sourceTimesBasisFunction);
}

/*double ModalBasis::SourceTimesBasisFunction(double x)
{
	this->_sourceFunction(x) * pow(x, p);
}*/



ModalBasis::~ModalBasis()
{
}
