#pragma once
#include <functional>

class Poisson1D
{
private:
	int _n;
	//double(*_sourceFunction)(double);
	std::function<double(double)> _sourceFunction;
	double* _x;
public:
	Poisson1D(int n, std::function<double(double)> sourceFunction);
	~Poisson1D();
	void DiscretizeDG(int maxPolynomialDegree, int penalizationCoefficient);
};

