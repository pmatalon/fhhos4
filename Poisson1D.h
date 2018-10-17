#pragma once
#include <functional>
#include "CartesianGrid1D.h"

class Poisson1D
{
private:
	std::function<double(double)> _sourceFunction;
	CartesianGrid1D* _grid;
public:
	Poisson1D(int n, std::function<double(double)> sourceFunction);
	~Poisson1D();
	void DiscretizeDG(int maxPolynomialDegree, int penalizationCoefficient);
};

