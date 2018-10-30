#pragma once
#include <functional>
#include "FunctionalBasisWithNumbers.h"
#include "CartesianGrid1D.h"
#include "Element.h"
using namespace std;

class Poisson1D
{
private:
	std::function<double(double)> _sourceFunction;
public:
	Poisson1D(std::function<double(double)> sourceFunction);
	~Poisson1D();
	void DiscretizeDG(CartesianGrid1D* grid, FunctionalBasisWithNumbers* basis, int penalizationCoefficient, string outputDirectory);
};

