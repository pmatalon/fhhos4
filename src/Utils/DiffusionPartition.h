#pragma once
#include <functional>
#include "../Mesh/Point.h"
using namespace std;

class DiffusionPartition
{
public:
	function<bool(Point)> IsInPart1;
	double Kappa1 = 1;
	double Kappa2 = 1;

	DiffusionPartition(function<bool(Point)> isInPart1, double kappa1, double kappa2)
	{
		this->IsInPart1 = isInPart1;
		this->Kappa1 = kappa1;
		this->Kappa2 = kappa2;
	}

	double Coefficient(Point p)
	{
		double kappa = this->IsInPart1(p) ? this->Kappa1 : this->Kappa2;
		return kappa;
	}
};