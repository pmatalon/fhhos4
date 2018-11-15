#pragma once
#include "IMesh.h"
#include "Interval.h"

class CartesianGrid1DNew : public IMesh
{
private:
public:
	CartesianGrid1DNew(BigNumber n) : IMesh(1, n)
	{
		this->Elements.reserve(n * n);
		double h = (double)1 / n;
		for (BigNumber k = 0; k < n + 1; k++)
			this->Elements.push_back(new Interval(k, k/h, (k+1)/h));
	}
};