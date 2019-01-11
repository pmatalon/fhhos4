#pragma once
#include "IMesh.h"
#include "Interval.h"

class CartesianGrid1D : public IMesh
{
private:
public:
	CartesianGrid1D(BigNumber n) : IMesh(1, n)
	{
		this->Elements.reserve(n);
		this->Interfaces.reserve(n + 1);
		double h = (double)1 / n;

		for (BigNumber k = 0; k < n + 1; k++)
		{
			Element1DInterface* point = new Element1DInterface(k, k * h);
			this->Interfaces.push_back(point);
		}

		for (BigNumber k = 0; k < n; k++)
		{
			Element1DInterface* leftPoint = dynamic_cast<Element1DInterface*>(this->Interfaces[k]);
			Element1DInterface* rightPoint = dynamic_cast<Element1DInterface*>(this->Interfaces[k+1]);
			Interval* element = new Interval(k, leftPoint, rightPoint);
			this->Elements.push_back(element);
		}

		for (BigNumber k = 0; k < n + 1; k++)
		{
			Element1DInterface* point = dynamic_cast<Element1DInterface*>(this->Interfaces[k]);

			if (k == 0)
			{
				point->IsDomainBoundary = true;
				point->Element1 = this->Elements[k];
			}
			else if (k == n)
			{
				point->IsDomainBoundary = true;
				point->Element2 = this->Elements[k - 1];
			}
			else
			{
				point->Element1 = this->Elements[k - 1];
				point->Element2 = this->Elements[k];
			}
		}
	}
};