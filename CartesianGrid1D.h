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
		this->Faces.reserve(n + 1);
		double h = (double)1 / n;

		for (BigNumber k = 0; k < n + 1; k++)
		{
			Face1D* point = new Face1D(k, k * h);
			this->Faces.push_back(point);
		}

		for (BigNumber k = 0; k < n; k++)
		{
			Face1D* leftPoint = dynamic_cast<Face1D*>(this->Faces[k]);
			Face1D* rightPoint = dynamic_cast<Face1D*>(this->Faces[k+1]);
			Interval* element = new Interval(k, leftPoint, rightPoint);
			this->Elements.push_back(element);
		}

		for (BigNumber k = 0; k < n + 1; k++)
		{
			Face1D* point = dynamic_cast<Face1D*>(this->Faces[k]);

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