#pragma once
#include "Mesh.h"
#include "Interval.h"
#include "PointFace.h"

class CartesianGrid1D : public Mesh<1>
{
private:
public:
	CartesianGrid1D(BigNumber n) : Mesh(n)
	{
		this->Elements.reserve(n);
		this->Faces.reserve(n + 1);
		double h = (double)1 / n;

		for (BigNumber k = 0; k < n + 1; k++)
		{
			PointFace* point = new PointFace(k, k * h);
			this->Faces.push_back(point);
		}

		for (BigNumber k = 0; k < n; k++)
		{
			PointFace* leftPoint = static_cast<PointFace*>(this->Faces[k]);
			PointFace* rightPoint = static_cast<PointFace*>(this->Faces[k+1]);
			Interval* element = new Interval(k, leftPoint->X, rightPoint->X, leftPoint, rightPoint);
			this->Elements.push_back(element);
		}

		for (BigNumber k = 0; k < n + 1; k++)
		{
			Face<1>* point = this->Faces[k];

			if (k == 0)
			{
				point->IsDomainBoundary = true;
				this->BoundaryFaces.push_back(point);
				point->Element1 = this->Elements[k];
			}
			else if (k == n)
			{
				point->IsDomainBoundary = true;
				this->BoundaryFaces.push_back(point);
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