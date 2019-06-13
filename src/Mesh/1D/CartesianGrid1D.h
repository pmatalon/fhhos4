#pragma once
#include "../Mesh.h"
#include "Interval.h"
#include "InterfacePoint.h"

class CartesianGrid1D : public Mesh<1>
{
public:
	CartesianGrid1D(BigNumber n) : Mesh(n)
	{
		this->Elements.reserve(n);
		this->Faces.reserve(n + 1);
		double h = (double)1 / n;

		for (BigNumber k = 0; k < n + 1; k++)
		{
			InterfacePoint* point = new InterfacePoint(k, k * h);
			this->Faces.push_back(point);
		}

		for (BigNumber k = 0; k < n; k++)
		{
			InterfacePoint* leftPoint = dynamic_cast<InterfacePoint*>(this->Faces[k]);
			InterfacePoint* rightPoint = dynamic_cast<InterfacePoint*>(this->Faces[k+1]);
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

	void BuildCoarserMesh()
	{
		cout << "Error: BuildCoarserMesh not implemented!" << endl;
		exit(EXIT_FAILURE);
	}
};