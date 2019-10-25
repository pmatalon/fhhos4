#pragma once
#include "../Mesh.h"
#include "Interval.h"
#include "InterfacePoint.h"

class CartesianGrid1D : public Mesh<1>
{
public:
	BigNumber N;

	CartesianGrid1D(BigNumber n) : Mesh()
	{
		this->N = n;
		this->Elements.reserve(n);
		this->Faces.reserve(n + 1);
		double h = (double)1 / n;

		for (BigNumber k = 0; k < n + 1; k++)
		{
			Vertex* vertex = new Vertex(k, k * h);
			this->Vertices.push_back(vertex);
			InterfacePoint* point = new InterfacePoint(k, vertex);
			this->Faces.push_back(point);
		}

		for (BigNumber k = 0; k < n; k++)
		{
			InterfacePoint* leftFace = dynamic_cast<InterfacePoint*>(this->Faces[k]);
			InterfacePoint* rightFace = dynamic_cast<InterfacePoint*>(this->Faces[k+1]);
			Interval* element = new Interval(k, leftFace->V, rightFace->V);
			element->SetLeftInterface(leftFace);
			element->SetRightInterface(rightFace);
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
				point->Element1 = this->Elements[k - 1];
			}
			else
			{
				this->InteriorFaces.push_back(point);
				point->Element1 = this->Elements[k - 1];
				point->Element2 = this->Elements[k];
			}
		}
	}

	string Description()
	{
		return "Uniform, " + to_string(this->N) + " subdivisions";
	}

	string FileNamePart()
	{
		return "n" + to_string(this->N);
	}

	double H()
	{
		return (double)1 / this->N;
	}

	void CoarsenMesh(CoarseningStrategy strategy)
	{
		cout << "Error: CoarsenMesh not implemented!" << endl;
		exit(EXIT_FAILURE);
	}
};