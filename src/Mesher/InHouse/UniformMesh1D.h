#pragma once
#include "../../Mesh/Mesh.h"
#include "../../Mesh/1D/Interval.h"
#include "../../Mesh/1D/InterfacePoint.h"
#include "../SegmentGeometry.h"

class UniformMesh1D : public Mesh<1>
{
public:
	BigNumber N;

	UniformMesh1D(BigNumber n) : Mesh()
	{
		this->N = n;

		this->PhysicalParts = SegmentGeometry::PhysicalParts();
		PhysicalGroup<1>* leftPart = this->PhysicalParts[0];
		PhysicalGroup<1>* rightPart = this->PhysicalParts[1];

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
			element->PhysicalPart = k < n / 2 ? leftPart : rightPart;
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
		return "segment-inhouse-cart-n" + to_string(this->N);
	}
	string GeometryDescription()
	{
		return "1D unit segment";
	}

	double H()
	{
		return (double)1 / this->N;
	}

	double Regularity() override
	{
		return -1;
	}

	Mesh<1>* Copy() override
	{
		return new UniformMesh1D(this->N);
	}
};