#pragma once
#include <vector>
#include "Element.h"
#include "Rectangle.h"
#include "IntervalFace.h"
#include "Mesh.h"
using namespace std;

class CartesianGrid2D : public Mesh<2>
{
public:
	BigNumber Ny;

	CartesianGrid2D(BigNumber nx, BigNumber ny) : Mesh(nx)
	{
		// nx = ny falls down to square elements

		this->Ny = ny;

		//----------//
		// Elements //
		//----------//

		this->Elements.reserve(nx * ny);
		double hx = 1 / (double)nx;
		double hy = 1 / (double)ny;
		for (BigNumber i = 0; i < ny; ++i)
		{
			for (BigNumber j = 0; j < nx; ++j)
			{
				Rectangle* rectangle = new Rectangle(i*nx + j, j * hx, i * hy, hx, hy);
				this->Elements.push_back(rectangle);
			}
		}

		//-------//
		// Faces //
		//-------//

		this->Faces.reserve(nx * (ny + 1) + ny * (nx + 1));
		BigNumber numberInterface = 0;

		for (BigNumber j = 0; j < nx; ++j)
		{
			// South boundary
			Rectangle* rectangle = dynamic_cast<Rectangle*>(this->Elements[j]);
			IntervalFace* southBoundary = new IntervalFace(numberInterface++, rectangle->BottomLeftCorner, hx, rectangle, CartesianShapeOrientation::Horizontal);
			this->Faces.push_back(southBoundary);
			this->BoundaryFaces.push_back(southBoundary);
			dynamic_cast<Rectangle*>(this->Elements[j])->SetSouthInterface(southBoundary);

			// North boundary
			rectangle = dynamic_cast<Rectangle*>(this->Elements[(ny - 1)*nx + j]);
			IntervalFace* northBoundary = new IntervalFace(numberInterface++, rectangle->TopLeftCorner, hx, rectangle, CartesianShapeOrientation::Horizontal);
			this->Faces.push_back(northBoundary);
			this->BoundaryFaces.push_back(northBoundary);
			dynamic_cast<Rectangle*>(this->Elements[(ny - 1)*nx + j])->SetNorthInterface(northBoundary);
		}

		for (BigNumber i = 0; i < ny; ++i)
		{
			// West boundary
			Rectangle* rectangle = dynamic_cast<Rectangle*>(this->Elements[i*nx]);
			IntervalFace* westBoundary = new IntervalFace(numberInterface++, rectangle->BottomLeftCorner, hy, rectangle, CartesianShapeOrientation::Vertical);
			this->Faces.push_back(westBoundary);
			this->BoundaryFaces.push_back(westBoundary);
			dynamic_cast<Rectangle*>(this->Elements[i*nx])->SetWestInterface(westBoundary);

			// East boundary
			rectangle = dynamic_cast<Rectangle*>(this->Elements[i*nx + nx - 1]);
			IntervalFace* eastBoundary = new IntervalFace(numberInterface++, rectangle->BottomRightCorner, hy, rectangle, CartesianShapeOrientation::Vertical);
			this->Faces.push_back(eastBoundary);
			this->BoundaryFaces.push_back(eastBoundary);
			dynamic_cast<Rectangle*>(this->Elements[i*nx + nx - 1])->SetEastInterface(eastBoundary);
		}

		for (BigNumber i = 0; i < ny; i++)
		{
			for (BigNumber j = 0; j < nx; j++)
			{
				Rectangle* element = dynamic_cast<Rectangle*>(this->Elements[i*nx + j]);
				if (j != nx - 1)
				{
					// East
					Rectangle* eastNeighbour = dynamic_cast<Rectangle*>(this->Elements[i*nx + j + 1]);
					IntervalFace* interface = new IntervalFace(numberInterface++, eastNeighbour->BottomLeftCorner, hy, element, eastNeighbour, CartesianShapeOrientation::Vertical);
					this->Faces.push_back(interface);
					this->InteriorFaces.push_back(interface);
					element->SetEastInterface(interface);
					eastNeighbour->SetWestInterface(interface);
				}
				if (i != ny - 1)
				{
					// North
					Rectangle* northNeighbour = dynamic_cast<Rectangle*>(this->Elements[(i + 1)*nx + j]);
					IntervalFace* interface = new IntervalFace(numberInterface++, northNeighbour->BottomLeftCorner, hx, element, northNeighbour, CartesianShapeOrientation::Horizontal);
					this->Faces.push_back(interface);
					this->InteriorFaces.push_back(interface);
					element->SetNorthInterface(interface);
					northNeighbour->SetSouthInterface(interface);
				}
			}
		}

	}

	void BuildCoarserMesh()
	{
		BigNumber nx = this->N;
		BigNumber ny = this->Ny;

		if (nx % 2 != 0 || ny % 2 != 0)
		{
			cout << "Error: impossible to build coarser mesh. Nx and Ny must be even: Nx = " << nx << ", Ny = " << ny << "." << endl;
		}
		else
		{
			CartesianGrid2D* coarserMesh = new CartesianGrid2D(nx / 2, ny / 2);

			for (BigNumber i = 0; i < ny; ++i)
			{
				for (BigNumber j = 0; j < nx; ++j)
				{
					Rectangle* fineElement = dynamic_cast<Rectangle*>(this->Elements[i*nx + j]);
					auto coarseElement = coarserMesh->Elements[(i / 2) * coarserMesh->N + j / 2];
					//cout << "fineElement=" << fineElement->Number << " --> coarse=" << coarseElement->Number << endl;
					coarseElement->FinerElements.push_back(fineElement);
					fineElement->CoarserElement = coarseElement;
					if (i % 2 == 0 && !fineElement->NorthFace->IsDomainBoundary)
						coarseElement->FinerFacesRemoved.push_back(fineElement->NorthFace);
					if (j % 2 == 0 && !fineElement->EastFace->IsDomainBoundary)
						coarseElement->FinerFacesRemoved.push_back(fineElement->EastFace);
				}
			}

			/*for (auto e : coarserMesh->Elements)
			{
				cout << "coarseElement " << e->Number << " removed faces: " ;
				for (auto f : e->FinerFacesRemoved)
					cout << f->Number << " ";
				cout << endl;
			}*/

			this->CoarserMesh = coarserMesh;
		}
	}

};