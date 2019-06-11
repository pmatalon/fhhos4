#pragma once
#include <vector>
#include "Element.h"
#include "Parallelepiped.h"
#include "RectangularFace.h"
#include "Mesh.h"

using namespace std;

class CartesianGrid3D : public Mesh<3>
{
public:
	BigNumber Ny;
	BigNumber Nz;

	CartesianGrid3D(BigNumber nx, BigNumber ny, BigNumber nz) : Mesh(nx)
	{
		// nx = ny = nz falls down to cubic elements

		this->Ny = ny;
		this->Nz = nz;

		//----------//
		// Elements //
		//----------//

		this->Elements.reserve(nx * ny * nz);
		double hx = 1 / (double)nx;
		double hy = 1 / (double)ny;
		double hz = 1 / (double)nz;

		for (BigNumber iz = 0; iz < nz; ++iz)
		{
			for (BigNumber iy = 0; iy < ny; ++iy)
			{
				for (BigNumber ix = 0; ix < nx; ++ix)
				{
					Parallelepiped* element = new Parallelepiped(index(ix, iy, iz), ix * hx, iy * hy, iz * hz, hx, hy, hz);
					this->Elements.push_back(element);
				}
			}
		}

		//-------//
		// Faces //
		//-------//

		this->Faces.reserve(nx * (ny + 1) * (nz + 1) + ny * (nx + 1) * (nz + 1) + nz * (nx + 1) * (ny + 1));
		BigNumber numberInterface = 0;

		for (BigNumber iy = 0; iy < ny; ++iy)
		{
			for (BigNumber ix = 0; ix < nx; ++ix)
			{
				// Bottom boundary
				Parallelepiped* element = dynamic_cast<Parallelepiped*>(this->Elements[index(ix, iy, 0)]);
				RectangularFace* bottomBoundary = new RectangularFace(numberInterface++, element->BackLeftBottomCorner, hx, hy, element, CartesianShapeOrientation::InXOY);
				this->Faces.push_back(bottomBoundary);
				this->BoundaryFaces.push_back(bottomBoundary);
				element->SetBottomFace(bottomBoundary);

				// Top boundary
				element = dynamic_cast<Parallelepiped*>(this->Elements[index(ix, iy, nz - 1)]);
				RectangularFace* topBoundary = new RectangularFace(numberInterface++, element->BackLeftTopCorner, hx, hy, element, CartesianShapeOrientation::InXOY);
				this->Faces.push_back(topBoundary);
				this->BoundaryFaces.push_back(topBoundary);
				element->SetTopFace(topBoundary);
			}
		}

		for (BigNumber iz = 0; iz < nz; ++iz)
		{
			for (BigNumber ix = 0; ix < nx; ++ix)
			{
				// Left boundary
				Parallelepiped* element = dynamic_cast<Parallelepiped*>(this->Elements[index(ix, 0, iz)]);
				RectangularFace* leftBoundary = new RectangularFace(numberInterface++, element->BackLeftBottomCorner, hx, hz, element, CartesianShapeOrientation::InXOZ);
				this->Faces.push_back(leftBoundary);
				this->BoundaryFaces.push_back(leftBoundary);
				element->SetLeftFace(leftBoundary);

				// Right boundary
				element = dynamic_cast<Parallelepiped*>(this->Elements[index(ix, ny - 1, iz)]);
				RectangularFace* rightBoundary = new RectangularFace(numberInterface++, element->BackRightBottomCorner, hx, hz, element, CartesianShapeOrientation::InXOZ);
				this->Faces.push_back(rightBoundary);
				this->BoundaryFaces.push_back(rightBoundary);
				element->SetRightFace(rightBoundary);
			}
		}

		for (BigNumber iz = 0; iz < nz; ++iz)
		{
			for (BigNumber iy = 0; iy < ny; ++iy)
			{
				// Back boundary
				Parallelepiped* element = dynamic_cast<Parallelepiped*>(this->Elements[index(0, iy, iz)]);
				RectangularFace* backBoundary = new RectangularFace(numberInterface++, element->BackLeftBottomCorner, hy, hz, element, CartesianShapeOrientation::InYOZ);
				this->Faces.push_back(backBoundary);
				this->BoundaryFaces.push_back(backBoundary);
				element->SetBackFace(backBoundary);

				// Front boundary
				element = dynamic_cast<Parallelepiped*>(this->Elements[index(nx - 1, iy, iz)]);
				RectangularFace* frontBoundary = new RectangularFace(numberInterface++, element->FrontLeftBottomCorner, hy, hz, element, CartesianShapeOrientation::InYOZ);
				this->Faces.push_back(frontBoundary);
				this->BoundaryFaces.push_back(frontBoundary);
				element->SetFrontFace(frontBoundary);
			}
		}

		for (BigNumber iz = 0; iz < nz; ++iz)
		{
			for (BigNumber iy = 0; iy < ny; ++iy)
			{
				for (BigNumber ix = 0; ix < nx; ++ix)
				{
					Parallelepiped* element = dynamic_cast<Parallelepiped*>(this->Elements[index(ix, iy, iz)]);
					if (ix != nx - 1)
					{
						// Front
						Parallelepiped* frontNeighbour = dynamic_cast<Parallelepiped*>(this->Elements[index(ix+1, iy, iz)]);
						RectangularFace* interface = new RectangularFace(numberInterface++, element->FrontLeftBottomCorner, hy, hz, element, frontNeighbour, CartesianShapeOrientation::InYOZ);
						this->Faces.push_back(interface);
						this->InteriorFaces.push_back(interface);
						element->SetFrontFace(interface);
						frontNeighbour->SetBackFace(interface);
					}
					if (iy != ny - 1)
					{
						// Right
						Parallelepiped* rightNeighbour = dynamic_cast<Parallelepiped*>(this->Elements[index(ix, iy+1, iz)]);
						RectangularFace* interface = new RectangularFace(numberInterface++, element->BackRightBottomCorner, hx, hz, element, rightNeighbour, CartesianShapeOrientation::InXOZ);
						this->Faces.push_back(interface);
						this->InteriorFaces.push_back(interface);
						element->SetRightFace(interface);
						rightNeighbour->SetLeftFace(interface);
					}
					if (iz != nz - 1)
					{
						// Top
						Parallelepiped* topNeighbour = dynamic_cast<Parallelepiped*>(this->Elements[index(ix, iy, iz+1)]);
						RectangularFace* interface = new RectangularFace(numberInterface++, element->BackLeftTopCorner, hx, hy, element, topNeighbour, CartesianShapeOrientation::InXOY);
						this->Faces.push_back(interface);
						this->InteriorFaces.push_back(interface);
						element->SetTopFace(interface);
						topNeighbour->SetBottomFace(interface);
					}
				}
			}
		}

	}

	void BuildCoarserMesh()
	{
		cout << "Error: BuildCoarserMesh not implemented!" << endl;
		exit(EXIT_FAILURE);
	}

private:
	inline BigNumber index(BigNumber ix, BigNumber iy, BigNumber iz)
	{
		BigNumber nx = this->N;
		BigNumber ny = this->Ny;
		return iz * ny*nx + iy * nx + ix;
	}
};