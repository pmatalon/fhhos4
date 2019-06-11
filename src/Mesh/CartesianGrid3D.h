#pragma once
#include <vector>
#include "Element.h"
#include "Cube.h"
#include "SquareFace.h"
#include "Mesh.h"

using namespace std;

class CartesianGrid3D : public Mesh<3>
{
public:
	CartesianGrid3D(BigNumber n) : Mesh(n)
	{
		//----------//
		// Elements //
		//----------//

		this->Elements.reserve(n * n * n);
		double h = 1 / (double)n;

		for (BigNumber iz = 0; iz < n; ++iz)
		{
			for (BigNumber iy = 0; iy < n; ++iy)
			{
				for (BigNumber ix = 0; ix < n; ++ix)
				{
					Cube* cube = new Cube(index(ix, iy, iz), (double)ix / n, (double)iy / n, (double)iz / n, h);
					this->Elements.push_back(cube);
				}
			}
		}

		//-------//
		// Faces //
		//-------//

		this->Faces.reserve(2*n*n * (n + 1));
		BigNumber numberInterface = 0;

		for (BigNumber iy = 0; iy < n; ++iy)
		{
			for (BigNumber ix = 0; ix < n; ++ix)
			{
				// Bottom boundary
				Cube* cube = dynamic_cast<Cube*>(this->Elements[index(ix, iy, 0)]);
				SquareFace* bottomBoundary = new SquareFace(numberInterface++, cube->BackLeftBottomCorner, h, cube, CartesianShapeOrientation::InXOY);
				this->Faces.push_back(bottomBoundary);
				this->BoundaryFaces.push_back(bottomBoundary);
				cube->SetBottomFace(bottomBoundary);

				// Top boundary
				cube = dynamic_cast<Cube*>(this->Elements[index(ix, iy, n - 1)]);
				SquareFace* topBoundary = new SquareFace(numberInterface++, cube->BackLeftTopCorner, h, this->Elements[index(ix, iy, n-1)], CartesianShapeOrientation::InXOY);
				this->Faces.push_back(topBoundary);
				this->BoundaryFaces.push_back(topBoundary);
				cube->SetTopFace(topBoundary);
			}
		}

		for (BigNumber iz = 0; iz < n; ++iz)
		{
			for (BigNumber ix = 0; ix < n; ++ix)
			{
				// Left boundary
				Cube* cube = dynamic_cast<Cube*>(this->Elements[index(ix, 0, iz)]);
				SquareFace* leftBoundary = new SquareFace(numberInterface++, cube->BackLeftBottomCorner, h, cube, CartesianShapeOrientation::InXOZ);
				this->Faces.push_back(leftBoundary);
				this->BoundaryFaces.push_back(leftBoundary);
				cube->SetLeftFace(leftBoundary);

				// Right boundary
				cube = dynamic_cast<Cube*>(this->Elements[index(ix, n - 1, iz)]);
				SquareFace* rightBoundary = new SquareFace(numberInterface++, cube->BackRightBottomCorner, h, cube, CartesianShapeOrientation::InXOZ);
				this->Faces.push_back(rightBoundary);
				this->BoundaryFaces.push_back(rightBoundary);
				cube->SetRightFace(rightBoundary);
			}
		}

		for (BigNumber iz = 0; iz < n; ++iz)
		{
			for (BigNumber iy = 0; iy < n; ++iy)
			{
				// Back boundary
				Cube* cube = dynamic_cast<Cube*>(this->Elements[index(0, iy, iz)]);
				SquareFace* backBoundary = new SquareFace(numberInterface++, cube->BackLeftBottomCorner, h, cube, CartesianShapeOrientation::InYOZ);
				this->Faces.push_back(backBoundary);
				this->BoundaryFaces.push_back(backBoundary);
				cube->SetBackFace(backBoundary);

				// Front boundary
				cube = dynamic_cast<Cube*>(this->Elements[index(n - 1, iy, iz)]);
				SquareFace* frontBoundary = new SquareFace(numberInterface++, cube->FrontLeftBottomCorner, h, cube, CartesianShapeOrientation::InYOZ);
				this->Faces.push_back(frontBoundary);
				this->BoundaryFaces.push_back(frontBoundary);
				cube->SetFrontFace(frontBoundary);
			}
		}

		for (BigNumber iz = 0; iz < n; ++iz)
		{
			for (BigNumber iy = 0; iy < n; ++iy)
			{
				for (BigNumber ix = 0; ix < n; ++ix)
				{
					Cube* element = dynamic_cast<Cube*>(this->Elements[index(ix, iy, iz)]);
					if (ix != n - 1)
					{
						// Front
						Cube* frontNeighbour = dynamic_cast<Cube*>(this->Elements[index(ix+1, iy, iz)]);
						SquareFace* interface = new SquareFace(numberInterface++, element->FrontLeftBottomCorner, h, element, frontNeighbour, CartesianShapeOrientation::InYOZ);
						this->Faces.push_back(interface);
						this->InteriorFaces.push_back(interface);
						element->SetFrontFace(interface);
						frontNeighbour->SetBackFace(interface);
					}
					if (iy != n - 1)
					{
						// Right
						Cube* rightNeighbour = dynamic_cast<Cube*>(this->Elements[index(ix, iy+1, iz)]);
						SquareFace* interface = new SquareFace(numberInterface++, element->BackRightBottomCorner, h, element, rightNeighbour, CartesianShapeOrientation::InXOZ);
						this->Faces.push_back(interface);
						this->InteriorFaces.push_back(interface);
						element->SetRightFace(interface);
						rightNeighbour->SetLeftFace(interface);
					}
					if (iz != n - 1)
					{
						// Top
						Cube* topNeighbour = dynamic_cast<Cube*>(this->Elements[index(ix, iy, iz+1)]);
						SquareFace* interface = new SquareFace(numberInterface++, element->BackLeftTopCorner, h, element, topNeighbour, CartesianShapeOrientation::InXOY);
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
	BigNumber index(BigNumber ix, BigNumber iy, BigNumber iz)
	{
		BigNumber n = this->N;
		return iz * n*n + iy * n + ix;
	}
};