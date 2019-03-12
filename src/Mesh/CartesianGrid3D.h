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
				SquareFace* bottomBoundary = new SquareFace(numberInterface++, h, this->Elements[index(ix, iy, 0)]);
				this->Faces.push_back(bottomBoundary);
				this->BoundaryFaces.push_back(bottomBoundary);
				dynamic_cast<Cube*>(this->Elements[index(ix, iy, 0)])->SetBottomInterface(bottomBoundary);

				// Top boundary
				SquareFace* topBoundary = new SquareFace(numberInterface++, h, this->Elements[index(ix, iy, n-1)]);
				this->Faces.push_back(topBoundary);
				this->BoundaryFaces.push_back(topBoundary);
				dynamic_cast<Cube*>(this->Elements[index(ix, iy, n - 1)])->SetTopInterface(topBoundary);
			}
		}

		for (BigNumber iz = 0; iz < n; ++iz)
		{
			for (BigNumber ix = 0; ix < n; ++ix)
			{
				// Front boundary
				SquareFace* frontBoundary = new SquareFace(numberInterface++, h, this->Elements[index(ix, 0, iz)]);
				this->Faces.push_back(frontBoundary);
				this->BoundaryFaces.push_back(frontBoundary);
				dynamic_cast<Cube*>(this->Elements[index(ix, 0, iz)])->SetFrontInterface(frontBoundary);

				// Back boundary
				SquareFace* backBoundary = new SquareFace(numberInterface++, h, this->Elements[index(ix, n-1, iz)]);
				this->Faces.push_back(backBoundary);
				this->BoundaryFaces.push_back(backBoundary);
				dynamic_cast<Cube*>(this->Elements[index(ix, n - 1, iz)])->SetBackInterface(backBoundary);
			}
		}

		for (BigNumber iz = 0; iz < n; ++iz)
		{
			for (BigNumber iy = 0; iy < n; ++iy)
			{
				// Left boundary
				SquareFace* leftBoundary = new SquareFace(numberInterface++, h, this->Elements[index(0, iy, iz)]);
				this->Faces.push_back(leftBoundary);
				this->BoundaryFaces.push_back(leftBoundary);
				dynamic_cast<Cube*>(this->Elements[index(0, iy, iz)])->SetLeftInterface(leftBoundary);

				// Right boundary
				SquareFace* rightBoundary = new SquareFace(numberInterface++, h, this->Elements[index(n-1, iy, iz)]);
				this->Faces.push_back(rightBoundary);
				this->BoundaryFaces.push_back(rightBoundary);
				dynamic_cast<Cube*>(this->Elements[index(n - 1, iy, iz)])->SetRightInterface(rightBoundary);
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
						// Right
						Cube* rightNeighbour = dynamic_cast<Cube*>(this->Elements[index(ix+1, iy, iz)]);
						SquareFace* interface = new SquareFace(numberInterface++, h, element, rightNeighbour);
						this->Faces.push_back(interface);
						this->InteriorFaces.push_back(interface);
						element->SetRightInterface(interface);
						rightNeighbour->SetLeftInterface(interface);
					}
					if (iy != n - 1)
					{
						// Back
						Cube* backNeighbour = dynamic_cast<Cube*>(this->Elements[index(ix, iy+1, iz)]);
						SquareFace* interface = new SquareFace(numberInterface++, h, element, backNeighbour);
						this->Faces.push_back(interface);
						this->InteriorFaces.push_back(interface);
						element->SetBackInterface(interface);
						backNeighbour->SetFrontInterface(interface);
					}
					if (iz != n - 1)
					{
						// Top
						Cube* topNeighbour = dynamic_cast<Cube*>(this->Elements[index(ix, iy, iz+1)]);
						SquareFace* interface = new SquareFace(numberInterface++, h, element, topNeighbour);
						this->Faces.push_back(interface);
						this->InteriorFaces.push_back(interface);
						element->SetTopInterface(interface);
						topNeighbour->SetBottomInterface(interface);
					}
				}
			}
		}

	}

private:
	BigNumber index(BigNumber ix, BigNumber iy, BigNumber iz)
	{
		BigNumber n = this->N;
		return iz * n*n + iy * n + ix;
	}
};