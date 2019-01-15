#pragma once
#include <vector>
#include "Element.h"
#include "Cube.h"
#include "Face3D.h"
#include "IMesh.h"

using namespace std;

class CartesianGrid3D : public IMesh
{
public:
	CartesianGrid3D(BigNumber n) : IMesh(3, n)
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

		//------------//
		// Faces //
		//------------//

		this->Faces.reserve(2*n*n * (n + 1));
		BigNumber numberInterface = 0;

		for (BigNumber iy = 0; iy < n; ++iy)
		{
			for (BigNumber ix = 0; ix < n; ++ix)
			{
				// Bottom boundary
				Face3D* bottomBoundary = new Face3D(numberInterface++, this->Elements[index(ix, iy, 0)]);
				this->Faces.push_back(bottomBoundary);
				//this->BoundaryInterfaces.push_back(southBoundary);
				dynamic_cast<Cube*>(this->Elements[index(ix, iy, 0)])->SetBottomInterface(bottomBoundary);

				// Top boundary
				Face3D* topBoundary = new Face3D(numberInterface++, this->Elements[index(ix, iy, n-1)]);
				this->Faces.push_back(topBoundary);
				//this->BoundaryInterfaces.push_back(topBoundary);
				dynamic_cast<Cube*>(this->Elements[index(ix, iy, n - 1)])->SetTopInterface(topBoundary);
			}
		}

		for (BigNumber iz = 0; iz < n; ++iz)
		{
			for (BigNumber ix = 0; ix < n; ++ix)
			{
				// Front boundary
				Face3D* frontBoundary = new Face3D(numberInterface++, this->Elements[index(ix, 0, iz)]);
				this->Faces.push_back(frontBoundary);
				//this->BoundaryInterfaces.push_back(frontBoundary);
				dynamic_cast<Cube*>(this->Elements[index(ix, 0, iz)])->SetFrontInterface(frontBoundary);

				// Back boundary
				Face3D* backBoundary = new Face3D(numberInterface++, this->Elements[index(ix, n-1, iz)]);
				this->Faces.push_back(backBoundary);
				//this->BoundaryInterfaces.push_back(backBoundary);
				dynamic_cast<Cube*>(this->Elements[index(ix, n - 1, iz)])->SetBackInterface(backBoundary);
			}
		}

		for (BigNumber iz = 0; iz < n; ++iz)
		{
			for (BigNumber iy = 0; iy < n; ++iy)
			{
				// Left boundary
				Face3D* leftBoundary = new Face3D(numberInterface++, this->Elements[index(0, iy, iz)]);
				this->Faces.push_back(leftBoundary);
				//this->BoundaryInterfaces.push_back(leftBoundary);
				dynamic_cast<Cube*>(this->Elements[index(0, iy, iz)])->SetLeftInterface(leftBoundary);

				// Right boundary
				Face3D* rightBoundary = new Face3D(numberInterface++, this->Elements[index(n-1, iy, iz)]);
				this->Faces.push_back(rightBoundary);
				//this->BoundaryInterfaces.push_back(rightBoundary);
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
						Face3D* interface = new Face3D(numberInterface++, element, rightNeighbour);
						this->Faces.push_back(interface);
						element->SetRightInterface(interface);
						rightNeighbour->SetLeftInterface(interface);
					}
					if (iy != n - 1)
					{
						// Back
						Cube* backNeighbour = dynamic_cast<Cube*>(this->Elements[index(ix, iy+1, iz)]);
						Face3D* interface = new Face3D(numberInterface++, element, backNeighbour);
						this->Faces.push_back(interface);
						element->SetBackInterface(interface);
						backNeighbour->SetFrontInterface(interface);
					}
					if (iz != n - 1)
					{
						// Top
						Cube* topNeighbour = dynamic_cast<Cube*>(this->Elements[index(ix, iy, iz+1)]);
						Face3D* interface = new Face3D(numberInterface++, element, topNeighbour);
						this->Faces.push_back(interface);
						element->SetTopInterface(interface);
						topNeighbour->SetBottomInterface(interface);
					}
				}
			}
		}
	}
	
	~CartesianGrid3D() override
	{
	}
private:
	BigNumber index(BigNumber ix, BigNumber iy, BigNumber iz)
	{
		BigNumber n = this->N;
		return iz * n*n + iy * n + ix;
	}
};