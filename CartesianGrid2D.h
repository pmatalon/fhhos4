#pragma once
#include <vector>
#include "Element.h"
#include "Square.h"
#include "Face2D.h"
#include "IMesh.h"

using namespace std;

class CartesianGrid2D : public IMesh
{
public:
	CartesianGrid2D(BigNumber n) : IMesh(2, n)
	{
		//----------//
		// Elements //
		//----------//

		this->Elements.reserve(n * n);
		double h = 1 / (double)n;
		for (BigNumber i = 0; i < n; ++i)
		{
			for (BigNumber j = 0; j < n; ++j)
			{
				Square* square = new Square(i*n + j, (double)j / n, (double)i / n, h);
				this->Elements.push_back(square);
			}
		}

		//------------//
		// Faces //
		//------------//

		this->Faces.reserve(n * (n + 3));
		BigNumber numberInterface = 0;

		for (BigNumber j = 0; j < n; ++j)
		{
			// South boundary
			Face2D* southBoundary = new Face2D(numberInterface++, this->Elements[j]);
			this->Faces.push_back(southBoundary);
			//this->BoundaryInterfaces.push_back(southBoundary);
			dynamic_cast<Square*>(this->Elements[j])->SetSouthInterface(southBoundary);

			// North boundary
			Face2D* northBoundary = new Face2D(numberInterface++, this->Elements[(n-1)*n + j]);
			this->Faces.push_back(northBoundary);
			//this->BoundaryInterfaces.push_back(northBoundary);
			dynamic_cast<Square*>(this->Elements[(n - 1)*n + j])->SetNorthInterface(northBoundary);
		}

		for (BigNumber i = 0; i < n; ++i)
		{
			// West boundary
			Face2D* westBoundary = new Face2D(numberInterface++, this->Elements[i*n]);
			this->Faces.push_back(westBoundary);
			//this->BoundaryInterfaces.push_back(westBoundary);
			dynamic_cast<Square*>(this->Elements[i*n])->SetWestInterface(westBoundary);

			// East boundary
			Face2D* eastBoundary = new Face2D(numberInterface++, this->Elements[i*n + n-1]);
			this->Faces.push_back(eastBoundary);
			//this->BoundaryInterfaces.push_back(eastBoundary);
			dynamic_cast<Square*>(this->Elements[i*n + n - 1])->SetEastInterface(eastBoundary);
		}

		for (BigNumber i = 0; i < n; i++)
		{
			for (BigNumber j = 0; j < n; j++)
			{
				Square* element = dynamic_cast<Square*>(this->Elements[i*n + j]);
				if (j != n - 1)
				{
					// East
					Square* eastNeighbour = dynamic_cast<Square*>(this->Elements[i*n + j + 1]);
					Face2D* interface = new Face2D(numberInterface++, element, eastNeighbour);
					this->Faces.push_back(interface);
					element->SetEastInterface(interface);
					eastNeighbour->SetWestInterface(interface);
				}
				if (i != n - 1)
				{
					// North
					Square* northNeighbour = dynamic_cast<Square*>(this->Elements[(i+1)*n + j]);
					Face2D* interface = new Face2D(numberInterface++, element, northNeighbour);
					this->Faces.push_back(interface);
					element->SetNorthInterface(interface);
					northNeighbour->SetSouthInterface(interface);
				}
			}
		}
	}

	~CartesianGrid2D() override
	{
	}

};