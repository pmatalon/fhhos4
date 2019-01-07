#pragma once
#include <vector>
#include "Element.h"
#include "Square.h"
#include "Element2DInterface.h"
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
		// Interfaces //
		//------------//

		this->Interfaces.reserve(n * (n + 3));
		BigNumber numberInterface = 0;

		for (BigNumber j = 0; j < n; ++j)
		{
			// South boundary
			Element2DInterface* southBoundary = new Element2DInterface(numberInterface++, this->Elements[j]);
			this->Interfaces.push_back(southBoundary);
			//this->BoundaryInterfaces.push_back(southBoundary);
			dynamic_cast<Square*>(this->Elements[j])->SetSouthInterface(southBoundary);

			// North boundary
			Element2DInterface* northBoundary = new Element2DInterface(numberInterface++, this->Elements[(n-1)*n + j]);
			this->Interfaces.push_back(northBoundary);
			//this->BoundaryInterfaces.push_back(northBoundary);
			dynamic_cast<Square*>(this->Elements[(n - 1)*n + j])->SetNorthInterface(northBoundary);
		}

		for (BigNumber i = 0; i < n; ++i)
		{
			// West boundary
			Element2DInterface* westBoundary = new Element2DInterface(numberInterface++, this->Elements[i*n]);
			this->Interfaces.push_back(westBoundary);
			//this->BoundaryInterfaces.push_back(westBoundary);
			dynamic_cast<Square*>(this->Elements[i*n])->SetWestInterface(westBoundary);

			// East boundary
			Element2DInterface* eastBoundary = new Element2DInterface(numberInterface++, this->Elements[i*n + n-1]);
			this->Interfaces.push_back(eastBoundary);
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
					Element2DInterface* interface = new Element2DInterface(numberInterface++, element, eastNeighbour);
					this->Interfaces.push_back(interface);
					element->SetEastInterface(interface);
					eastNeighbour->SetWestInterface(interface);
				}
				if (i != n - 1)
				{
					// North
					Square* northNeighbour = dynamic_cast<Square*>(this->Elements[(i+1)*n + j]);
					Element2DInterface* interface = new Element2DInterface(numberInterface++, element, northNeighbour);
					this->Interfaces.push_back(interface);
					element->SetNorthInterface(interface);
					northNeighbour->SetSouthInterface(interface);
				}
			}
		}
	}

	~CartesianGrid2D() override
	{
		for (size_t i = 0; i < this->Elements.size(); ++i)
			delete this->Elements[i];
		this->Elements.clear();

		for (size_t i = 0; i < this->Interfaces.size(); ++i)
			delete this->Interfaces[i];
		this->Interfaces.clear();
	}

};