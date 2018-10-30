#pragma once
#include <vector>
#include "Element.h"
#include "Square.h"
#include "Element2DInterface.h"

using namespace std;

class CartesianGrid2D
{
public:
	BigNumber N;
	vector<Square*> Elements;
	vector<Element2DInterface*> Interfaces;
	vector<Element2DInterface*> BoundaryInterfaces;
	/*int NElements;
	int NInterfaces;
	int NBoudaryInterfaces;*/
private:
public:
	CartesianGrid2D(BigNumber n)
	{
		this->N = n;

		// [0,1]^2 descretized in 0, 1/n, 2/n, n/n

		//----------//
		// Elements //
		//----------//

		//this->NElements = n * n;
		this->Elements.reserve(n * n);
		for (BigNumber i = 0; i < n; ++i)
		{
			for (BigNumber j = 0; j < n; ++j)
			{
				Square* square = new Square(i*n + j, (double)j / n, (double)i / n, 1 / (double)n);
				this->Elements.push_back(square);
			}
		}

		//------------//
		// Interfaces //
		//------------//

		this->Interfaces.reserve(n * (n + 3));

		for (BigNumber j = 0; j < n; ++j)
		{
			// South boundary
			Element2DInterface* southBoundary = new Element2DInterface(this->Elements[j]);
			this->Interfaces.push_back(southBoundary);
			this->BoundaryInterfaces.push_back(southBoundary);
			this->Elements[j]->SetSouthInterface(southBoundary);

			// North boundary
			Element2DInterface* northBoundary = new Element2DInterface(this->Elements[(n-1)*n + j]);
			this->Interfaces.push_back(northBoundary);
			this->BoundaryInterfaces.push_back(northBoundary);
			this->Elements[(n - 1)*n + j]->SetNorthInterface(northBoundary);
		}

		for (BigNumber i = 0; i < n; ++i)
		{
			// West boundary
			Element2DInterface* westBoundary = new Element2DInterface(this->Elements[i*n]);
			this->Interfaces.push_back(westBoundary);
			this->BoundaryInterfaces.push_back(westBoundary);
			this->Elements[i*n]->SetWestInterface(westBoundary);

			// East boundary
			Element2DInterface* eastBoundary = new Element2DInterface(this->Elements[i*n + n-1]);
			this->Interfaces.push_back(eastBoundary);
			this->BoundaryInterfaces.push_back(eastBoundary);
			this->Elements[i*n + n - 1]->SetEastInterface(eastBoundary);
		}

		for (BigNumber i = 0; i < n; i++)
		{
			for (BigNumber j = 0; j < n; j++)
			{
				Square* element = this->Elements[i*n + j];
				if (j != n - 1)
				{
					// East
					Square* eastNeighbour = this->Elements[i*n + j + 1];
					Element2DInterface* interface = new Element2DInterface(element, eastNeighbour);
					this->Interfaces.push_back(interface);
					element->SetEastInterface(interface);
					eastNeighbour->SetWestInterface(interface);
				}
				if (i != n - 1)
				{
					// North
					Square* northNeighbour = this->Elements[(i+1)*n + j];
					Element2DInterface* interface = new Element2DInterface(element, northNeighbour);
					this->Interfaces.push_back(interface);
					element->SetNorthInterface(interface);
					northNeighbour->SetSouthInterface(interface);
				}
			}
		}
	}

	~CartesianGrid2D()
	{
		for (size_t i = 0; i < this->Elements.size(); ++i)
			delete this->Elements[i];
		this->Elements.clear();

		for (size_t i = 0; i < this->Interfaces.size(); ++i)
			delete this->Interfaces[i];
		this->Interfaces.clear();

		for (size_t i = 0; i < this->BoundaryInterfaces.size(); ++i)
			delete this->BoundaryInterfaces[i];
		this->BoundaryInterfaces.clear();
	}

};