#pragma once
#include <vector>
#include "Element.h"
#include "Square.h"
#include "IntervalFace.h"
#include "Mesh.h"

using namespace std;

class CartesianGrid2D : public Mesh<2>
{
public:
	CartesianGrid2D(BigNumber n) : Mesh(n)
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

		//-------//
		// Faces //
		//-------//

		this->Faces.reserve(n * (n + 3));
		BigNumber numberInterface = 0;

		for (BigNumber j = 0; j < n; ++j)
		{
			// South boundary
			Square* square = dynamic_cast<Square*>(this->Elements[j]);
			IntervalFace* southBoundary = new IntervalFace(numberInterface++, square->BottomLeftCorner, h, square, CartesianShapeOrientation::Horizontal);
			this->Faces.push_back(southBoundary);
			this->BoundaryFaces.push_back(southBoundary);
			dynamic_cast<Square*>(this->Elements[j])->SetSouthInterface(southBoundary);

			// North boundary
			square = dynamic_cast<Square*>(this->Elements[(n - 1)*n + j]);
			IntervalFace* northBoundary = new IntervalFace(numberInterface++, square->TopLeftCorner, h, square, CartesianShapeOrientation::Horizontal);
			this->Faces.push_back(northBoundary);
			this->BoundaryFaces.push_back(northBoundary);
			dynamic_cast<Square*>(this->Elements[(n - 1)*n + j])->SetNorthInterface(northBoundary);
		}

		for (BigNumber i = 0; i < n; ++i)
		{
			// West boundary
			Square* square = dynamic_cast<Square*>(this->Elements[i*n]);
			IntervalFace* westBoundary = new IntervalFace(numberInterface++, square->BottomLeftCorner, h, square, CartesianShapeOrientation::Vertical);
			this->Faces.push_back(westBoundary);
			this->BoundaryFaces.push_back(westBoundary);
			dynamic_cast<Square*>(this->Elements[i*n])->SetWestInterface(westBoundary);

			// East boundary
			square = dynamic_cast<Square*>(this->Elements[i*n + n - 1]);
			IntervalFace* eastBoundary = new IntervalFace(numberInterface++, square->BottomRightCorner, h, square, CartesianShapeOrientation::Vertical);
			this->Faces.push_back(eastBoundary);
			this->BoundaryFaces.push_back(eastBoundary);
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
					IntervalFace* interface = new IntervalFace(numberInterface++, eastNeighbour->BottomLeftCorner, h, element, eastNeighbour, CartesianShapeOrientation::Vertical);
					this->Faces.push_back(interface);
					this->InteriorFaces.push_back(interface);
					element->SetEastInterface(interface);
					eastNeighbour->SetWestInterface(interface);
				}
				if (i != n - 1)
				{
					// North
					Square* northNeighbour = dynamic_cast<Square*>(this->Elements[(i+1)*n + j]);
					IntervalFace* interface = new IntervalFace(numberInterface++, northNeighbour->BottomLeftCorner, h, element, northNeighbour, CartesianShapeOrientation::Horizontal);
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
		BigNumber n = this->N;
		if (n % 2 != 0)
		{
			cout << "Error: impossible to build coarser mesh." << endl;
		}
		else
		{
			CartesianGrid2D* coarserMesh = new CartesianGrid2D(n / 2);

			for (BigNumber i = 0; i < n; ++i)
			{
				for (BigNumber j = 0; j < n; ++j)
				{
					Square* fineElement = dynamic_cast<Square*>(this->Elements[i*n + j]);
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