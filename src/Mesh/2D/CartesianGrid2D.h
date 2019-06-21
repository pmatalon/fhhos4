#pragma once
#include <vector>
#include "Rectangle.h"
#include "RectangularPolygon.h"
#include "CartesianPolygonalMesh2D.h"
#include "Edge.h"
#include "../Mesh.h"
using namespace std;

class CartesianGrid2D : public Mesh<2>
{
public:
	BigNumber Nx;
	BigNumber Ny;

	CartesianGrid2D(BigNumber nx, BigNumber ny) : Mesh()
	{
		// nx = ny falls down to square elements
		this->Nx = nx;
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
			Edge* southBoundary = new Edge(numberInterface++, rectangle->BottomLeftCorner, hx, rectangle, CartesianShapeOrientation::Horizontal);
			this->Faces.push_back(southBoundary);
			this->BoundaryFaces.push_back(southBoundary);
			dynamic_cast<Rectangle*>(this->Elements[j])->SetSouthInterface(southBoundary);

			// North boundary
			rectangle = dynamic_cast<Rectangle*>(this->Elements[(ny - 1)*nx + j]);
			Edge* northBoundary = new Edge(numberInterface++, rectangle->TopLeftCorner, hx, rectangle, CartesianShapeOrientation::Horizontal);
			this->Faces.push_back(northBoundary);
			this->BoundaryFaces.push_back(northBoundary);
			dynamic_cast<Rectangle*>(this->Elements[(ny - 1)*nx + j])->SetNorthInterface(northBoundary);
		}

		for (BigNumber i = 0; i < ny; ++i)
		{
			// West boundary
			Rectangle* rectangle = dynamic_cast<Rectangle*>(this->Elements[i*nx]);
			Edge* westBoundary = new Edge(numberInterface++, rectangle->BottomLeftCorner, hy, rectangle, CartesianShapeOrientation::Vertical);
			this->Faces.push_back(westBoundary);
			this->BoundaryFaces.push_back(westBoundary);
			dynamic_cast<Rectangle*>(this->Elements[i*nx])->SetWestInterface(westBoundary);

			// East boundary
			rectangle = dynamic_cast<Rectangle*>(this->Elements[i*nx + nx - 1]);
			Edge* eastBoundary = new Edge(numberInterface++, rectangle->BottomRightCorner, hy, rectangle, CartesianShapeOrientation::Vertical);
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
					Edge* interface = new Edge(numberInterface++, eastNeighbour->BottomLeftCorner, hy, element, eastNeighbour, CartesianShapeOrientation::Vertical);
					this->Faces.push_back(interface);
					this->InteriorFaces.push_back(interface);
					element->SetEastInterface(interface);
					eastNeighbour->SetWestInterface(interface);
				}
				if (i != ny - 1)
				{
					// North
					Rectangle* northNeighbour = dynamic_cast<Rectangle*>(this->Elements[(i + 1)*nx + j]);
					Edge* interface = new Edge(numberInterface++, northNeighbour->BottomLeftCorner, hx, element, northNeighbour, CartesianShapeOrientation::Horizontal);
					this->Faces.push_back(interface);
					this->InteriorFaces.push_back(interface);
					element->SetNorthInterface(interface);
					northNeighbour->SetSouthInterface(interface);
				}
			}
		}

	}

	string Description()
	{
		return "Subdivisions in each cartesian direction: " + to_string(this->Nx) + " x " + to_string(this->Ny);
	}

	string FileNamePart()
	{
		return "n" + to_string(this->Nx);
	}

	double H()
	{
		return 1 / (double)this->Nx;
	}

	void CoarsenMesh(CoarseningStrategy strategy)
	{
		if (strategy == CoarseningStrategy::AgglomerationAndMergeColinearFaces)
			CoarsenByAgglomerationAndMergeColinearFaces();
		else if (strategy == CoarseningStrategy::AgglomerationAndKeepFineFaces)
			CoarsenByAgglomerationAndKeepFineFaces();
		else
			assert(false && "Coarsening strategy not implemented!");
	}

	void CoarsenByAgglomerationAndMergeColinearFaces()
	{
		BigNumber nx = this->Nx;
		BigNumber ny = this->Ny;

		if (nx % 2 != 0 || ny % 2 != 0)
		{
			cout << "Error: impossible to build coarse mesh. Nx and Ny must be even: Nx = " << nx << ", Ny = " << ny << "." << endl;
		}
		else
		{
			CartesianGrid2D* coarseMesh = new CartesianGrid2D(nx / 2, ny / 2);

			for (BigNumber i = 0; i < ny; ++i)
			{
				for (BigNumber j = 0; j < nx; ++j)
				{
					Rectangle* fineElement = dynamic_cast<Rectangle*>(this->Elements[i*nx + j]);
					auto coarseElement = coarseMesh->Elements[(i / 2) * coarseMesh->Nx + j / 2];

					coarseElement->FinerElements.push_back(fineElement);
					fineElement->CoarserElement = coarseElement;
					if (i % 2 == 0 && !fineElement->NorthFace->IsDomainBoundary)
						coarseElement->FinerFacesRemoved.push_back(fineElement->NorthFace);
					if (j % 2 == 0 && !fineElement->EastFace->IsDomainBoundary)
						coarseElement->FinerFacesRemoved.push_back(fineElement->EastFace);
				}
			}

			/*for (auto e : coarseMesh->Elements)
			{
				cout << "coarseElement " << e->Number << " removed faces: " ;
				for (auto f : e->FinerFacesRemoved)
					cout << f->Number << " ";
				cout << endl;
			}*/

			this->CoarseMesh = coarseMesh;
		}
	}

	void CoarsenByAgglomerationAndKeepFineFaces()
	{
		BigNumber nx = this->Nx;
		BigNumber ny = this->Ny;

		if (nx % 2 != 0 || ny % 2 != 0)
		{
			cout << "Error: impossible to build coarse mesh. Nx and Ny must be even: Nx = " << nx << ", Ny = " << ny << "." << endl;
		}
		else
		{
			CartesianPolygonalMesh2D* coarseMesh = new CartesianPolygonalMesh2D();
			coarseMesh->Nx = nx / 2;
			coarseMesh->Ny = ny / 2;

			// Elements //
			BigNumber faceNumber = 0;
			for (BigNumber i = 0; i < ny/2; ++i)
			{
				for (BigNumber j = 0; j < nx/2; ++j)
				{
					Rectangle* bottomLeftElement = dynamic_cast<Rectangle*>(this->Elements[2 * i * nx + 2 * j]);
					Rectangle* bottomRightElement = dynamic_cast<Rectangle*>(this->Elements[2 * i * nx + 2 * j + 1]);
					Rectangle* topLeftElement = dynamic_cast<Rectangle*>(this->Elements[(2 * i + 1) * nx + 2 * j]);
					Rectangle* topRightElement = dynamic_cast<Rectangle*>(this->Elements[(2 * i + 1) * nx + 2 * j + 1]);

					// Coarse element
					RectangularPolygon* coarseElement = new RectangularPolygon(i*nx/2 + j, bottomLeftElement->BottomLeftCorner, bottomLeftElement->WidthX + bottomRightElement->WidthX, bottomLeftElement->WidthY + topLeftElement->WidthY);
					coarseMesh->Elements.push_back(coarseElement);

					coarseElement->FinerElements.push_back(bottomLeftElement);
					coarseElement->FinerElements.push_back(bottomRightElement);
					coarseElement->FinerElements.push_back(topLeftElement);
					coarseElement->FinerElements.push_back(topRightElement);

					bottomLeftElement->CoarserElement = coarseElement;
					bottomRightElement->CoarserElement = coarseElement;
					topLeftElement->CoarserElement = coarseElement;
					topRightElement->CoarserElement = coarseElement;

					// West faces
					if (j == 0)
					{
						Face<2>* bottomWestFace = bottomLeftElement->WestFace->CreateSameGeometricFace(faceNumber++, coarseElement);
						Face<2>* topWestFace = topLeftElement->WestFace->CreateSameGeometricFace(faceNumber++, coarseElement);
						coarseElement->AddWestFace(bottomWestFace);
						coarseElement->AddWestFace(topWestFace);
						coarseMesh->AddFace(bottomWestFace);
						coarseMesh->AddFace(topWestFace);
					}
					else
					{
						RectangularPolygon* leftNeighbour = dynamic_cast<RectangularPolygon*>(coarseMesh->Elements[i * nx/2 + j - 1]);
						coarseElement->SetWestFacesFromNeighbour(leftNeighbour);
					}

					// North faces
					Face<2>* leftNorthFace = topLeftElement->NorthFace->CreateSameGeometricFace(faceNumber++, coarseElement);
					Face<2>* rightNorthFace = topRightElement->NorthFace->CreateSameGeometricFace(faceNumber++, coarseElement);
					coarseElement->AddNorthFace(leftNorthFace);
					coarseElement->AddNorthFace(rightNorthFace);
					coarseMesh->AddFace(leftNorthFace);
					coarseMesh->AddFace(rightNorthFace);

					// East faces
					Face<2>* topEastFace = topRightElement->EastFace->CreateSameGeometricFace(faceNumber++, coarseElement);
					Face<2>* bottomEastFace = bottomRightElement->EastFace->CreateSameGeometricFace(faceNumber++, coarseElement);
					coarseElement->AddEastFace(topEastFace);
					coarseElement->AddEastFace(bottomEastFace);
					coarseMesh->AddFace(topEastFace);
					coarseMesh->AddFace(bottomEastFace);

					// South faces
					if (i == 0)
					{
						Face<2>* rightSouthFace = bottomRightElement->SouthFace->CreateSameGeometricFace(faceNumber++, coarseElement);
						Face<2>* leftSouthFace = bottomLeftElement->SouthFace->CreateSameGeometricFace(faceNumber++, coarseElement);
						coarseElement->AddSouthFace(rightSouthFace);
						coarseElement->AddSouthFace(leftSouthFace);
						coarseMesh->AddFace(rightSouthFace);
						coarseMesh->AddFace(leftSouthFace);
					}
					else
					{
						RectangularPolygon* bottomNeighbour = dynamic_cast<RectangularPolygon*>(coarseMesh->Elements[(i - 1) * nx/2 + j]);
						coarseElement->SetSouthFacesFromNeighbour(bottomNeighbour);
					}

					// Finer faces removed
					coarseElement->FinerFacesRemoved.push_back(bottomLeftElement->NorthFace);
					coarseElement->FinerFacesRemoved.push_back(topLeftElement->EastFace);
					coarseElement->FinerFacesRemoved.push_back(topRightElement->SouthFace);
					coarseElement->FinerFacesRemoved.push_back(bottomRightElement->WestFace);
				}
			}

			this->CoarseMesh = coarseMesh;
		}
	}

};