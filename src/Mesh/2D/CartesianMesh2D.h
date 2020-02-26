#pragma once
#include <vector>
#include "Rectangle.h"
#include "RectangularPolygon.h"
#include "CartesianPolygonalMesh2D.h"
#include "CartesianEdge.h"
#include "../PolyhedralMesh.h"
using namespace std;

class CartesianMesh2D : public PolyhedralMesh<2>
{
public:
	BigNumber Nx;
	BigNumber Ny;

	CartesianMesh2D(BigNumber nx, BigNumber ny) : PolyhedralMesh()
	{
		// nx = ny falls down to square elements
		this->Nx = nx;
		this->Ny = ny;

		double hx = 1 / (double)nx;
		double hy = 1 / (double)ny;

		//----------//
		// Vertices //
		//----------//

		this->Vertices.reserve((nx + 1) * (ny + 1));
		for (BigNumber iy = 0; iy < ny + 1; ++iy)
		{
			for (BigNumber ix = 0; ix < nx + 1; ++ix)
			{
				BigNumber number = indexV(ix, iy);
				Vertex* vertex = new Vertex(number, ix * hx, iy * hy);
				this->Vertices.push_back(vertex);
			}
		}

		//----------//
		// Elements //
		//----------//

		this->Elements.reserve(nx * ny);
		for (BigNumber iy = 0; iy < ny; ++iy)
		{
			for (BigNumber ix = 0; ix < nx; ++ix)
			{
				BigNumber number = index(ix, iy);
				Vertex* bottomLeftCorner  = Vertices[indexV(ix,     iy    )];
				Vertex* topLeftCorner     = Vertices[indexV(ix,     iy + 1)];
				Vertex* topRightCorner    = Vertices[indexV(ix + 1, iy + 1)];
				Vertex* bottomRightCorner = Vertices[indexV(ix + 1, iy    )];
				Rectangle* rectangle = new Rectangle(number, bottomLeftCorner, topLeftCorner, topRightCorner, bottomRightCorner);
				this->Elements.push_back(rectangle);
			}
		}

		//-------//
		// Faces //
		//-------//

		this->Faces.reserve(nx * (ny + 1) + ny * (nx + 1));
		BigNumber numberInterface = 0;

		for (BigNumber ix = 0; ix < nx; ++ix)
		{
			// South boundary
			Rectangle* rectangle = dynamic_cast<Rectangle*>(this->Elements[index(ix, 0)]);
			CartesianEdge* southBoundary = new CartesianEdge(numberInterface++, rectangle->BottomLeftCorner, rectangle->BottomRightCorner, rectangle, CartesianShapeOrientation::Horizontal);
			this->Faces.push_back(southBoundary);
			this->BoundaryFaces.push_back(southBoundary);
			rectangle->SetSouthInterface(southBoundary);

			// North boundary
			rectangle = dynamic_cast<Rectangle*>(this->Elements[index(ix, ny - 1)]);
			CartesianEdge* northBoundary = new CartesianEdge(numberInterface++, rectangle->TopLeftCorner, rectangle->TopRightCorner, rectangle, CartesianShapeOrientation::Horizontal);
			this->Faces.push_back(northBoundary);
			this->BoundaryFaces.push_back(northBoundary);
			rectangle->SetNorthInterface(northBoundary);
		}

		for (BigNumber iy = 0; iy < ny; ++iy)
		{
			// West boundary
			Rectangle* rectangle = dynamic_cast<Rectangle*>(this->Elements[index(0, iy)]);
			CartesianEdge* westBoundary = new CartesianEdge(numberInterface++, rectangle->BottomLeftCorner, rectangle->TopLeftCorner, rectangle, CartesianShapeOrientation::Vertical);
			this->Faces.push_back(westBoundary);
			this->BoundaryFaces.push_back(westBoundary);
			rectangle->SetWestInterface(westBoundary);

			// East boundary
			rectangle = dynamic_cast<Rectangle*>(this->Elements[index(nx-1, iy)]);
			CartesianEdge* eastBoundary = new CartesianEdge(numberInterface++, rectangle->BottomRightCorner, rectangle->TopRightCorner, rectangle, CartesianShapeOrientation::Vertical);
			this->Faces.push_back(eastBoundary);
			this->BoundaryFaces.push_back(eastBoundary);
			rectangle->SetEastInterface(eastBoundary);
		}

		for (BigNumber iy = 0; iy < ny; iy++)
		{
			for (BigNumber ix = 0; ix < nx; ix++)
			{
				Rectangle* element = dynamic_cast<Rectangle*>(this->Elements[index(ix, iy)]);
				if (ix != nx - 1)
				{
					// East
					Rectangle* eastNeighbour = dynamic_cast<Rectangle*>(this->Elements[index(ix + 1, iy)]);
					CartesianEdge* interface = new CartesianEdge(numberInterface++, eastNeighbour->BottomLeftCorner, eastNeighbour->TopLeftCorner, element, eastNeighbour, CartesianShapeOrientation::Vertical);
					this->Faces.push_back(interface);
					this->InteriorFaces.push_back(interface);
					element->SetEastInterface(interface);
					eastNeighbour->SetWestInterface(interface);
				}
				if (iy != ny - 1)
				{
					// North
					Rectangle* northNeighbour = dynamic_cast<Rectangle*>(this->Elements[index(ix, iy + 1)]);
					CartesianEdge* interface = new CartesianEdge(numberInterface++, northNeighbour->BottomLeftCorner, northNeighbour->BottomRightCorner, element, northNeighbour, CartesianShapeOrientation::Horizontal);
					this->Faces.push_back(interface);
					this->InteriorFaces.push_back(interface);
					element->SetNorthInterface(interface);
					northNeighbour->SetSouthInterface(interface);
				}
			}
		}

	}

private:
	inline BigNumber indexV(BigNumber x, BigNumber y)
	{
		return y * (Nx + 1) + x;
	}
	inline BigNumber index(BigNumber x, BigNumber y)
	{
		return y * Nx + x;
	}

public:
	string Description() override
	{
		return "Cartesian " + to_string(this->Nx) + " x " + to_string(this->Ny);
	}

	string FileNamePart() override
	{
		return "cart-n" + to_string(this->Nx);
	}

	double H() override
	{
		return 1.0 / this->Nx;
	}

	double Regularity() override
	{
		return min((double)this->Nx, (double)this->Ny) / max((double)this->Nx, (double)this->Ny);
	}

	void CoarsenMesh(CoarseningStrategy strategy) override
	{
		if (strategy == CoarseningStrategy::StandardCoarsening)
			StandardCoarsening();
		else if (strategy == CoarseningStrategy::AgglomerationCoarsening)
			CoarsenByAgglomerationAndKeepFineFaces();
		else
			assert(false && "Coarsening strategy not implemented!");
		this->CoarseMesh->SetDiffusionCoefficient(this->_diffusionPartition);
		this->CoarseMesh->SetBoundaryConditions(this->_boundaryConditions);
	}

	void RefineMesh(CoarseningStrategy strategy) override
	{
		Utils::FatalError("Refinement strategy not implemented!");
	}

	void StandardCoarsening()
	{
		BigNumber nx = this->Nx;
		BigNumber ny = this->Ny;

		if (nx % 2 != 0 || ny % 2 != 0)
		{
			cout << "Error: impossible to build coarse mesh. Nx and Ny must be even: Nx = " << nx << ", Ny = " << ny << "." << endl;
		}
		else
		{
			CartesianMesh2D* coarseMesh = new CartesianMesh2D(nx / 2, ny / 2);
			coarseMesh->ComesFrom.CS = CoarseningStrategy::StandardCoarsening;
			coarseMesh->ComesFrom.nFineElementsByCoarseElement = 4;
			coarseMesh->ComesFrom.nFineFacesAddedByCoarseElement = 4;
			coarseMesh->ComesFrom.nFineFacesByKeptCoarseFace = 2;

			for (BigNumber i = 0; i < ny; ++i)
			{
				for (BigNumber j = 0; j < nx; ++j)
				{
					Rectangle* fineElement = dynamic_cast<Rectangle*>(this->Elements[i*nx + j]);
					Rectangle* coarseElement = dynamic_cast<Rectangle*>(coarseMesh->Elements[(i / 2) * coarseMesh->Nx + j / 2]);

					coarseElement->FinerElements.push_back(fineElement);
					fineElement->CoarserElement = coarseElement;
					if (i % 2 == 0 && !fineElement->NorthFace->IsDomainBoundary)
					{
						fineElement->NorthFace->IsRemovedOnCoarserGrid = true;
						coarseElement->FinerFacesRemoved.push_back(fineElement->NorthFace);
						coarseElement->SouthFace->FinerFaces.push_back(fineElement->SouthFace);
						fineElement->SouthFace->CoarseFace = coarseElement->SouthFace;
					}
					if (i == ny - 1)
					{
						coarseElement->NorthFace->FinerFaces.push_back(fineElement->NorthFace);
						fineElement->NorthFace->CoarseFace = coarseElement->NorthFace;
					}

					if (j % 2 == 0 && !fineElement->EastFace->IsDomainBoundary)
					{
						fineElement->EastFace->IsRemovedOnCoarserGrid = true;
						coarseElement->FinerFacesRemoved.push_back(fineElement->EastFace);
						coarseElement->WestFace->FinerFaces.push_back(fineElement->WestFace);
						fineElement->WestFace->CoarseFace = coarseElement->WestFace;
					}
					if (j == nx - 1)
					{
						coarseElement->EastFace->FinerFaces.push_back(fineElement->EastFace);
						fineElement->EastFace->CoarseFace = coarseElement->EastFace;
					}
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
			coarseMesh->ComesFrom.CS = CoarseningStrategy::AgglomerationCoarsening;
			coarseMesh->ComesFrom.nFineElementsByCoarseElement = 4;
			coarseMesh->ComesFrom.nFineFacesAddedByCoarseElement = 4;
			coarseMesh->ComesFrom.nFineFacesByKeptCoarseFace = 1;

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
					RectangularPolygon* coarseElement = new RectangularPolygon(i*nx/2 + j, bottomLeftElement->BottomLeftCorner, topLeftElement->TopLeftCorner, topRightElement->TopRightCorner, bottomRightElement->BottomRightCorner);
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
						bottomWestFace->FinerFaces.push_back(bottomLeftElement->WestFace);
						Face<2>* topWestFace = topLeftElement->WestFace->CreateSameGeometricFace(faceNumber++, coarseElement);
						topWestFace->FinerFaces.push_back(topLeftElement->WestFace);
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
					leftNorthFace->FinerFaces.push_back(topLeftElement->NorthFace);
					Face<2>* rightNorthFace = topRightElement->NorthFace->CreateSameGeometricFace(faceNumber++, coarseElement);
					rightNorthFace->FinerFaces.push_back(topRightElement->NorthFace);
					coarseElement->AddNorthFace(leftNorthFace);
					coarseElement->AddNorthFace(rightNorthFace);
					coarseMesh->AddFace(leftNorthFace);
					coarseMesh->AddFace(rightNorthFace);

					// East faces
					Face<2>* topEastFace = topRightElement->EastFace->CreateSameGeometricFace(faceNumber++, coarseElement);
					topEastFace->FinerFaces.push_back(topRightElement->EastFace);
					Face<2>* bottomEastFace = bottomRightElement->EastFace->CreateSameGeometricFace(faceNumber++, coarseElement);
					bottomEastFace->FinerFaces.push_back(bottomRightElement->EastFace);
					coarseElement->AddEastFace(topEastFace);
					coarseElement->AddEastFace(bottomEastFace);
					coarseMesh->AddFace(topEastFace);
					coarseMesh->AddFace(bottomEastFace);

					// South faces
					if (i == 0)
					{
						Face<2>* rightSouthFace = bottomRightElement->SouthFace->CreateSameGeometricFace(faceNumber++, coarseElement);
						rightSouthFace->FinerFaces.push_back(bottomRightElement->SouthFace);
						Face<2>* leftSouthFace = bottomLeftElement->SouthFace->CreateSameGeometricFace(faceNumber++, coarseElement);
						leftSouthFace->FinerFaces.push_back(bottomLeftElement->SouthFace);
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
					bottomLeftElement->NorthFace->IsRemovedOnCoarserGrid = true;
					topLeftElement->EastFace->IsRemovedOnCoarserGrid = true;
					topRightElement->SouthFace->IsRemovedOnCoarserGrid = true;
					bottomRightElement->WestFace->IsRemovedOnCoarserGrid = true;
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