#pragma once
#include "../PolyhedralMesh.h"
#include "RectangularPolygonalElement.h"
using namespace std;

class Square_CartesianPolygonalMesh : public PolyhedralMesh<2>
{
public:
	BigNumber Nx;
	BigNumber Ny;

	Square_CartesianPolygonalMesh() : PolyhedralMesh() {}

	Square_CartesianPolygonalMesh(BigNumber nx, BigNumber ny) : PolyhedralMesh()
	{
		// nx = ny falls down to square elements
		this->Nx = nx;
		this->Ny = ny;

		double hx = 1.0 / nx;
		double hy = 1.0 / ny;

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
				RectangularPolygonalElement* rectangle = new RectangularPolygonalElement(number, bottomLeftCorner, topLeftCorner, topRightCorner, bottomRightCorner);
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
			RectangularPolygonalElement* rectangle = dynamic_cast<RectangularPolygonalElement*>(this->Elements[index(ix, 0)]);
			CartesianEdge* southBoundary = new CartesianEdge(numberInterface++, rectangle->BottomLeftCorner, rectangle->BottomRightCorner, rectangle, CartesianShapeOrientation::Horizontal);
			this->Faces.push_back(southBoundary);
			this->BoundaryFaces.push_back(southBoundary);
			rectangle->AddSouthFace(southBoundary);

			// North boundary
			rectangle = dynamic_cast<RectangularPolygonalElement*>(this->Elements[index(ix, ny - 1)]);
			CartesianEdge* northBoundary = new CartesianEdge(numberInterface++, rectangle->TopLeftCorner, rectangle->TopRightCorner, rectangle, CartesianShapeOrientation::Horizontal);
			this->Faces.push_back(northBoundary);
			this->BoundaryFaces.push_back(northBoundary);
			rectangle->AddNorthFace(northBoundary);
		}

		for (BigNumber iy = 0; iy < ny; ++iy)
		{
			// West boundary
			RectangularPolygonalElement* rectangle = dynamic_cast<RectangularPolygonalElement*>(this->Elements[index(0, iy)]);
			CartesianEdge* westBoundary = new CartesianEdge(numberInterface++, rectangle->BottomLeftCorner, rectangle->TopLeftCorner, rectangle, CartesianShapeOrientation::Vertical);
			this->Faces.push_back(westBoundary);
			this->BoundaryFaces.push_back(westBoundary);
			rectangle->AddWestFace(westBoundary);

			// East boundary
			rectangle = dynamic_cast<RectangularPolygonalElement*>(this->Elements[index(nx - 1, iy)]);
			CartesianEdge* eastBoundary = new CartesianEdge(numberInterface++, rectangle->BottomRightCorner, rectangle->TopRightCorner, rectangle, CartesianShapeOrientation::Vertical);
			this->Faces.push_back(eastBoundary);
			this->BoundaryFaces.push_back(eastBoundary);
			rectangle->AddEastFace(eastBoundary);
		}

		for (BigNumber iy = 0; iy < ny; iy++)
		{
			for (BigNumber ix = 0; ix < nx; ix++)
			{
				RectangularPolygonalElement* element = dynamic_cast<RectangularPolygonalElement*>(this->Elements[index(ix, iy)]);
				if (ix != nx - 1)
				{
					// East
					RectangularPolygonalElement* eastNeighbour = dynamic_cast<RectangularPolygonalElement*>(this->Elements[index(ix + 1, iy)]);
					CartesianEdge* interface = new CartesianEdge(numberInterface++, eastNeighbour->BottomLeftCorner, eastNeighbour->TopLeftCorner, element, eastNeighbour, CartesianShapeOrientation::Vertical);
					this->Faces.push_back(interface);
					this->InteriorFaces.push_back(interface);
					element->AddEastFace(interface);
					eastNeighbour->AddWestFace(interface);
				}
				if (iy != ny - 1)
				{
					// North
					RectangularPolygonalElement* northNeighbour = dynamic_cast<RectangularPolygonalElement*>(this->Elements[index(ix, iy + 1)]);
					CartesianEdge* interface = new CartesianEdge(numberInterface++, northNeighbour->BottomLeftCorner, northNeighbour->BottomRightCorner, element, northNeighbour, CartesianShapeOrientation::Horizontal);
					this->Faces.push_back(interface);
					this->InteriorFaces.push_back(interface);
					element->AddNorthFace(interface);
					northNeighbour->AddSouthFace(interface);
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
		return "CartesianPolygonalMesh: " + to_string(this->Nx) + " x " + to_string(this->Ny);
	}

	string FileNamePart() override
	{
		return "square-cartpoly-n" + to_string(this->Nx);
	}

	string GeometryDescription() override
	{
		return "Square";
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
		if (strategy == CoarseningStrategy::AgglomerationCoarsening)
			CoarsenByAgglomerationAndKeepFineFaces();
		else
			PolyhedralMesh<2>::CoarsenMesh(strategy);

		this->CoarseMesh->SetDiffusionCoefficient(this->_diffusionPartition);
		this->CoarseMesh->SetBoundaryConditions(this->_boundaryConditions);
	}

	void CoarsenByAgglomerationAndKeepFineFaces()
	{
		BigNumber nx = this->Nx;
		BigNumber ny = this->Ny;

		if (nx % 2 != 0 || ny % 2 != 0)
		{
			cout << "Error: impossible to build coarse mesh. Nx and Ny must be even: Nx = " << nx << ", Ny = " << ny << "." << endl;
			return;
		}

		Square_CartesianPolygonalMesh* coarseMesh = new Square_CartesianPolygonalMesh();
		coarseMesh->Nx = nx / 2;
		coarseMesh->Ny = ny / 2;
		coarseMesh->ComesFrom.CS = CoarseningStrategy::AgglomerationCoarsening;
		coarseMesh->ComesFrom.nFineElementsByCoarseElement = 4;
		coarseMesh->ComesFrom.nFineFacesAddedByCoarseElement = 4;
		coarseMesh->ComesFrom.nFineFacesByKeptCoarseFace = 1;

		BigNumber faceNumber = 0;
		for (BigNumber i = 0; i < ny / 2; ++i)
		{
			for (BigNumber j = 0; j < nx / 2; ++j)
			{
				RectangularPolygonalElement* bottomLeftElement = dynamic_cast<RectangularPolygonalElement*>(this->Elements[2 * i * nx + 2 * j]);
				RectangularPolygonalElement* bottomRightElement = dynamic_cast<RectangularPolygonalElement*>(this->Elements[2 * i * nx + 2 * j + 1]);
				RectangularPolygonalElement* topLeftElement = dynamic_cast<RectangularPolygonalElement*>(this->Elements[(2 * i + 1) * nx + 2 * j]);
				RectangularPolygonalElement* topRightElement = dynamic_cast<RectangularPolygonalElement*>(this->Elements[(2 * i + 1) * nx + 2 * j + 1]);

				// Coarse element
				RectangularPolygonalElement* coarseElement = new RectangularPolygonalElement(i*nx / 2 + j, bottomLeftElement->BottomLeftCorner, topLeftElement->TopLeftCorner, topRightElement->TopRightCorner, bottomRightElement->BottomRightCorner);
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
					for (Face<2>* f : bottomLeftElement->WestFaces)
					{
						Face<2>* newFace = f->CreateSameGeometricFace(faceNumber++, coarseElement);
						f->CoarseFace = newFace;
						newFace->FinerFaces.push_back(f);
						coarseElement->AddWestFace(newFace);
						coarseMesh->AddFace(newFace);
					}
					for (Face<2>* f : topLeftElement->WestFaces)
					{
						Face<2>* newFace = f->CreateSameGeometricFace(faceNumber++, coarseElement);
						f->CoarseFace = newFace;
						newFace->FinerFaces.push_back(f);
						coarseElement->AddWestFace(newFace);
						coarseMesh->AddFace(newFace);
					}
				}
				else
				{
					RectangularPolygonalElement* leftNeighbour = dynamic_cast<RectangularPolygonalElement*>(coarseMesh->Elements[i * nx / 2 + j - 1]);
					coarseElement->SetWestFacesFromNeighbour(leftNeighbour);
				}

				// North faces
				for (Face<2>* f : topLeftElement->NorthFaces)
				{
					Face<2>* newFace = f->CreateSameGeometricFace(faceNumber++, coarseElement);
					f->CoarseFace = newFace;
					newFace->FinerFaces.push_back(f);
					coarseElement->AddNorthFace(newFace);
					coarseMesh->AddFace(newFace);
				}
				for (Face<2>* f : topRightElement->NorthFaces)
				{
					Face<2>* newFace = f->CreateSameGeometricFace(faceNumber++, coarseElement);
					f->CoarseFace = newFace;
					newFace->FinerFaces.push_back(f);
					coarseElement->AddNorthFace(newFace);
					coarseMesh->AddFace(newFace);
				}

				// East faces
				for (Face<2>* f : topRightElement->EastFaces)
				{
					Face<2>* newFace = f->CreateSameGeometricFace(faceNumber++, coarseElement);
					f->CoarseFace = newFace;
					newFace->FinerFaces.push_back(f);
					coarseElement->AddEastFace(newFace);
					coarseMesh->AddFace(newFace);
				}
				for (Face<2>* f : bottomRightElement->EastFaces)
				{
					Face<2>* newFace = f->CreateSameGeometricFace(faceNumber++, coarseElement);
					f->CoarseFace = newFace;
					newFace->FinerFaces.push_back(f);
					coarseElement->AddEastFace(newFace);
					coarseMesh->AddFace(newFace);
				}

				// South faces
				if (i == 0)
				{
					for (Face<2>* f : bottomRightElement->SouthFaces)
					{
						Face<2>* newFace = f->CreateSameGeometricFace(faceNumber++, coarseElement);
						f->CoarseFace = newFace;
						newFace->FinerFaces.push_back(f);
						coarseElement->AddSouthFace(newFace);
						coarseMesh->AddFace(newFace);
					}
					for (Face<2>* f : bottomLeftElement->SouthFaces)
					{
						Face<2>* newFace = f->CreateSameGeometricFace(faceNumber++, coarseElement);
						f->CoarseFace = newFace;
						newFace->FinerFaces.push_back(f);
						coarseElement->AddSouthFace(newFace);
						coarseMesh->AddFace(newFace);
					}
				}
				else
				{
					RectangularPolygonalElement* bottomNeighbour = dynamic_cast<RectangularPolygonalElement*>(coarseMesh->Elements[(i - 1) * nx / 2 + j]);
					coarseElement->SetSouthFacesFromNeighbour(bottomNeighbour);
				}

				// Finer faces removed
				for (Face<2>* f : bottomLeftElement->NorthFaces)
				{
					f->IsRemovedOnCoarserGrid = true;
					coarseElement->FinerFacesRemoved.push_back(f);
				}
				for (Face<2>* f : topLeftElement->EastFaces)
				{
					f->IsRemovedOnCoarserGrid = true;
					coarseElement->FinerFacesRemoved.push_back(f);
				}
				for (Face<2>* f : topRightElement->SouthFaces)
				{
					f->IsRemovedOnCoarserGrid = true;
					coarseElement->FinerFacesRemoved.push_back(f);
				}
				for (Face<2>* f : bottomRightElement->WestFaces)
				{
					f->IsRemovedOnCoarserGrid = true;
					coarseElement->FinerFacesRemoved.push_back(f);
				}
			}
		}

		this->CoarseMesh = coarseMesh;
	}

};