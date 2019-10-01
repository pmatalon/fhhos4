#pragma once
#include "../Mesh.h"
#include "RectangularPolygon.h"
using namespace std;

class CartesianPolygonalMesh2D : public Mesh<2>
{
public:
	BigNumber Nx;
	BigNumber Ny;

	CartesianPolygonalMesh2D() : Mesh() {}

	CartesianPolygonalMesh2D(BigNumber nx, BigNumber ny) : Mesh()
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
				RectangularPolygon* rectangle = new RectangularPolygon(number, bottomLeftCorner, topLeftCorner, topRightCorner, bottomRightCorner);
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
			RectangularPolygon* rectangle = dynamic_cast<RectangularPolygon*>(this->Elements[index(ix, 0)]);
			Edge* southBoundary = new Edge(numberInterface++, rectangle->BottomLeftCorner, rectangle->BottomRightCorner, rectangle, CartesianShapeOrientation::Horizontal);
			this->Faces.push_back(southBoundary);
			this->BoundaryFaces.push_back(southBoundary);
			rectangle->AddSouthFace(southBoundary);

			// North boundary
			rectangle = dynamic_cast<RectangularPolygon*>(this->Elements[index(ix, ny - 1)]);
			Edge* northBoundary = new Edge(numberInterface++, rectangle->TopLeftCorner, rectangle->TopRightCorner, rectangle, CartesianShapeOrientation::Horizontal);
			this->Faces.push_back(northBoundary);
			this->BoundaryFaces.push_back(northBoundary);
			rectangle->AddNorthFace(northBoundary);
		}

		for (BigNumber iy = 0; iy < ny; ++iy)
		{
			// West boundary
			RectangularPolygon* rectangle = dynamic_cast<RectangularPolygon*>(this->Elements[index(0, iy)]);
			Edge* westBoundary = new Edge(numberInterface++, rectangle->BottomLeftCorner, rectangle->TopLeftCorner, rectangle, CartesianShapeOrientation::Vertical);
			this->Faces.push_back(westBoundary);
			this->BoundaryFaces.push_back(westBoundary);
			rectangle->AddWestFace(westBoundary);

			// East boundary
			rectangle = dynamic_cast<RectangularPolygon*>(this->Elements[index(nx - 1, iy)]);
			Edge* eastBoundary = new Edge(numberInterface++, rectangle->BottomRightCorner, rectangle->TopRightCorner, rectangle, CartesianShapeOrientation::Vertical);
			this->Faces.push_back(eastBoundary);
			this->BoundaryFaces.push_back(eastBoundary);
			rectangle->AddEastFace(eastBoundary);
		}

		for (BigNumber iy = 0; iy < ny; iy++)
		{
			for (BigNumber ix = 0; ix < nx; ix++)
			{
				RectangularPolygon* element = dynamic_cast<RectangularPolygon*>(this->Elements[index(ix, iy)]);
				if (ix != nx - 1)
				{
					// East
					RectangularPolygon* eastNeighbour = dynamic_cast<RectangularPolygon*>(this->Elements[index(ix + 1, iy)]);
					Edge* interface = new Edge(numberInterface++, eastNeighbour->BottomLeftCorner, eastNeighbour->TopLeftCorner, element, eastNeighbour, CartesianShapeOrientation::Vertical);
					this->Faces.push_back(interface);
					this->InteriorFaces.push_back(interface);
					element->AddEastFace(interface);
					eastNeighbour->AddWestFace(interface);
				}
				if (iy != ny - 1)
				{
					// North
					RectangularPolygon* northNeighbour = dynamic_cast<RectangularPolygon*>(this->Elements[index(ix, iy + 1)]);
					Edge* interface = new Edge(numberInterface++, northNeighbour->BottomLeftCorner, northNeighbour->BottomRightCorner, element, northNeighbour, CartesianShapeOrientation::Horizontal);
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
	string Description()
	{
		return "CartesianPolygonalMesh: Subdivisions in each cartesian direction: " + to_string(this->Nx) + " x " + to_string(this->Ny);
	}

	string FileNamePart()
	{
		return "n" + to_string(this->Nx);
	}

	double H()
	{
		return 1 / (double)this->Nx;
	}

	virtual void CoarsenMesh(CoarseningStrategy strategy)
	{
		/*if (strategy == CoarseningStrategy::Standard)
			CoarsenByAgglomerationAndMergeColinearFaces();
		else*/ if (strategy == CoarseningStrategy::Agglomeration)
			CoarsenByAgglomerationAndKeepFineFaces();
		else
			assert(false && "Coarsening strategy not implemented!");

		this->CoarseMesh->SetDiffusionCoefficient(this->_diffusionPartition);
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

			BigNumber faceNumber = 0;
			for (BigNumber i = 0; i < ny / 2; ++i)
			{
				for (BigNumber j = 0; j < nx / 2; ++j)
				{
					RectangularPolygon* bottomLeftElement = dynamic_cast<RectangularPolygon*>(this->Elements[2 * i * nx + 2 * j]);
					RectangularPolygon* bottomRightElement = dynamic_cast<RectangularPolygon*>(this->Elements[2 * i * nx + 2 * j + 1]);
					RectangularPolygon* topLeftElement = dynamic_cast<RectangularPolygon*>(this->Elements[(2 * i + 1) * nx + 2 * j]);
					RectangularPolygon* topRightElement = dynamic_cast<RectangularPolygon*>(this->Elements[(2 * i + 1) * nx + 2 * j + 1]);

					// Coarse element
					RectangularPolygon* coarseElement = new RectangularPolygon(i*nx / 2 + j, bottomLeftElement->BottomLeftCorner, topLeftElement->TopLeftCorner, topRightElement->TopRightCorner, bottomRightElement->BottomRightCorner);
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
							newFace->FinerFaces.push_back(f);
							coarseElement->AddWestFace(newFace);
							coarseMesh->AddFace(newFace);
						}
						for (Face<2>* f : topLeftElement->WestFaces)
						{
							Face<2>* newFace = f->CreateSameGeometricFace(faceNumber++, coarseElement);
							newFace->FinerFaces.push_back(f);
							coarseElement->AddWestFace(newFace);
							coarseMesh->AddFace(newFace);
						}
					}
					else
					{
						RectangularPolygon* leftNeighbour = dynamic_cast<RectangularPolygon*>(coarseMesh->Elements[i * nx / 2 + j - 1]);
						coarseElement->SetWestFacesFromNeighbour(leftNeighbour);
					}

					// North faces
					for (Face<2>* f : topLeftElement->NorthFaces)
					{
						Face<2>* newFace = f->CreateSameGeometricFace(faceNumber++, coarseElement);
						newFace->FinerFaces.push_back(f);
						coarseElement->AddNorthFace(newFace);
						coarseMesh->AddFace(newFace);
					}
					for (Face<2>* f : topRightElement->NorthFaces)
					{
						Face<2>* newFace = f->CreateSameGeometricFace(faceNumber++, coarseElement);
						newFace->FinerFaces.push_back(f);
						coarseElement->AddNorthFace(newFace);
						coarseMesh->AddFace(newFace);
					}

					// East faces
					for (Face<2>* f : topRightElement->EastFaces)
					{
						Face<2>* newFace = f->CreateSameGeometricFace(faceNumber++, coarseElement);
						newFace->FinerFaces.push_back(f);
						coarseElement->AddEastFace(newFace);
						coarseMesh->AddFace(newFace);
					}
					for (Face<2>* f : bottomRightElement->EastFaces)
					{
						Face<2>* newFace = f->CreateSameGeometricFace(faceNumber++, coarseElement);
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
							newFace->FinerFaces.push_back(f);
							coarseElement->AddSouthFace(newFace);
							coarseMesh->AddFace(newFace);
						}
						for (Face<2>* f : bottomLeftElement->SouthFaces)
						{
							Face<2>* newFace = f->CreateSameGeometricFace(faceNumber++, coarseElement);
							newFace->FinerFaces.push_back(f);
							coarseElement->AddSouthFace(newFace);
							coarseMesh->AddFace(newFace);
						}
					}
					else
					{
						RectangularPolygon* bottomNeighbour = dynamic_cast<RectangularPolygon*>(coarseMesh->Elements[(i - 1) * nx / 2 + j]);
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
	}

};