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
				RectangularPolygon* rectangle = new RectangularPolygon(i*nx + j, DomPoint(j * hx, i * hy), hx, hy);
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
			RectangularPolygon* rectangle = dynamic_cast<RectangularPolygon*>(this->Elements[j]);
			Edge* southBoundary = new Edge(numberInterface++, rectangle->BottomLeftCorner, hx, rectangle, CartesianShapeOrientation::Horizontal);
			this->Faces.push_back(southBoundary);
			this->BoundaryFaces.push_back(southBoundary);
			dynamic_cast<RectangularPolygon*>(this->Elements[j])->AddSouthFace(southBoundary);

			// North boundary
			rectangle = dynamic_cast<RectangularPolygon*>(this->Elements[(ny - 1)*nx + j]);
			Edge* northBoundary = new Edge(numberInterface++, rectangle->TopLeftCorner, hx, rectangle, CartesianShapeOrientation::Horizontal);
			this->Faces.push_back(northBoundary);
			this->BoundaryFaces.push_back(northBoundary);
			dynamic_cast<RectangularPolygon*>(this->Elements[(ny - 1)*nx + j])->AddNorthFace(northBoundary);
		}

		for (BigNumber i = 0; i < ny; ++i)
		{
			// West boundary
			RectangularPolygon* rectangle = dynamic_cast<RectangularPolygon*>(this->Elements[i*nx]);
			Edge* westBoundary = new Edge(numberInterface++, rectangle->BottomLeftCorner, hy, rectangle, CartesianShapeOrientation::Vertical);
			this->Faces.push_back(westBoundary);
			this->BoundaryFaces.push_back(westBoundary);
			dynamic_cast<RectangularPolygon*>(this->Elements[i*nx])->AddWestFace(westBoundary);

			// East boundary
			rectangle = dynamic_cast<RectangularPolygon*>(this->Elements[i*nx + nx - 1]);
			Edge* eastBoundary = new Edge(numberInterface++, rectangle->BottomRightCorner, hy, rectangle, CartesianShapeOrientation::Vertical);
			this->Faces.push_back(eastBoundary);
			this->BoundaryFaces.push_back(eastBoundary);
			dynamic_cast<RectangularPolygon*>(this->Elements[i*nx + nx - 1])->AddEastFace(eastBoundary);
		}

		for (BigNumber i = 0; i < ny; i++)
		{
			for (BigNumber j = 0; j < nx; j++)
			{
				RectangularPolygon* element = dynamic_cast<RectangularPolygon*>(this->Elements[i*nx + j]);
				if (j != nx - 1)
				{
					// East
					RectangularPolygon* eastNeighbour = dynamic_cast<RectangularPolygon*>(this->Elements[i*nx + j + 1]);
					Edge* interface = new Edge(numberInterface++, eastNeighbour->BottomLeftCorner, hy, element, eastNeighbour, CartesianShapeOrientation::Vertical);
					this->Faces.push_back(interface);
					this->InteriorFaces.push_back(interface);
					element->AddEastFace(interface);
					eastNeighbour->AddWestFace(interface);
				}
				if (i != ny - 1)
				{
					// North
					RectangularPolygon* northNeighbour = dynamic_cast<RectangularPolygon*>(this->Elements[(i + 1)*nx + j]);
					Edge* interface = new Edge(numberInterface++, northNeighbour->BottomLeftCorner, hx, element, northNeighbour, CartesianShapeOrientation::Horizontal);
					this->Faces.push_back(interface);
					this->InteriorFaces.push_back(interface);
					element->AddNorthFace(interface);
					northNeighbour->AddSouthFace(interface);
				}
			}
		}
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
					RectangularPolygon* coarseElement = new RectangularPolygon(i*nx / 2 + j, bottomLeftElement->BottomLeftCorner, bottomLeftElement->WidthX + bottomRightElement->WidthX, bottomLeftElement->WidthY + topLeftElement->WidthY);
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