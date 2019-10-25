#pragma once
#include <vector>
#include "Triangle.h"
#include "../Mesh.h"
using namespace std;

class TriangularMesh : public Mesh<2>
{
public:
	BigNumber Nx;
	BigNumber Ny;

	TriangularMesh(BigNumber nx, BigNumber ny) : Mesh()
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
				Triangle* triangle = new Triangle(number, bottomLeftCorner, bottomRightCorner, topLeftCorner);
				this->Elements.push_back(triangle);
				Triangle* triangle2 = new Triangle(number + 1, topLeftCorner, bottomRightCorner, topRightCorner);
				this->Elements.push_back(triangle2);
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
			Triangle* triangle = dynamic_cast<Triangle*>(this->Elements[index(ix, 0)]);
			Edge* southBoundary = new Edge(numberInterface++, triangle->V1, triangle->V2, triangle);
			this->Faces.push_back(southBoundary);
			this->BoundaryFaces.push_back(southBoundary);
			triangle->AddFace(southBoundary);

			// North boundary
			triangle = dynamic_cast<Triangle*>(this->Elements[index(ix, ny - 1) + 1]);
			Edge* northBoundary = new Edge(numberInterface++, triangle->V1, triangle->V3, triangle);
			this->Faces.push_back(northBoundary);
			this->BoundaryFaces.push_back(northBoundary);
			triangle->AddFace(northBoundary);
		}

		for (BigNumber iy = 0; iy < ny; ++iy)
		{
			// West boundary
			Triangle* triangle = dynamic_cast<Triangle*>(this->Elements[index(0, iy)]);
			Edge* westBoundary = new Edge(numberInterface++, triangle->V3, triangle->V1, triangle);
			this->Faces.push_back(westBoundary);
			this->BoundaryFaces.push_back(westBoundary);
			triangle->AddFace(westBoundary);

			// East boundary
			triangle = dynamic_cast<Triangle*>(this->Elements[index(nx-1, iy)+1]);
			Edge* eastBoundary = new Edge(numberInterface++, triangle->V2, triangle->V3, triangle);
			this->Faces.push_back(eastBoundary);
			this->BoundaryFaces.push_back(eastBoundary);
			triangle->AddFace(eastBoundary);
		}

		for (BigNumber iy = 0; iy < ny; iy++)
		{
			for (BigNumber ix = 0; ix < nx; ix++)
			{
				Triangle* element1 = dynamic_cast<Triangle*>(this->Elements[index(ix, iy)]);
				Triangle* element2 = dynamic_cast<Triangle*>(this->Elements[index(ix, iy)+1]);

				Edge* interface = new Edge(numberInterface++, element1->V3, element1->V2, element1, element2);
				this->Faces.push_back(interface);
				this->InteriorFaces.push_back(interface);
				element1->AddFace(interface);
				element2->AddFace(interface);

				if (ix != nx - 1)
				{
					// East
					Triangle* eastNeighbour = dynamic_cast<Triangle*>(this->Elements[index(ix + 1, iy)]);
					interface = new Edge(numberInterface++, eastNeighbour->V3, eastNeighbour->V1, element2, eastNeighbour);
					this->Faces.push_back(interface);
					this->InteriorFaces.push_back(interface);
					element2->AddFace(interface);
					eastNeighbour->AddFace(interface);
				}
				if (iy != ny - 1)
				{
					// North
					Triangle* northNeighbour = dynamic_cast<Triangle*>(this->Elements[index(ix, iy + 1)]);
					interface = new Edge(numberInterface++, northNeighbour->V1, northNeighbour->V2, element2, northNeighbour);
					this->Faces.push_back(interface);
					this->InteriorFaces.push_back(interface);
					element2->AddFace(interface);
					northNeighbour->AddFace(interface);
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
		return (y * Nx + x) * 2;
	}

public:
	string Description()
	{
		return "Triangular " + to_string(this->Nx) + " x " + to_string(this->Ny);
	}

	string FileNamePart()
	{
		return "n" + to_string(this->Nx);
	}

	double H()
	{
		return sqrt(2) / (double)this->Nx;
	}

	void CoarsenMesh(CoarseningStrategy strategy)
	{
		//if (strategy == CoarseningStrategy::Standard)
			//CoarsenByAgglomerationAndMergeColinearFaces();
		//else if (strategy == CoarseningStrategy::Agglomeration)
			//CoarsenByAgglomerationAndKeepFineFaces();
		//else
			assert(false && "Coarsening strategy not implemented!");
		this->CoarseMesh->SetDiffusionCoefficient(this->_diffusionPartition);
		this->CoarseMesh->SetBoundaryConditions(this->_boundaryConditions);
	}

	/*void CoarsenByAgglomerationAndMergeColinearFaces()
	{
		BigNumber nx = this->Nx;
		BigNumber ny = this->Ny;

		if (nx % 2 != 0 || ny % 2 != 0)
		{
			cout << "Error: impossible to build coarse mesh. Nx and Ny must be even: Nx = " << nx << ", Ny = " << ny << "." << endl;
		}
		else
		{
			TriangularMesh* coarseMesh = new TriangularMesh(nx / 2, ny / 2);

			for (BigNumber i = 0; i < ny; ++i)
			{
				for (BigNumber j = 0; j < nx; ++j)
				{
					Triangle* fineElement = dynamic_cast<Triangle*>(this->Elements[i*nx + j]);
					Triangle* coarseElement = dynamic_cast<Triangle*>(coarseMesh->Elements[(i / 2) * coarseMesh->Nx + j / 2]);

					coarseElement->FinerElements.push_back(fineElement);
					fineElement->CoarserElement = coarseElement;
					if (i % 2 == 0 && !fineElement->NorthFace->IsDomainBoundary)
					{
						fineElement->NorthFace->IsRemovedOnCoarserGrid = true;
						coarseElement->FinerFacesRemoved.push_back(fineElement->NorthFace);
						coarseElement->SouthFace->FinerFaces.push_back(fineElement->SouthFace);
					}
					if (i == ny - 1)
						coarseElement->NorthFace->FinerFaces.push_back(fineElement->NorthFace);

					if (j % 2 == 0 && !fineElement->EastFace->IsDomainBoundary)
					{
						fineElement->EastFace->IsRemovedOnCoarserGrid = true;
						coarseElement->FinerFacesRemoved.push_back(fineElement->EastFace);
						coarseElement->WestFace->FinerFaces.push_back(fineElement->WestFace);
					}
					if (j == nx - 1)
						coarseElement->EastFace->FinerFaces.push_back(fineElement->EastFace);
				}
			}

			this->CoarseMesh = coarseMesh;
		}
	}*/

	/*void CoarsenByAgglomerationAndKeepFineFaces()
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
					Triangle* bottomLeftElement = dynamic_cast<Triangle*>(this->Elements[2 * i * nx + 2 * j]);
					Triangle* bottomRightElement = dynamic_cast<Triangle*>(this->Elements[2 * i * nx + 2 * j + 1]);
					Triangle* topLeftElement = dynamic_cast<Triangle*>(this->Elements[(2 * i + 1) * nx + 2 * j]);
					Triangle* topRightElement = dynamic_cast<Triangle*>(this->Elements[(2 * i + 1) * nx + 2 * j + 1]);

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
	}*/
};