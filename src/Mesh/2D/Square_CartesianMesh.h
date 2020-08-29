#pragma once
#include <vector>
#include "RectangularElement.h"
#include "RectangularPolygonalElement.h"
#include "Square_CartesianPolygonalMesh.h"
#include "CartesianEdge.h"
#include "../PolyhedralMesh.h"
using namespace std;

class Square_CartesianMesh : public PolyhedralMesh<2>
{
private:
public:
	BigNumber Nx;
	BigNumber Ny;
	bool With4Quadrants;

	Square_CartesianMesh(BigNumber nx, BigNumber ny, bool with4Quadrants = false, bool buildMesh = true) : PolyhedralMesh()
	{
		// nx = ny falls down to square elements
		this->Nx = nx;
		this->Ny = ny;
		this->With4Quadrants = with4Quadrants;

		if (with4Quadrants && (nx % 2 == 1 || ny % 2 == 1))
			Utils::FatalError("Building the mesh for a square with 4 quadrants requires the number of subdivisions in each direction to be even.");

		if (buildMesh)
			Build();
	}

	void Build()
	{
		BigNumber nx = this->Nx;
		BigNumber ny = this->Ny;
		double hx = 1 / (double)nx;
		double hy = 1 / (double)ny;

		// Boundary parts
		if (this->BoundaryParts.empty())
		{
			this->BoundaryParts.push_back(new BoundaryGroup(1, "bottomBoundary"));
			this->BoundaryParts.push_back(new BoundaryGroup(2, "rightBoundary"));
			this->BoundaryParts.push_back(new BoundaryGroup(3, "topBoundary"));
			this->BoundaryParts.push_back(new BoundaryGroup(4, "leftBoundary"));
		}
		BoundaryGroup* squareBottomBoundary = this->BoundaryParts[0];
		BoundaryGroup* squareRightBoundary = this->BoundaryParts[1];
		BoundaryGroup* squareTopBoundary = this->BoundaryParts[2];
		BoundaryGroup* squareLeftBoundary = this->BoundaryParts[3];

		//----------//
		// Vertices //
		//----------//

		this->Vertices.reserve((nx + 1) * (ny + 1));
		for (BigNumber iy = 0; iy < ny + 1; ++iy)
		{
			for (BigNumber ix = 0; ix < nx + 1; ++ix)
			{
				BigNumber number = indexV(ix, iy);
				MeshVertex<2>* vertex = new MeshVertex<2>(number, ix * hx, iy * hy);
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
				MeshVertex<2>* bottomLeftCorner  = static_cast<MeshVertex<2>*>(Vertices[indexV(ix,     iy    )]);
				MeshVertex<2>* topLeftCorner     = static_cast<MeshVertex<2>*>(Vertices[indexV(ix,     iy + 1)]);
				MeshVertex<2>* topRightCorner    = static_cast<MeshVertex<2>*>(Vertices[indexV(ix + 1, iy + 1)]);
				MeshVertex<2>* bottomRightCorner = static_cast<MeshVertex<2>*>(Vertices[indexV(ix + 1, iy    )]);
				RectangularElement* rectangle = new RectangularElement(number, bottomLeftCorner, topLeftCorner, topRightCorner, bottomRightCorner);
				this->Elements.push_back(rectangle);

				bottomLeftCorner->Elements.push_back(rectangle);
				topLeftCorner->Elements.push_back(rectangle);
				topRightCorner->Elements.push_back(rectangle);
				bottomRightCorner->Elements.push_back(rectangle);
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
			RectangularElement* rectangle = dynamic_cast<RectangularElement*>(this->Elements[index(ix, 0)]);
			CartesianEdge* southBoundary = new CartesianEdge(numberInterface++, rectangle->BottomLeftCorner, rectangle->BottomRightCorner, rectangle, CartesianShapeOrientation::Horizontal);
			this->Faces.push_back(southBoundary);
			this->BoundaryFaces.push_back(southBoundary);
			rectangle->SetSouthInterface(southBoundary);
			static_cast<MeshVertex<2>*>(rectangle->BottomLeftCorner)->Faces.push_back(southBoundary);
			static_cast<MeshVertex<2>*>(rectangle->BottomRightCorner)->Faces.push_back(southBoundary);
			southBoundary->BoundaryPart = squareBottomBoundary;

			// North boundary
			rectangle = dynamic_cast<RectangularElement*>(this->Elements[index(ix, ny - 1)]);
			CartesianEdge* northBoundary = new CartesianEdge(numberInterface++, rectangle->TopLeftCorner, rectangle->TopRightCorner, rectangle, CartesianShapeOrientation::Horizontal);
			this->Faces.push_back(northBoundary);
			this->BoundaryFaces.push_back(northBoundary);
			rectangle->SetNorthInterface(northBoundary);
			static_cast<MeshVertex<2>*>(rectangle->TopLeftCorner)->Faces.push_back(northBoundary);
			static_cast<MeshVertex<2>*>(rectangle->TopRightCorner)->Faces.push_back(northBoundary);
			northBoundary->BoundaryPart = squareTopBoundary;
		}

		for (BigNumber iy = 0; iy < ny; ++iy)
		{
			// West boundary
			RectangularElement* rectangle = dynamic_cast<RectangularElement*>(this->Elements[index(0, iy)]);
			CartesianEdge* westBoundary = new CartesianEdge(numberInterface++, rectangle->BottomLeftCorner, rectangle->TopLeftCorner, rectangle, CartesianShapeOrientation::Vertical);
			this->Faces.push_back(westBoundary);
			this->BoundaryFaces.push_back(westBoundary);
			rectangle->SetWestInterface(westBoundary);
			static_cast<MeshVertex<2>*>(rectangle->BottomLeftCorner)->Faces.push_back(westBoundary);
			static_cast<MeshVertex<2>*>(rectangle->TopLeftCorner)->Faces.push_back(westBoundary);
			westBoundary->BoundaryPart = squareLeftBoundary;

			// East boundary
			rectangle = dynamic_cast<RectangularElement*>(this->Elements[index(nx-1, iy)]);
			CartesianEdge* eastBoundary = new CartesianEdge(numberInterface++, rectangle->BottomRightCorner, rectangle->TopRightCorner, rectangle, CartesianShapeOrientation::Vertical);
			this->Faces.push_back(eastBoundary);
			this->BoundaryFaces.push_back(eastBoundary);
			rectangle->SetEastInterface(eastBoundary);
			static_cast<MeshVertex<2>*>(rectangle->BottomRightCorner)->Faces.push_back(eastBoundary);
			static_cast<MeshVertex<2>*>(rectangle->TopRightCorner)->Faces.push_back(eastBoundary);
			eastBoundary->BoundaryPart = squareRightBoundary;
		}

		for (BigNumber iy = 0; iy < ny; iy++)
		{
			for (BigNumber ix = 0; ix < nx; ix++)
			{
				RectangularElement* element = dynamic_cast<RectangularElement*>(this->Elements[index(ix, iy)]);
				if (ix != nx - 1)
				{
					// East
					RectangularElement* eastNeighbour = dynamic_cast<RectangularElement*>(this->Elements[index(ix + 1, iy)]);
					CartesianEdge* interface = new CartesianEdge(numberInterface++, eastNeighbour->BottomLeftCorner, eastNeighbour->TopLeftCorner, element, eastNeighbour, CartesianShapeOrientation::Vertical);
					this->Faces.push_back(interface);
					this->InteriorFaces.push_back(interface);
					element->SetEastInterface(interface);
					eastNeighbour->SetWestInterface(interface);
					static_cast<MeshVertex<2>*>(eastNeighbour->BottomLeftCorner)->Faces.push_back(interface);
					static_cast<MeshVertex<2>*>(eastNeighbour->TopLeftCorner)->Faces.push_back(interface);
				}
				if (iy != ny - 1)
				{
					// North
					RectangularElement* northNeighbour = dynamic_cast<RectangularElement*>(this->Elements[index(ix, iy + 1)]);
					CartesianEdge* interface = new CartesianEdge(numberInterface++, northNeighbour->BottomLeftCorner, northNeighbour->BottomRightCorner, element, northNeighbour, CartesianShapeOrientation::Horizontal);
					this->Faces.push_back(interface);
					this->InteriorFaces.push_back(interface);
					element->SetNorthInterface(interface);
					northNeighbour->SetSouthInterface(interface);
					static_cast<MeshVertex<2>*>(northNeighbour->BottomLeftCorner)->Faces.push_back(interface);
					static_cast<MeshVertex<2>*>(northNeighbour->BottomRightCorner)->Faces.push_back(interface);
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
		string geo = this->With4Quadrants ? "square4quadrants" : "square";
		return geo + "-inhouse-cart-n" + to_string(this->Nx);
	}

	string GeometryDescription() override
	{
		return this->With4Quadrants ? "Square with 4 quadrants" : "Square";
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
			PolyhedralMesh<2>::CoarsenMesh(strategy);

		this->CoarseMesh->SetDiffusionField(this->_diffusionField);
	}

	void StandardCoarsening()
	{
		BigNumber nx = this->Nx;
		BigNumber ny = this->Ny;

		if (nx % 2 != 0 || ny % 2 != 0)
		{
			cout << "Error: impossible to build coarse mesh. Nx and Ny must be even: Nx = " << nx << ", Ny = " << ny << "." << endl;
			return;
		}

		Square_CartesianMesh* coarseMesh = new Square_CartesianMesh(nx / 2, ny / 2, this->With4Quadrants, false);
		this->InitializeCoarsening(coarseMesh);
		coarseMesh->ComesFrom.CS = CoarseningStrategy::StandardCoarsening;
		coarseMesh->ComesFrom.nFineElementsByCoarseElement = 4;
		coarseMesh->ComesFrom.nFineFacesAddedByCoarseElement = 4;
		coarseMesh->ComesFrom.nFineFacesByKeptCoarseFace = 2;
		coarseMesh->Build();

		for (BigNumber i = 0; i < ny; ++i)
		{
			for (BigNumber j = 0; j < nx; ++j)
			{
				RectangularElement* fineElement = dynamic_cast<RectangularElement*>(this->Elements[i*nx + j]);
				RectangularElement* coarseElement = dynamic_cast<RectangularElement*>(coarseMesh->Elements[(i / 2) * coarseMesh->Nx + j / 2]);

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
		this->FinalizeCoarsening();
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

		Square_CartesianPolygonalMesh* coarseMesh = new Square_CartesianPolygonalMesh(this->With4Quadrants);
		this->InitializeCoarsening(coarseMesh);
		coarseMesh->Nx = nx / 2;
		coarseMesh->Ny = ny / 2;
		coarseMesh->ComesFrom.CS = CoarseningStrategy::AgglomerationCoarsening;
		coarseMesh->ComesFrom.nFineElementsByCoarseElement = 4;
		coarseMesh->ComesFrom.nFineFacesAddedByCoarseElement = 4;
		coarseMesh->ComesFrom.nFineFacesByKeptCoarseFace = 1;

		// Elements //
		BigNumber faceNumber = 0;
		for (BigNumber i = 0; i < ny / 2; ++i)
		{
			for (BigNumber j = 0; j < nx / 2; ++j)
			{
				RectangularElement* bottomLeftElement = dynamic_cast<RectangularElement*>(this->Elements[2 * i * nx + 2 * j]);
				RectangularElement* bottomRightElement = dynamic_cast<RectangularElement*>(this->Elements[2 * i * nx + 2 * j + 1]);
				RectangularElement* topLeftElement = dynamic_cast<RectangularElement*>(this->Elements[(2 * i + 1) * nx + 2 * j]);
				RectangularElement* topRightElement = dynamic_cast<RectangularElement*>(this->Elements[(2 * i + 1) * nx + 2 * j + 1]);

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
					Face<2>* bottomWestFace = bottomLeftElement->WestFace->CreateSameGeometricFace(faceNumber++, coarseElement);
					bottomLeftElement->WestFace->CoarseFace = bottomWestFace;
					bottomWestFace->FinerFaces.push_back(bottomLeftElement->WestFace);
					Face<2>* topWestFace = topLeftElement->WestFace->CreateSameGeometricFace(faceNumber++, coarseElement);
					topLeftElement->WestFace->CoarseFace = topWestFace;
					topWestFace->FinerFaces.push_back(topLeftElement->WestFace);
					coarseElement->AddWestFace(bottomWestFace);
					coarseElement->AddWestFace(topWestFace);
					coarseMesh->AddFace(bottomWestFace);
					coarseMesh->AddFace(topWestFace);
				}
				else
				{
					RectangularPolygonalElement* leftNeighbour = dynamic_cast<RectangularPolygonalElement*>(coarseMesh->Elements[i * nx / 2 + j - 1]);
					coarseElement->SetWestFacesFromNeighbour(leftNeighbour);
				}

				// North faces
				Face<2>* leftNorthFace = topLeftElement->NorthFace->CreateSameGeometricFace(faceNumber++, coarseElement);
				topLeftElement->NorthFace->CoarseFace = leftNorthFace;
				leftNorthFace->FinerFaces.push_back(topLeftElement->NorthFace);
				Face<2>* rightNorthFace = topRightElement->NorthFace->CreateSameGeometricFace(faceNumber++, coarseElement);
				topRightElement->NorthFace->CoarseFace = rightNorthFace;
				rightNorthFace->FinerFaces.push_back(topRightElement->NorthFace);
				coarseElement->AddNorthFace(leftNorthFace);
				coarseElement->AddNorthFace(rightNorthFace);
				coarseMesh->AddFace(leftNorthFace);
				coarseMesh->AddFace(rightNorthFace);

				// East faces
				Face<2>* topEastFace = topRightElement->EastFace->CreateSameGeometricFace(faceNumber++, coarseElement);
				topRightElement->EastFace->CoarseFace = topEastFace;
				topEastFace->FinerFaces.push_back(topRightElement->EastFace);
				Face<2>* bottomEastFace = bottomRightElement->EastFace->CreateSameGeometricFace(faceNumber++, coarseElement);
				bottomRightElement->EastFace->CoarseFace = bottomEastFace;
				bottomEastFace->FinerFaces.push_back(bottomRightElement->EastFace);
				coarseElement->AddEastFace(topEastFace);
				coarseElement->AddEastFace(bottomEastFace);
				coarseMesh->AddFace(topEastFace);
				coarseMesh->AddFace(bottomEastFace);

				// South faces
				if (i == 0)
				{
					Face<2>* rightSouthFace = bottomRightElement->SouthFace->CreateSameGeometricFace(faceNumber++, coarseElement);
					bottomRightElement->SouthFace->CoarseFace = rightSouthFace;
					rightSouthFace->FinerFaces.push_back(bottomRightElement->SouthFace);
					Face<2>* leftSouthFace = bottomLeftElement->SouthFace->CreateSameGeometricFace(faceNumber++, coarseElement);
					bottomLeftElement->SouthFace->CoarseFace = leftSouthFace;
					leftSouthFace->FinerFaces.push_back(bottomLeftElement->SouthFace);
					coarseElement->AddSouthFace(rightSouthFace);
					coarseElement->AddSouthFace(leftSouthFace);
					coarseMesh->AddFace(rightSouthFace);
					coarseMesh->AddFace(leftSouthFace);
				}
				else
				{
					RectangularPolygonalElement* bottomNeighbour = dynamic_cast<RectangularPolygonalElement*>(coarseMesh->Elements[(i - 1) * nx / 2 + j]);
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
		this->FinalizeCoarsening();
	}

};