#pragma once
#include <vector>
#include "Parallelepiped.h"
#include "RectangularFace.h"
#include "../Mesh.h"

using namespace std;

class CartesianGrid3D : public Mesh<3>
{
public:
	BigNumber Nx;
	BigNumber Ny;
	BigNumber Nz;

	CartesianGrid3D(BigNumber nx, BigNumber ny, BigNumber nz) : Mesh()
	{
		// nx = ny = nz falls down to cubic elements
		this->Nx = nx;
		this->Ny = ny;
		this->Nz = nz;

		double hx = 1 / (double)nx;
		double hy = 1 / (double)ny;
		double hz = 1 / (double)nz;

		//----------//
		// Vertices //
		//----------//

		this->Vertices.reserve((nx + 1) * (ny + 1) * (nz + 1));

		for (BigNumber iz = 0; iz < nz + 1; ++iz)
		{
			for (BigNumber iy = 0; iy < ny + 1; ++iy)
			{
				for (BigNumber ix = 0; ix < nx + 1; ++ix)
				{
					Vertex* vertex = new Vertex(indexV(ix, iy, iz), ix * hx, iy * hy, iz * hz);
					this->Vertices.push_back(vertex);
				}
			}
		}

		//----------//
		// Elements //
		//----------//

		this->Elements.reserve(nx * ny * nz);

		for (BigNumber iz = 0; iz < nz; ++iz)
		{
			for (BigNumber iy = 0; iy < ny; ++iy)
			{
				for (BigNumber ix = 0; ix < nx; ++ix)
				{
					Vertex* backLeftBottomCorner   = Vertices[indexV(ix,   iy,   iz)];
					Vertex* frontLeftBottomCorner  = Vertices[indexV(ix+1, iy,   iz)];
					Vertex* backRightBottomCorner  = Vertices[indexV(ix,   iy+1, iz)];
					Vertex* backLeftTopCorner      = Vertices[indexV(ix,   iy,   iz+1)];
					Vertex* frontLeftTopCorner     = Vertices[indexV(ix+1, iy,   iz+1)];
					Vertex* backRightTopCorner     = Vertices[indexV(ix,   iy+1, iz+1)]; 
					Vertex* frontRightBottomCorner = Vertices[indexV(ix+1, iy+1, iz)];
					Vertex* frontRightTopCorner    = Vertices[indexV(ix+1, iy+1, iz+1)];
					Parallelepiped* element = new Parallelepiped(index(ix, iy, iz), backLeftBottomCorner, frontLeftBottomCorner, backRightBottomCorner, backLeftTopCorner, frontLeftTopCorner, backRightTopCorner, frontRightBottomCorner, frontRightTopCorner);
					this->Elements.push_back(element);
				}
			}
		}

		//-------//
		// Faces //
		//-------//

		this->Faces.reserve(nx * (ny + 1) * (nz + 1) + ny * (nx + 1) * (nz + 1) + nz * (nx + 1) * (ny + 1));
		BigNumber numberInterface = 0;

		for (BigNumber iy = 0; iy < ny; ++iy)
		{
			for (BigNumber ix = 0; ix < nx; ++ix)
			{
				// Bottom boundary
				Parallelepiped* element = dynamic_cast<Parallelepiped*>(this->Elements[index(ix, iy, 0)]);
				RectangularFace* bottomBoundary = new RectangularFace(numberInterface++, element->BackLeftBottomCorner, element->FrontLeftBottomCorner, element->BackRightBottomCorner, element, CartesianShapeOrientation::InXOY);
				this->Faces.push_back(bottomBoundary);
				this->BoundaryFaces.push_back(bottomBoundary);
				element->SetBottomFace(bottomBoundary);

				// Top boundary
				element = dynamic_cast<Parallelepiped*>(this->Elements[index(ix, iy, nz - 1)]);
				RectangularFace* topBoundary = new RectangularFace(numberInterface++, element->BackLeftTopCorner, element->FrontLeftTopCorner, element->BackRightTopCorner, element, CartesianShapeOrientation::InXOY);
				this->Faces.push_back(topBoundary);
				this->BoundaryFaces.push_back(topBoundary);
				element->SetTopFace(topBoundary);
			}
		}

		for (BigNumber iz = 0; iz < nz; ++iz)
		{
			for (BigNumber ix = 0; ix < nx; ++ix)
			{
				// Left boundary
				Parallelepiped* element = dynamic_cast<Parallelepiped*>(this->Elements[index(ix, 0, iz)]);
				RectangularFace* leftBoundary = new RectangularFace(numberInterface++, element->BackLeftBottomCorner, element->FrontLeftBottomCorner, element->BackLeftTopCorner, element, CartesianShapeOrientation::InXOZ);
				this->Faces.push_back(leftBoundary);
				this->BoundaryFaces.push_back(leftBoundary);
				element->SetLeftFace(leftBoundary);

				// Right boundary
				element = dynamic_cast<Parallelepiped*>(this->Elements[index(ix, ny - 1, iz)]);
				RectangularFace* rightBoundary = new RectangularFace(numberInterface++, element->BackRightBottomCorner, element->FrontRightBottomCorner, element->BackRightTopCorner, element, CartesianShapeOrientation::InXOZ);
				this->Faces.push_back(rightBoundary);
				this->BoundaryFaces.push_back(rightBoundary);
				element->SetRightFace(rightBoundary);
			}
		}

		for (BigNumber iz = 0; iz < nz; ++iz)
		{
			for (BigNumber iy = 0; iy < ny; ++iy)
			{
				// Back boundary
				Parallelepiped* element = dynamic_cast<Parallelepiped*>(this->Elements[index(0, iy, iz)]);
				RectangularFace* backBoundary = new RectangularFace(numberInterface++, element->BackLeftBottomCorner, element->BackRightBottomCorner, element->BackLeftTopCorner, element, CartesianShapeOrientation::InYOZ);
				this->Faces.push_back(backBoundary);
				this->BoundaryFaces.push_back(backBoundary);
				element->SetBackFace(backBoundary);

				// Front boundary
				element = dynamic_cast<Parallelepiped*>(this->Elements[index(nx - 1, iy, iz)]);
				RectangularFace* frontBoundary = new RectangularFace(numberInterface++, element->FrontLeftBottomCorner, element->FrontRightBottomCorner, element->FrontLeftTopCorner, element, CartesianShapeOrientation::InYOZ);
				this->Faces.push_back(frontBoundary);
				this->BoundaryFaces.push_back(frontBoundary);
				element->SetFrontFace(frontBoundary);
			}
		}

		for (BigNumber iz = 0; iz < nz; ++iz)
		{
			for (BigNumber iy = 0; iy < ny; ++iy)
			{
				for (BigNumber ix = 0; ix < nx; ++ix)
				{
					Parallelepiped* element = dynamic_cast<Parallelepiped*>(this->Elements[index(ix, iy, iz)]);
					if (ix != nx - 1)
					{
						// Front
						Parallelepiped* frontNeighbour = dynamic_cast<Parallelepiped*>(this->Elements[index(ix+1, iy, iz)]);
						RectangularFace* interface = new RectangularFace(numberInterface++, element->FrontLeftBottomCorner, element->FrontRightBottomCorner, element->FrontLeftTopCorner, element, frontNeighbour, CartesianShapeOrientation::InYOZ);
						this->Faces.push_back(interface);
						this->InteriorFaces.push_back(interface);
						element->SetFrontFace(interface);
						frontNeighbour->SetBackFace(interface);
					}
					if (iy != ny - 1)
					{
						// Right
						Parallelepiped* rightNeighbour = dynamic_cast<Parallelepiped*>(this->Elements[index(ix, iy+1, iz)]);
						RectangularFace* interface = new RectangularFace(numberInterface++, element->BackRightBottomCorner, element->FrontRightBottomCorner, element->BackRightTopCorner, element, rightNeighbour, CartesianShapeOrientation::InXOZ);
						this->Faces.push_back(interface);
						this->InteriorFaces.push_back(interface);
						element->SetRightFace(interface);
						rightNeighbour->SetLeftFace(interface);
					}
					if (iz != nz - 1)
					{
						// Top
						Parallelepiped* topNeighbour = dynamic_cast<Parallelepiped*>(this->Elements[index(ix, iy, iz+1)]);
						RectangularFace* interface = new RectangularFace(numberInterface++, element->BackLeftTopCorner, element->FrontLeftTopCorner, element->BackRightTopCorner, element, topNeighbour, CartesianShapeOrientation::InXOY);
						this->Faces.push_back(interface);
						this->InteriorFaces.push_back(interface);
						element->SetTopFace(interface);
						topNeighbour->SetBottomFace(interface);
					}
				}
			}
		}

	}

	string Description()
	{
		return "Cartesian " + to_string(this->Nx) + " x " + to_string(this->Ny) + " x " + to_string(this->Nz);
	}

	string FileNamePart()
	{
		return "cart-n" + to_string(this->Nx);
	}

	double H()
	{
		return 1.0 / this->Nx;
	}

	double Regularity() override
	{
		return 1;
	}

	void CoarsenMesh(CoarseningStrategy strategy)
	{
		if (strategy == CoarseningStrategy::StandardCoarsening)
			StandardCoarsening();
		else
			Utils::FatalError("Coarsening strategy not implemented!");
		if (this->_diffusionPartition)
		{
			this->CoarseMesh->SetDiffusionCoefficient(this->_diffusionPartition);
			this->CoarseMesh->SetBoundaryConditions(this->_boundaryConditions);
		}
	}

	void RefineMesh(CoarseningStrategy strategy) override
	{
		Utils::FatalError("Refinement strategy not implemented!");
	}

private:
	void StandardCoarsening()
	{
		BigNumber nx = this->Nx;
		BigNumber ny = this->Ny;
		BigNumber nz = this->Nz;

		if (nx % 2 != 0 || ny % 2 != 0 || nz % 2 != 0)
		{
			cout << "Error: impossible to build coarse mesh. Nx, Ny and Nz must be even: Nx = " << nx << ", Ny = " << ny << ", Nz = " << nz << "." << endl;
		}
		else
		{
			CartesianGrid3D* coarseMesh = new CartesianGrid3D(nx / 2, ny / 2, nz / 2);
			coarseMesh->ComesFrom.CS = CoarseningStrategy::StandardCoarsening;
			coarseMesh->ComesFrom.nFineElementsByCoarseElement = 8;
			coarseMesh->ComesFrom.nFineFacesAddedByCoarseElement = 12;
			coarseMesh->ComesFrom.nFineFacesByKeptCoarseFace = 4;

			for (BigNumber iz = 0; iz < nz; ++iz)
			{
				for (BigNumber iy = 0; iy < ny; ++iy)
				{
					for (BigNumber ix = 0; ix < nx; ++ix)
					{
						Parallelepiped* fineElement = dynamic_cast<Parallelepiped*>(this->Elements[index(ix, iy, iz)]);
						Parallelepiped* coarseElement = dynamic_cast<Parallelepiped*>(coarseMesh->Elements[coarseMesh->index(ix / 2, iy / 2, iz / 2)]);

						coarseElement->FinerElements.push_back(fineElement);
						fineElement->CoarserElement = coarseElement;
						if (iz % 2 == 0 && !fineElement->TopFace->IsDomainBoundary)
						{
							fineElement->TopFace->IsRemovedOnCoarserGrid = true;
							coarseElement->FinerFacesRemoved.push_back(fineElement->TopFace);
							coarseElement->BottomFace->FinerFaces.push_back(fineElement->BottomFace);
							fineElement->BottomFace->CoarseFace = coarseElement->BottomFace;
						}
						if (iz == nz - 1)
						{
							coarseElement->TopFace->FinerFaces.push_back(fineElement->TopFace);
							fineElement->TopFace->CoarseFace = coarseElement->TopFace;
						}

						if (iy % 2 == 0 && !fineElement->RightFace->IsDomainBoundary)
						{
							fineElement->RightFace->IsRemovedOnCoarserGrid = true;
							coarseElement->FinerFacesRemoved.push_back(fineElement->RightFace);
							coarseElement->LeftFace->FinerFaces.push_back(fineElement->LeftFace);
							fineElement->LeftFace->CoarseFace = coarseElement->LeftFace;
						}
						if (iy == ny - 1)
						{
							coarseElement->RightFace->FinerFaces.push_back(fineElement->RightFace);
							fineElement->RightFace->CoarseFace = coarseElement->RightFace;
						}

						if (ix % 2 == 0 && !fineElement->FrontFace->IsDomainBoundary)
						{
							fineElement->FrontFace->IsRemovedOnCoarserGrid = true;
							coarseElement->FinerFacesRemoved.push_back(fineElement->FrontFace);
							coarseElement->BackFace->FinerFaces.push_back(fineElement->BackFace);
							fineElement->BackFace->CoarseFace = coarseElement->BackFace;
						}
						if (ix == nx - 1)
						{
							coarseElement->FrontFace->FinerFaces.push_back(fineElement->FrontFace);
							fineElement->FrontFace->CoarseFace = coarseElement->FrontFace;
						}
					}
				}
			}

			this->CoarseMesh = coarseMesh;
		}
	}

private:
	inline BigNumber indexV(BigNumber ix, BigNumber iy, BigNumber iz)
	{
		return iz * (Ny + 1)*(Nx + 1) + iy * (Nx + 1) + ix;
	}
	inline BigNumber index(BigNumber ix, BigNumber iy, BigNumber iz)
	{
		return iz * Ny*Nx + iy * Nx + ix;
	}
};