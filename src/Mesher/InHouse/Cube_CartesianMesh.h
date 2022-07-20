#pragma once
#include "../../Mesh/3D/RectangularFace.h"
#include "../../Mesh/PolyhedralMesh.h"
#include "../CubeGeometry.h"

using namespace std;

class Cube_CartesianMesh : public PolyhedralMesh<3>
{
public:
	BigNumber Nx;
	BigNumber Ny;
	BigNumber Nz;

	Cube_CartesianMesh(BigNumber nx, BigNumber ny, BigNumber nz, bool buildMesh = true) : PolyhedralMesh()
	{
		// nx = ny = nz falls down to cubic elements
		this->Nx = nx;
		this->Ny = ny;
		this->Nz = nz;

		if (buildMesh)
			Build();
	}

private:
	void Build()
	{
		BigNumber nx = this->Nx;
		BigNumber ny = this->Ny;
		BigNumber nz = this->Nz;

		double hx = 1 / (double)nx;
		double hy = 1 / (double)ny;
		double hz = 1 / (double)nz;

		// Physical parts
		PhysicalGroup<3>* domain = nullptr;
		if (this->PhysicalParts.empty())
			this->PhysicalParts = CubeGeometry::PhysicalParts();
		domain = this->PhysicalParts[0];

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
					MeshVertex<3>* vertex = new MeshVertex<3>(indexV(ix, iy, iz), ix * hx, iy * hy, iz * hz);
					this->Vertices.push_back(vertex);
				}
			}
		}

		//----------//
		// Elements //
		//----------//

		this->_parallelepipedElements.reserve(nx * ny * nz);
		this->Elements.reserve(nx * ny * nz);

		for (BigNumber iz = 0; iz < nz; ++iz)
		{
			for (BigNumber iy = 0; iy < ny; ++iy)
			{
				for (BigNumber ix = 0; ix < nx; ++ix)
				{
					MeshVertex<3>* backLeftBottomCorner   = static_cast<MeshVertex<3>*>(Vertices[indexV(ix,   iy,   iz)]);
					MeshVertex<3>* frontLeftBottomCorner  = static_cast<MeshVertex<3>*>(Vertices[indexV(ix+1, iy,   iz)]);
					MeshVertex<3>* backRightBottomCorner  = static_cast<MeshVertex<3>*>(Vertices[indexV(ix,   iy+1, iz)]);
					MeshVertex<3>* backLeftTopCorner      = static_cast<MeshVertex<3>*>(Vertices[indexV(ix,   iy,   iz+1)]);
					MeshVertex<3>* frontLeftTopCorner     = static_cast<MeshVertex<3>*>(Vertices[indexV(ix+1, iy,   iz+1)]);
					MeshVertex<3>* backRightTopCorner     = static_cast<MeshVertex<3>*>(Vertices[indexV(ix,   iy+1, iz+1)]);
					MeshVertex<3>* frontRightBottomCorner = static_cast<MeshVertex<3>*>(Vertices[indexV(ix+1, iy+1, iz)]);
					MeshVertex<3>* frontRightTopCorner    = static_cast<MeshVertex<3>*>(Vertices[indexV(ix+1, iy+1, iz+1)]);
					auto i = index(ix, iy, iz);
					this->_parallelepipedElements.emplace_back(i, backLeftBottomCorner, frontLeftBottomCorner, backRightBottomCorner, backLeftTopCorner, frontLeftTopCorner, backRightTopCorner, frontRightBottomCorner, frontRightTopCorner);
					ParallelepipedElement* element = &this->_parallelepipedElements.back();
					//ParallelepipedElement* element = new ParallelepipedElement(index(ix, iy, iz), backLeftBottomCorner, frontLeftBottomCorner, backRightBottomCorner, backLeftTopCorner, frontLeftTopCorner, backRightTopCorner, frontRightBottomCorner, frontRightTopCorner);
					element->PhysicalPart = domain;
					this->Elements.push_back(element);

					backLeftBottomCorner->Elements.push_back(element);
					frontLeftBottomCorner->Elements.push_back(element);
					backRightBottomCorner->Elements.push_back(element);
					backLeftTopCorner->Elements.push_back(element);
					frontLeftTopCorner->Elements.push_back(element);
					backRightTopCorner->Elements.push_back(element);
					frontRightBottomCorner->Elements.push_back(element);
					frontRightTopCorner->Elements.push_back(element);
				}
			}
		}

		//-------//
		// Faces //
		//-------//

		BigNumber nFaces = nx * (ny + 1) * (nz + 1) + ny * (nx + 1) * (nz + 1) + nz * (nx + 1) * (ny + 1);
		this->_rectangularFaces.reserve(nFaces);
		this->Faces.reserve(nFaces);
		BigNumber numberInterface = 0;

		for (BigNumber iy = 0; iy < ny; ++iy)
		{
			for (BigNumber ix = 0; ix < nx; ++ix)
			{
				// Bottom boundary
				ParallelepipedElement* element = dynamic_cast<ParallelepipedElement*>(this->Elements[index(ix, iy, 0)]);
				this->_rectangularFaces.emplace_back(numberInterface++, element->BackLeftBottomCorner, element->FrontLeftBottomCorner, element->BackRightBottomCorner, element->FrontRightBottomCorner, element, CartesianShapeOrientation::InXOY);
				RectangularFace* bottomBoundary = &this->_rectangularFaces.back();
				this->Faces.push_back(bottomBoundary);
				this->BoundaryFaces.push_back(bottomBoundary);
				element->SetBottomFace(bottomBoundary);

				static_cast<MeshVertex<3>*>(element->BackLeftBottomCorner)->Faces.push_back(bottomBoundary);
				static_cast<MeshVertex<3>*>(element->FrontLeftBottomCorner)->Faces.push_back(bottomBoundary);
				static_cast<MeshVertex<3>*>(element->FrontRightBottomCorner)->Faces.push_back(bottomBoundary);
				static_cast<MeshVertex<3>*>(element->BackLeftBottomCorner)->Faces.push_back(bottomBoundary);

				// Top boundary
				element = dynamic_cast<ParallelepipedElement*>(this->Elements[index(ix, iy, nz - 1)]);
				this->_rectangularFaces.emplace_back(numberInterface++, element->BackLeftTopCorner, element->FrontLeftTopCorner, element->BackRightTopCorner, element->FrontRightTopCorner, element, CartesianShapeOrientation::InXOY);
				RectangularFace* topBoundary = &this->_rectangularFaces.back();
				this->Faces.push_back(topBoundary);
				this->BoundaryFaces.push_back(topBoundary);
				element->SetTopFace(topBoundary);

				static_cast<MeshVertex<3>*>(element->BackLeftTopCorner)->Faces.push_back(topBoundary);
				static_cast<MeshVertex<3>*>(element->FrontLeftTopCorner)->Faces.push_back(topBoundary);
				static_cast<MeshVertex<3>*>(element->FrontRightTopCorner)->Faces.push_back(topBoundary);
				static_cast<MeshVertex<3>*>(element->BackLeftTopCorner)->Faces.push_back(topBoundary);
			}
		}

		for (BigNumber iz = 0; iz < nz; ++iz)
		{
			for (BigNumber ix = 0; ix < nx; ++ix)
			{
				// Left boundary
				ParallelepipedElement* element = dynamic_cast<ParallelepipedElement*>(this->Elements[index(ix, 0, iz)]);
				this->_rectangularFaces.emplace_back(numberInterface++, element->BackLeftBottomCorner, element->FrontLeftBottomCorner, element->BackLeftTopCorner, element->FrontLeftTopCorner, element, CartesianShapeOrientation::InXOZ);
				RectangularFace* leftBoundary = &this->_rectangularFaces.back();
				this->Faces.push_back(leftBoundary);
				this->BoundaryFaces.push_back(leftBoundary);
				element->SetLeftFace(leftBoundary);

				static_cast<MeshVertex<3>*>(element->BackLeftTopCorner)->Faces.push_back(leftBoundary);
				static_cast<MeshVertex<3>*>(element->FrontLeftTopCorner)->Faces.push_back(leftBoundary);
				static_cast<MeshVertex<3>*>(element->BackLeftBottomCorner)->Faces.push_back(leftBoundary);
				static_cast<MeshVertex<3>*>(element->FrontLeftBottomCorner)->Faces.push_back(leftBoundary);

				// Right boundary
				element = dynamic_cast<ParallelepipedElement*>(this->Elements[index(ix, ny - 1, iz)]);
				this->_rectangularFaces.emplace_back(numberInterface++, element->BackRightBottomCorner, element->FrontRightBottomCorner, element->BackRightTopCorner, element->FrontRightTopCorner, element, CartesianShapeOrientation::InXOZ);
				RectangularFace* rightBoundary = &this->_rectangularFaces.back();
				this->Faces.push_back(rightBoundary);
				this->BoundaryFaces.push_back(rightBoundary);
				element->SetRightFace(rightBoundary);

				static_cast<MeshVertex<3>*>(element->BackRightTopCorner)->Faces.push_back(rightBoundary);
				static_cast<MeshVertex<3>*>(element->FrontRightTopCorner)->Faces.push_back(rightBoundary);
				static_cast<MeshVertex<3>*>(element->BackRightBottomCorner)->Faces.push_back(rightBoundary);
				static_cast<MeshVertex<3>*>(element->FrontRightBottomCorner)->Faces.push_back(rightBoundary);
			}
		}

		for (BigNumber iz = 0; iz < nz; ++iz)
		{
			for (BigNumber iy = 0; iy < ny; ++iy)
			{
				// Back boundary
				ParallelepipedElement* element = dynamic_cast<ParallelepipedElement*>(this->Elements[index(0, iy, iz)]);
				this->_rectangularFaces.emplace_back(numberInterface++, element->BackLeftBottomCorner, element->BackRightBottomCorner, element->BackLeftTopCorner, element->BackRightTopCorner, element, CartesianShapeOrientation::InYOZ);
				RectangularFace* backBoundary = &this->_rectangularFaces.back();
				this->Faces.push_back(backBoundary);
				this->BoundaryFaces.push_back(backBoundary);
				element->SetBackFace(backBoundary);

				static_cast<MeshVertex<3>*>(element->BackLeftTopCorner)->Faces.push_back(backBoundary);
				static_cast<MeshVertex<3>*>(element->BackLeftBottomCorner)->Faces.push_back(backBoundary);
				static_cast<MeshVertex<3>*>(element->BackRightBottomCorner)->Faces.push_back(backBoundary);
				static_cast<MeshVertex<3>*>(element->BackRightTopCorner)->Faces.push_back(backBoundary);

				// Front boundary
				element = dynamic_cast<ParallelepipedElement*>(this->Elements[index(nx - 1, iy, iz)]);
				this->_rectangularFaces.emplace_back(numberInterface++, element->FrontLeftBottomCorner, element->FrontRightBottomCorner, element->FrontLeftTopCorner, element->FrontRightTopCorner, element, CartesianShapeOrientation::InYOZ);
				RectangularFace* frontBoundary = &this->_rectangularFaces.back();
				this->Faces.push_back(frontBoundary);
				this->BoundaryFaces.push_back(frontBoundary);
				element->SetFrontFace(frontBoundary);

				static_cast<MeshVertex<3>*>(element->FrontLeftTopCorner)->Faces.push_back(frontBoundary);
				static_cast<MeshVertex<3>*>(element->FrontLeftBottomCorner)->Faces.push_back(frontBoundary);
				static_cast<MeshVertex<3>*>(element->FrontRightBottomCorner)->Faces.push_back(frontBoundary);
				static_cast<MeshVertex<3>*>(element->FrontRightTopCorner)->Faces.push_back(frontBoundary);
			}
		}

		for (BigNumber iz = 0; iz < nz; ++iz)
		{
			for (BigNumber iy = 0; iy < ny; ++iy)
			{
				for (BigNumber ix = 0; ix < nx; ++ix)
				{
					ParallelepipedElement* element = dynamic_cast<ParallelepipedElement*>(this->Elements[index(ix, iy, iz)]);
					if (ix != nx - 1)
					{
						// Front
						ParallelepipedElement* frontNeighbour = dynamic_cast<ParallelepipedElement*>(this->Elements[index(ix+1, iy, iz)]);
						this->_rectangularFaces.emplace_back(numberInterface++, element->FrontLeftBottomCorner, element->FrontRightBottomCorner, element->FrontLeftTopCorner, element->FrontRightTopCorner, element, frontNeighbour, CartesianShapeOrientation::InYOZ);
						RectangularFace* interface = &this->_rectangularFaces.back();
						this->Faces.push_back(interface);
						this->InteriorFaces.push_back(interface);
						element->SetFrontFace(interface);
						frontNeighbour->SetBackFace(interface);

						static_cast<MeshVertex<3>*>(element->FrontLeftTopCorner)->Faces.push_back(interface);
						static_cast<MeshVertex<3>*>(element->FrontLeftBottomCorner)->Faces.push_back(interface);
						static_cast<MeshVertex<3>*>(element->FrontRightBottomCorner)->Faces.push_back(interface);
						static_cast<MeshVertex<3>*>(element->FrontRightTopCorner)->Faces.push_back(interface);
					}
					if (iy != ny - 1)
					{
						// Right
						ParallelepipedElement* rightNeighbour = dynamic_cast<ParallelepipedElement*>(this->Elements[index(ix, iy+1, iz)]);
						this->_rectangularFaces.emplace_back(numberInterface++, element->BackRightBottomCorner, element->FrontRightBottomCorner, element->BackRightTopCorner, element->FrontRightTopCorner, element, rightNeighbour, CartesianShapeOrientation::InXOZ);
						RectangularFace* interface = &this->_rectangularFaces.back();
						this->Faces.push_back(interface);
						this->InteriorFaces.push_back(interface);
						element->SetRightFace(interface);
						rightNeighbour->SetLeftFace(interface);

						static_cast<MeshVertex<3>*>(element->FrontRightTopCorner)->Faces.push_back(interface);
						static_cast<MeshVertex<3>*>(element->FrontRightBottomCorner)->Faces.push_back(interface);
						static_cast<MeshVertex<3>*>(element->BackRightBottomCorner)->Faces.push_back(interface);
						static_cast<MeshVertex<3>*>(element->BackRightTopCorner)->Faces.push_back(interface);
					}
					if (iz != nz - 1)
					{
						// Top
						ParallelepipedElement* topNeighbour = dynamic_cast<ParallelepipedElement*>(this->Elements[index(ix, iy, iz+1)]);
						this->_rectangularFaces.emplace_back(numberInterface++, element->BackLeftTopCorner, element->FrontLeftTopCorner, element->BackRightTopCorner, element->FrontRightTopCorner, element, topNeighbour, CartesianShapeOrientation::InXOY);
						RectangularFace* interface = &this->_rectangularFaces.back();
						this->Faces.push_back(interface);
						this->InteriorFaces.push_back(interface);
						element->SetTopFace(interface);
						topNeighbour->SetBottomFace(interface);

						static_cast<MeshVertex<3>*>(element->FrontRightTopCorner)->Faces.push_back(interface);
						static_cast<MeshVertex<3>*>(element->FrontLeftTopCorner)->Faces.push_back(interface);
						static_cast<MeshVertex<3>*>(element->BackRightTopCorner)->Faces.push_back(interface);
						static_cast<MeshVertex<3>*>(element->BackLeftTopCorner)->Faces.push_back(interface);
					}
				}
			}
		}

	}

public:
	string Description() override
	{
		return "Cartesian " + to_string(this->Nx) + " x " + to_string(this->Ny) + " x " + to_string(this->Nz);
	}
	string FileNamePart() override
	{
		return "cube-cart-n" + to_string(this->Nx);
	}
	string GeometryDescription() override
	{
		return "Cube";
	}

	double H() override
	{
		return 1.0 / this->Nx;
	}

	double Regularity() override
	{
		return 1;
	}

	double AverageH() override
	{
		return H();
	}

	void CoarsenMesh(H_CoarsStgy elemCoarseningStgy, FaceCoarseningStrategy faceCoarseningStgy, double coarseningFactor)
	{
		if (elemCoarseningStgy == H_CoarsStgy::StandardCoarsening)
		{
			if (faceCoarseningStgy == FaceCoarseningStrategy::InterfaceCollapsing)
				StandardCoarsening();
			else
				Utils::FatalError("Unmanaged face coarsening strategy");
		}
		else
			Mesh<3>::CoarsenMesh(elemCoarseningStgy, faceCoarseningStgy, coarseningFactor);
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
			Cube_CartesianMesh* coarseMesh = new Cube_CartesianMesh(nx / 2, ny / 2, nz / 2, false);
			this->InitializeCoarsening(coarseMesh);
			coarseMesh->ComesFrom.CS = H_CoarsStgy::StandardCoarsening;
			coarseMesh->ComesFrom.nFineElementsByCoarseElement = 8;
			coarseMesh->ComesFrom.nFineFacesAddedByCoarseElement = 12;
			coarseMesh->ComesFrom.nFineFacesByKeptCoarseFace = 4;
			coarseMesh->Build();

			for (BigNumber iz = 0; iz < nz; ++iz)
			{
				for (BigNumber iy = 0; iy < ny; ++iy)
				{
					for (BigNumber ix = 0; ix < nx; ++ix)
					{
						ParallelepipedElement* fineElement = dynamic_cast<ParallelepipedElement*>(this->Elements[index(ix, iy, iz)]);
						ParallelepipedElement* coarseElement = dynamic_cast<ParallelepipedElement*>(coarseMesh->Elements[coarseMesh->index(ix / 2, iy / 2, iz / 2)]);

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
			this->FinalizeCoarsening();
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

public:
	Mesh<3>* Copy() override
	{
		return new Cube_CartesianMesh(this->Nx, this->Ny, this->Nz);
	}
};