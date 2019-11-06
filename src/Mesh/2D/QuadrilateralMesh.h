#pragma once
#include "Quadrilateral.h"
#include "../Mesh.h"
using namespace std;

class QuadrilateralMesh : public Mesh<2>
{
public:
	BigNumber Nx;
	BigNumber Ny;

	QuadrilateralMesh(BigNumber nx, BigNumber ny, double xShiftAsFraction) : Mesh()
	{
		// nx = ny falls down to square elements
		this->Nx = nx;
		this->Ny = ny;

		double hx = 1.0 / nx;
		double hy = 1.0 / ny;

		assert(xShiftAsFraction >= 0 && xShiftAsFraction < 1);
		double xShift = xShiftAsFraction * hx;

		//----------//
		// Vertices //
		//----------//

		this->Vertices.reserve((nx + 1) * (ny + 1));
		for (BigNumber iy = 0; iy < ny + 1; ++iy)
		{
			for (BigNumber ix = 0; ix < nx + 1; ++ix)
			{
				BigNumber number = indexV(ix, iy);
				double x = ix * hx;
				double y = iy * hy;
				if (ix != 0 && ix != nx && iy % 2 == 1)
					x += xShift;
				Vertex* vertex = new Vertex(number, x, y);
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
				Quadrilateral* quad = new Quadrilateral(number, bottomLeftCorner, bottomRightCorner, topRightCorner, topLeftCorner);
				this->Elements.push_back(quad);
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
			Quadrilateral* quad = dynamic_cast<Quadrilateral*>(this->Elements[index(ix, 0)]);
			Edge* southBoundary = new Edge(numberInterface++, quad->V1(), quad->V2(), quad);
			this->Faces.push_back(southBoundary);
			this->BoundaryFaces.push_back(southBoundary);
			quad->AddFace(southBoundary);

			// North boundary
			quad = dynamic_cast<Quadrilateral*>(this->Elements[index(ix, ny - 1)]);
			Edge* northBoundary = new Edge(numberInterface++, quad->V4(), quad->V3(), quad);
			this->Faces.push_back(northBoundary);
			this->BoundaryFaces.push_back(northBoundary);
			quad->AddFace(northBoundary);
		}

		for (BigNumber iy = 0; iy < ny; ++iy)
		{
			// West boundary
			Quadrilateral* quad = dynamic_cast<Quadrilateral*>(this->Elements[index(0, iy)]);
			Edge* westBoundary = new Edge(numberInterface++, quad->V1(), quad->V4(), quad);
			this->Faces.push_back(westBoundary);
			this->BoundaryFaces.push_back(westBoundary);
			quad->AddFace(westBoundary);

			// East boundary
			quad = dynamic_cast<Quadrilateral*>(this->Elements[index(nx-1, iy)]);
			Edge* eastBoundary = new Edge(numberInterface++, quad->V2(), quad->V3(), quad);
			this->Faces.push_back(eastBoundary);
			this->BoundaryFaces.push_back(eastBoundary);
			quad->AddFace(eastBoundary);
		}

		for (BigNumber iy = 0; iy < ny; iy++)
		{
			for (BigNumber ix = 0; ix < nx; ix++)
			{
				Quadrilateral* element = dynamic_cast<Quadrilateral*>(this->Elements[index(ix, iy)]);
				if (ix != nx - 1)
				{
					// East
					Quadrilateral* eastNeighbour = dynamic_cast<Quadrilateral*>(this->Elements[index(ix + 1, iy)]);
					Edge* interface = new Edge(numberInterface++, eastNeighbour->V1(), eastNeighbour->V4(), element, eastNeighbour);
					this->Faces.push_back(interface);
					this->InteriorFaces.push_back(interface);
					element->AddFace(interface);
					eastNeighbour->AddFace(interface);
				}
				if (iy != ny - 1)
				{
					// North
					Quadrilateral* northNeighbour = dynamic_cast<Quadrilateral*>(this->Elements[index(ix, iy + 1)]);
					Edge* interface = new Edge(numberInterface++, northNeighbour->V1(), northNeighbour->V2(), element, northNeighbour);
					this->Faces.push_back(interface);
					this->InteriorFaces.push_back(interface);
					element->AddFace(interface);
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
		return y * Nx + x;
	}

public:
	string Description()
	{
		return "Quadrilateral " + to_string(this->Nx) + " x " + to_string(this->Ny);
	}

	string FileNamePart()
	{
		return "n" + to_string(this->Nx);
	}

	double H()
	{
		return 1.0 / this->Nx;
	}

	void CoarsenMesh(CoarseningStrategy strategy)
	{
		/*if (strategy == CoarseningStrategy::Standard)
			CoarsenByAgglomerationAndMergeColinearFaces();
		//else if (strategy == CoarseningStrategy::Agglomeration)
			//CoarsenByAgglomerationAndKeepFineFaces();
		else*/
			assert(false && "Coarsening strategy not implemented!");
		this->CoarseMesh->SetDiffusionCoefficient(this->_diffusionPartition);
		this->CoarseMesh->SetBoundaryConditions(this->_boundaryConditions);
	}

};