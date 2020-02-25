#pragma once
#include "Polygon.h"
#include "../PolyhedralMesh.h"
using namespace std;

class QuadrilateralAsPolygonalMesh : public PolyhedralMesh<2>
{
public:
	BigNumber Nx;
	BigNumber Ny;

	QuadrilateralAsPolygonalMesh(BigNumber nx, BigNumber ny, double xShiftAsFraction) : PolyhedralMesh()
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
				vector<Vertex*> vertices = { bottomLeftCorner, bottomRightCorner, topRightCorner, topLeftCorner };
				Polygon* p = new Polygon(number, vertices);
				this->Elements.push_back(p);
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
			Polygon* p = dynamic_cast<Polygon*>(this->Elements[index(ix, 0)]);
			Edge* southBoundary = new Edge(numberInterface++, p->Vertices()[0], p->Vertices()[1], p);
			this->Faces.push_back(southBoundary);
			this->BoundaryFaces.push_back(southBoundary);
			p->AddFace(southBoundary);

			// North boundary
			p = dynamic_cast<Polygon*>(this->Elements[index(ix, ny - 1)]);
			Edge* northBoundary = new Edge(numberInterface++, p->Vertices()[3], p->Vertices()[2], p);
			this->Faces.push_back(northBoundary);
			this->BoundaryFaces.push_back(northBoundary);
			p->AddFace(northBoundary);
		}

		for (BigNumber iy = 0; iy < ny; ++iy)
		{
			// West boundary
			Polygon* p = dynamic_cast<Polygon*>(this->Elements[index(0, iy)]);
			Edge* westBoundary = new Edge(numberInterface++, p->Vertices()[0], p->Vertices()[3], p);
			this->Faces.push_back(westBoundary);
			this->BoundaryFaces.push_back(westBoundary);
			p->AddFace(westBoundary);

			// East boundary
			p = dynamic_cast<Polygon*>(this->Elements[index(nx-1, iy)]);
			Edge* eastBoundary = new Edge(numberInterface++, p->Vertices()[1], p->Vertices()[2], p);
			this->Faces.push_back(eastBoundary);
			this->BoundaryFaces.push_back(eastBoundary);
			p->AddFace(eastBoundary);
		}

		for (BigNumber iy = 0; iy < ny; iy++)
		{
			for (BigNumber ix = 0; ix < nx; ix++)
			{
				Polygon* element = dynamic_cast<Polygon*>(this->Elements[index(ix, iy)]);
				if (ix != nx - 1)
				{
					// East
					Polygon* eastNeighbour = dynamic_cast<Polygon*>(this->Elements[index(ix + 1, iy)]);
					Edge* interface = new Edge(numberInterface++, eastNeighbour->Vertices()[0], eastNeighbour->Vertices()[3], element, eastNeighbour);
					this->Faces.push_back(interface);
					this->InteriorFaces.push_back(interface);
					element->AddFace(interface);
					eastNeighbour->AddFace(interface);
				}
				if (iy != ny - 1)
				{
					// North
					Polygon* northNeighbour = dynamic_cast<Polygon*>(this->Elements[index(ix, iy + 1)]);
					Edge* interface = new Edge(numberInterface++, northNeighbour->Vertices()[0], northNeighbour->Vertices()[1], element, northNeighbour);
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
	string Description() override
	{
		return "Quadrilaterals as polygons " + to_string(this->Nx) + " x " + to_string(this->Ny);
	}

	string FileNamePart() override
	{
		return "n" + to_string(this->Nx);
	}

	double H() override
	{
		return this->Elements[1]->Diameter();
	}

	double Regularity() override
	{
		return this->Elements[1]->Regularity();
	}

};