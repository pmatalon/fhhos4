#pragma once
#include <vector>
#include "Triangle.h"
#include "../PolyhedralMesh.h"
using namespace std;

class LowerTriangle : public Triangle
{
public:
	LowerTriangle(int number, Vertex* v1, Vertex* v2, Vertex* v3) : Triangle(number, v1, v2, v3), Element(number) {}

	Edge* GetEdge(Vertex* v1, Vertex* v2)
	{
		for (Face<2>* f : this->Faces)
		{
			Edge* e = dynamic_cast<Edge*>(f);
			if ((e->Vertex1() == v1 && e->Vertex2() == v2) || (e->Vertex1() == v2 && e->Vertex2() == v1))
				return e;
		}
		assert(false);
	}

	Edge* SouthEdge() { return GetEdge(this->V1(), this->V2());	}
	Edge* ObliqueEdge() { return GetEdge(this->V2(), this->V3()); }
	Edge* WestEdge() { return GetEdge(this->V3(), this->V1()); }
};

class UpperTriangle : public Triangle
{
public:
	UpperTriangle(int number, Vertex* v1, Vertex* v2, Vertex* v3) : Triangle(number, v1, v2, v3), Element(number) {}

	Edge* GetEdge(Vertex* v1, Vertex* v2)
	{
		for (Face<2>* f : this->Faces)
		{
			Edge* e = dynamic_cast<Edge*>(f);
			if ((e->Vertex1() == v1 && e->Vertex2() == v2) || (e->Vertex1() == v2 && e->Vertex2() == v1))
				return e;
		}
		assert(false);
	}

	Edge* ObliqueEdge() { return GetEdge(this->V1(), this->V2()); }
	Edge* EastEdge() { return GetEdge(this->V2(), this->V3()); }
	Edge* NorthEdge() { return GetEdge(this->V3(), this->V1()); }
};

class TriangularMesh : public PolyhedralMesh<2>
{
public:
	BigNumber Nx;
	BigNumber Ny;

	TriangularMesh(BigNumber nx, BigNumber ny) : PolyhedralMesh()
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

		this->Elements.reserve(2 * nx * ny);
		for (BigNumber iy = 0; iy < ny; ++iy)
		{
			for (BigNumber ix = 0; ix < nx; ++ix)
			{
				BigNumber number = index(ix, iy);
				Vertex* bottomLeftCorner  = Vertices[indexV(ix,     iy    )];
				Vertex* topLeftCorner     = Vertices[indexV(ix,     iy + 1)];
				Vertex* topRightCorner    = Vertices[indexV(ix + 1, iy + 1)];
				Vertex* bottomRightCorner = Vertices[indexV(ix + 1, iy    )];
				Triangle* lowerTriangle = new LowerTriangle(number, bottomLeftCorner, bottomRightCorner, topLeftCorner);
				this->Elements.push_back(lowerTriangle);
				Triangle* upperTriangle = new UpperTriangle(number + 1, topLeftCorner, bottomRightCorner, topRightCorner);
				this->Elements.push_back(upperTriangle);
			}
		}

		//-------//
		// Faces //
		//-------//

		this->Faces.reserve(nx * (ny + 1) + ny * (nx + 1) + nx * ny);
		BigNumber numberInterface = 0;

		for (BigNumber ix = 0; ix < nx; ++ix)
		{
			// South boundary
			Triangle* lowerTriangle = dynamic_cast<Triangle*>(this->Elements[index(ix, 0)]);
			Edge* southBoundary = new Edge(numberInterface++, lowerTriangle->V1(), lowerTriangle->V2(), lowerTriangle);
			this->Faces.push_back(southBoundary);
			this->BoundaryFaces.push_back(southBoundary);
			lowerTriangle->AddFace(southBoundary);

			// North boundary
			Triangle* upperTriangle = dynamic_cast<Triangle*>(this->Elements[index(ix, ny - 1) + 1]);
			Edge* northBoundary = new Edge(numberInterface++, upperTriangle->V1(), upperTriangle->V3(), upperTriangle);
			this->Faces.push_back(northBoundary);
			this->BoundaryFaces.push_back(northBoundary);
			upperTriangle->AddFace(northBoundary);
		}

		for (BigNumber iy = 0; iy < ny; ++iy)
		{
			// West boundary
			Triangle* lowerTriangle = dynamic_cast<Triangle*>(this->Elements[index(0, iy)]);
			Edge* westBoundary = new Edge(numberInterface++, lowerTriangle->V3(), lowerTriangle->V1(), lowerTriangle);
			this->Faces.push_back(westBoundary);
			this->BoundaryFaces.push_back(westBoundary);
			lowerTriangle->AddFace(westBoundary);

			// East boundary
			Triangle* upperTriangle = dynamic_cast<Triangle*>(this->Elements[index(nx-1, iy)+1]);
			Edge* eastBoundary = new Edge(numberInterface++, upperTriangle->V2(), upperTriangle->V3(), upperTriangle);
			this->Faces.push_back(eastBoundary);
			this->BoundaryFaces.push_back(eastBoundary);
			upperTriangle->AddFace(eastBoundary);
		}

		for (BigNumber iy = 0; iy < ny; iy++)
		{
			for (BigNumber ix = 0; ix < nx; ix++)
			{
				Triangle* lowerTriangle = dynamic_cast<Triangle*>(this->Elements[index(ix, iy)]);
				Triangle* upperTriangle = dynamic_cast<Triangle*>(this->Elements[index(ix, iy)+1]);

				Edge* interface = new Edge(numberInterface++, lowerTriangle->V3(), lowerTriangle->V2(), lowerTriangle, upperTriangle);
				this->Faces.push_back(interface);
				this->InteriorFaces.push_back(interface);
				lowerTriangle->AddFace(interface);
				upperTriangle->AddFace(interface);

				if (ix != nx - 1)
				{
					// East
					Triangle* eastNeighbour = dynamic_cast<Triangle*>(this->Elements[index(ix + 1, iy)]);
					interface = new Edge(numberInterface++, eastNeighbour->V3(), eastNeighbour->V1(), upperTriangle, eastNeighbour);
					this->Faces.push_back(interface);
					this->InteriorFaces.push_back(interface);
					upperTriangle->AddFace(interface);
					eastNeighbour->AddFace(interface);
				}
				if (iy != ny - 1)
				{
					// North
					Triangle* northNeighbour = dynamic_cast<Triangle*>(this->Elements[index(ix, iy + 1)]);
					interface = new Edge(numberInterface++, northNeighbour->V1(), northNeighbour->V2(), upperTriangle, northNeighbour);
					this->Faces.push_back(interface);
					this->InteriorFaces.push_back(interface);
					upperTriangle->AddFace(interface);
					northNeighbour->AddFace(interface);
				}
			}
		}

	}

	void SanityCheck() override
	{
		Mesh::SanityCheck();
		RefFunction refX = [](RefPoint p) { return p.X; };
		DomFunction domX = [](DomPoint p) { return p.X; };

		for (auto e : this->Elements)
		{
			Triangle* t = dynamic_cast<Triangle*>(e);
			assert(t->Faces.size() == 3);

			//assert(t->Shape()->DetJacobian() == 1/t->Shape()->InverseJacobianTranspose().determinant());

			for (auto f : t->Faces)
			{
				auto n = t->OuterNormalVector(f);
				Edge* edge = dynamic_cast<Edge*>(f);
				assert(n.dot(Vect<2>(edge->Vertex1(), edge->Vertex2())) == 0);
			}
			
			assert(t->ConvertToDomain(t->ConvertToReference(*(t->V1()))) == *(t->V1()));
			assert(t->ConvertToDomain(t->ConvertToReference(*(t->V2()))) == *(t->V2()));
			assert(t->ConvertToDomain(t->ConvertToReference(*(t->V3()))) == *(t->V3()));

			RefPoint ref1 = t->ConvertToReference(*(t->V1()));
			assert(ref1 == RefPoint(0, 0) || ref1 == RefPoint(1, 0) || ref1 == RefPoint(0, 1));
			RefPoint ref2 = t->ConvertToReference(*(t->V2()));
			assert((ref2 == RefPoint(0, 0) || ref2 == RefPoint(1, 0) || ref2 == RefPoint(0, 1)) && ref2 != ref1);
			RefPoint ref3 = t->ConvertToReference(*(t->V3()));
			assert((ref3 == RefPoint(0, 0) || ref3 == RefPoint(1, 0) || ref3 == RefPoint(0, 1)) && ref3 != ref1 && ref3 != ref2);
		}

		ReferenceTriangle refTriangle;
		for (int degree = 1; degree < 5; degree++)
		{
			double integral = refTriangle.Integral(refX, degree);
			assert(abs(integral - 1.0 / 6.0) < 1e-15);
		}

		BigNumber number = 0;
		Vertex bottomLeft(number,  0, 0);
		Vertex bottomRight(number, 1, 0);
		Vertex topRight(number,    1, 1);
		Vertex topLeft(number,     0, 1);

		Triangle lower(number, &bottomLeft, &bottomRight, &topLeft);
		Triangle upper(number, &topLeft, &bottomRight, &topRight);

		/*double lowerIntegral = lower.Integral(domX, 1);
		assert(abs(lowerIntegral - 1.0 / 6) < 1e-14);
		double upperIntegral = upper.Integral(domX, 1);
		assert(abs(upperIntegral - 1.0 / 3) < 1e-14);*/
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
	string Description() override
	{
		return "Triangular " + to_string(this->Nx) + " x " + to_string(this->Ny);
	}

	string FileNamePart() override
	{
		return "tri-n" + to_string(this->Nx);
	}

	double H() override
	{
		return sqrt(2) / (double)this->Nx;
	}

	void CoarsenMesh(CoarseningStrategy strategy) override
	{
		if (strategy == CoarseningStrategy::Standard)
			Standard();
		//else if (strategy == CoarseningStrategy::Agglomeration)
			//CoarsenByAgglomerationAndKeepFineFaces();
		else
			assert(false && "Coarsening strategy not implemented!");
		this->CoarseMesh->SetDiffusionCoefficient(this->_diffusionPartition);
		this->CoarseMesh->SetBoundaryConditions(this->_boundaryConditions);
	}

	void Standard()
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
			coarseMesh->ComesFrom.CS = CoarseningStrategy::Standard;
			coarseMesh->ComesFrom.nFineElementsByCoarseElement = 4;
			coarseMesh->ComesFrom.nFineFacesAddedByCoarseElement = 3;
			coarseMesh->ComesFrom.nFineFacesByKeptCoarseFace = 2;

			for (BigNumber i = 0; i < ny; ++i)
			{
				for (BigNumber j = 0; j < nx; ++j)
				{
					LowerTriangle* fineLower = dynamic_cast<LowerTriangle*>(this->Elements[index(i, j)]);
					UpperTriangle* fineUpper = dynamic_cast<UpperTriangle*>(this->Elements[index(i, j)+1]);

					LowerTriangle* coarseLower = dynamic_cast<LowerTriangle*>(coarseMesh->Elements[coarseMesh->index(i / 2, j / 2)]);
					UpperTriangle* coarseUpper = dynamic_cast<UpperTriangle*>(coarseMesh->Elements[coarseMesh->index(i / 2, j / 2)+1]);

					if (j % 2 == 0) // on an even row
					{
						coarseLower->FinerElements.push_back(fineLower);
						fineLower->CoarserElement = coarseLower;
						//Edge* coarseSouthEdge = coarseLower->SouthEdge();
						//cout << "south edge = " << *coarseSouthEdge << endl;
						coarseLower->SouthEdge()->FinerFaces.push_back(fineLower->SouthEdge());

						if (i % 2 == 0) // on an even column
						{
							coarseLower->FinerElements.push_back(fineUpper);
							fineUpper->CoarserElement = coarseLower;
							
							fineUpper->ObliqueEdge()->IsRemovedOnCoarserGrid = true;
							fineUpper->EastEdge()->IsRemovedOnCoarserGrid = true;
							fineUpper->NorthEdge()->IsRemovedOnCoarserGrid = true;
							coarseLower->FinerFacesRemoved.push_back(fineUpper->ObliqueEdge());
							coarseLower->FinerFacesRemoved.push_back(fineUpper->EastEdge());
							coarseLower->FinerFacesRemoved.push_back(fineUpper->NorthEdge());

							coarseLower->WestEdge()->FinerFaces.push_back(fineLower->WestEdge());
						}
						else
						{
							coarseUpper->FinerElements.push_back(fineUpper);
							fineUpper->CoarserElement = coarseUpper;

							coarseLower->ObliqueEdge()->FinerFaces.push_back(fineLower->ObliqueEdge());
						}
					}
					else // on an odd row
					{
						coarseUpper->FinerElements.push_back(fineUpper);
						fineUpper->CoarserElement = coarseUpper;
						if (i % 2 == 0) // on an even col
						{
							coarseLower->FinerElements.push_back(fineLower);
							fineLower->CoarserElement = coarseLower;

							coarseLower->WestEdge()->FinerFaces.push_back(fineLower->WestEdge());
							coarseLower->ObliqueEdge()->FinerFaces.push_back(fineLower->ObliqueEdge());
						}
						else // on an odd col
						{
							coarseUpper->FinerElements.push_back(fineLower);
							fineLower->CoarserElement = coarseUpper;

							fineLower->SouthEdge()->IsRemovedOnCoarserGrid = true;
							fineLower->ObliqueEdge()->IsRemovedOnCoarserGrid = true;
							fineLower->WestEdge()->IsRemovedOnCoarserGrid = true;
							coarseUpper->FinerFacesRemoved.push_back(fineLower->SouthEdge());
							coarseUpper->FinerFacesRemoved.push_back(fineLower->ObliqueEdge());
							coarseUpper->FinerFacesRemoved.push_back(fineLower->WestEdge());
						}
					}
					if (j == nx - 1) // last row
						coarseUpper->NorthEdge()->FinerFaces.push_back(fineUpper->NorthEdge());
					if (i == ny - 1) // last col
						coarseUpper->EastEdge()->FinerFaces.push_back(fineUpper->EastEdge());
				}
			}

			this->CoarseMesh = coarseMesh;
			//SanityCheckStandardCoarsening();
		}
	}

	void SanityCheckStandardCoarsening()
	{
		for (Element<2>* ce : this->CoarseMesh->Elements)
		{
			assert(ce->FinerElements.size() == 4);
			assert(ce->FinerFacesRemoved.size() == 3);
		}

		for (Face<2>* cf : this->CoarseMesh->Faces)
			assert(cf->FinerFaces.size() == 2);
	}

};