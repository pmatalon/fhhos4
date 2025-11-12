#pragma once
#include <vector>
#include "../../Mesh/2D/RectangularElement.h"
#include "../../Mesh/2D/RectangularPolygonalElement.h"
#include "Square_CartesianPolygonalMesh.h"
#include "../../Mesh/2D/CartesianEdge.h"
#include "../../Mesh/PolyhedralMesh.h"
#include "../SquareGeometry.h"
#include "../Square4quadrantsGeometry.h"
using namespace std;

class Square_CartesianPolyStripeMesh : public PolyhedralMesh<2>
{
public:
	BigNumber Nx_l; //high reso stat
	BigNumber Ny_l; //low reso stat

	BigNumber Nx_r; //high reso stat
	BigNumber Ny_r; //low reso stat
	bool With4Quadrants;

	// changeable parameters
	int ky = 16; // factor relating high and low reso along y-axis
	int n_breaks; //= 2;
	double r_length = 0.1;

	// derived quantities

	int n_r_stripes; //= 1 + n_breaks / 2;
	int n_l_stripes; //= (n_breaks + 1) / 2;
	double l_length; //= (1.0 - n_r_stripes * r_length) / n_l_stripes;
	int n_stripes; //= n_breaks + 1;

	vector<double> h_breaks = {0.0};
	vector<int> r_indices = {};
	vector<int> l_indices = {};

	vector<BigNumber> stripe_vertex_breaks = {0};
	vector<BigNumber> stripe_element_breaks = {0};

	Square_CartesianPolyStripeMesh(BigNumber nx, BigNumber ny, const int nb_stripes, bool with4Quadrants = false, bool buildMesh = true) : PolyhedralMesh()
	{
		this->n_stripes = nb_stripes;
		this->n_breaks = nb_stripes - 1;
		this->n_r_stripes = 1 + n_breaks / 2;
		this->n_l_stripes = (n_breaks + 1) / 2;
		this->l_length = (1.0 - n_r_stripes * r_length) / n_l_stripes;

		assert(nb_stripes > 0);
		assert(l_length > 0.0);

		for (int i = 0; i < n_stripes; i++)
		{
			if (i % 2 == 0)
			{
				// r-stripe
				h_breaks.push_back(h_breaks.back() + r_length);
				r_indices.push_back(i);
			}
			else
			{
				// l-stripe
				h_breaks.push_back(h_breaks.back() + l_length);
				l_indices.push_back(i);
			}
		}

		this->Nx_l = nx;
		this->Ny_l = ny;

		this->Nx_r = Nx_l;
		this->Ny_r = ky * ny;
		this->With4Quadrants = with4Quadrants;

		for (int stripe_counter = 0; stripe_counter < n_stripes; stripe_counter++)
		{
			if (stripe_counter % 2 == 0)
			{
				// r-stripe
				stripe_vertex_breaks.push_back(stripe_vertex_breaks.back() + (Nx_r + 1) * (Ny_r + 1));
				stripe_element_breaks.push_back(stripe_element_breaks.back() + Nx_r * Ny_r);
			}
			else
			{
				// l-stripe
				if (stripe_counter == n_stripes - 1)
				{
					// include vertices on east boundary
					stripe_vertex_breaks.push_back(stripe_vertex_breaks.back() + Nx_l * (Ny_l + 1));
				}
				else
				{
					stripe_vertex_breaks.push_back(stripe_vertex_breaks.back() + (Nx_l - 1) * (Ny_l + 1));
				}

				stripe_element_breaks.push_back(stripe_element_breaks.back() + Nx_l * Ny_l);
			}
		}

		if (with4Quadrants && (nx % 2 == 1 || ny % 2 == 1))
			Utils::FatalError("Building the mesh for a square with 4 quadrants requires the number of subdivisions in each direction to be even.");

		if (buildMesh)
			Build();
	}

	void Build()
	{
		BigNumber nx_l = this->Nx_l;
		BigNumber ny_l = this->Ny_l;
		BigNumber nx_r = this->Nx_r;
		BigNumber ny_r = this->Ny_r;

		double hy_l = 1.0 / ny_l;
		double hx_l = l_length / nx_l;
		double hy_r = 1.0 / ny_r;
		double hx_r = r_length / nx_r;

		// Physical parts
		PhysicalGroup<2>* domain = nullptr;
		PhysicalGroup<2>* quadrantBottomLeft = nullptr;
		PhysicalGroup<2>* quadrantBottomRight = nullptr;
		PhysicalGroup<2>* quadrantTopRight = nullptr;
		PhysicalGroup<2>* quadrantTopLeft = nullptr;

		if (this->With4Quadrants)
		{
			Utils::FatalError("Hetereogenity for Emil-mesh not yet implemented");
			if (this->PhysicalParts.empty())
				this->PhysicalParts = Square4quadrantsGeometry::PhysicalParts();
			quadrantBottomLeft = this->PhysicalParts[0];
			quadrantBottomRight = this->PhysicalParts[1];
			quadrantTopRight = this->PhysicalParts[2];
			quadrantTopLeft = this->PhysicalParts[3];
		}
		else
		{
			if (this->PhysicalParts.empty())
				this->PhysicalParts = SquareGeometry::PhysicalParts();
			domain = this->PhysicalParts[0];
		}

		// Boundary parts
		if (this->BoundaryParts.empty())
			this->BoundaryParts = SquareGeometry::BoundaryParts();
		BoundaryGroup* squareBottomBoundary = this->BoundaryParts[0];
		BoundaryGroup* squareRightBoundary = this->BoundaryParts[1];
		BoundaryGroup* squareTopBoundary = this->BoundaryParts[2];
		BoundaryGroup* squareLeftBoundary = this->BoundaryParts[3];

		//----------//
		// Vertices //
		//----------//

		//this->Vertices.reserve(stripe_vertex_breaks.back());
		this->Vertices.resize(stripe_vertex_breaks.back());

		cout << "vertices before adding them: " << Vertices.size() << endl;

		for (int r_index : r_indices)
		{
			for (BigNumber iy = 0; iy < ny_r + 1; ++iy)
			{
				for (BigNumber ix = 0; ix < nx_r + 1; ++ix)
				{
					auto vertexIndex = indexV(ix, iy, r_index);
					auto* vertex = new Vertex(vertexIndex, h_breaks[r_index] + ix * hx_r, iy * hy_r);
					//this->Vertices.push_back(vertex);
					this->Vertices[vertexIndex] = vertex;
				}
			}
		}

		for (int l_index : l_indices)
		{
			for (BigNumber iy = 0; iy < ny_l + 1; ++iy)
			{
				for (BigNumber ix = 1; ix < nx_l; ++ix)
				{
					auto vertexIndex = indexV(ix, iy, l_index);
					auto* vertex = new Vertex(vertexIndex, h_breaks[l_index] + ix * hx_l, iy * hy_l);
					//this->Vertices.push_back(vertex);
					this->Vertices[vertexIndex] = vertex;
				}

				if (l_index == n_stripes - 1)
				{
					// include vertices on left boundary
					BigNumber ix = nx_l;
					auto vertexIndex = indexV(ix, iy, l_index);
					auto* vertex = new Vertex(vertexIndex, r_length + ix * hx_l, iy * hy_l);
					//this->Vertices.push_back(vertex);
					this->Vertices[vertexIndex] = vertex;
				}
			}
		}

		if (Vertices.size()!=stripe_vertex_breaks.back())
			cout << "Vertices size vs what should be: " << Vertices.size() << "/" << stripe_vertex_breaks.back() << endl;


		//----------//
		// Elements //
		//----------//

		//this->Elements.reserve( stripe_element_breaks.back());
		this->Elements.resize(stripe_element_breaks.back());

		for (int r_index : r_indices)
		{
			for (BigNumber iy = 0; iy < ny_r; ++iy)
			{
				for (BigNumber ix = 0; ix < nx_r; ++ix)
				{
					Vertex* bottomLeftCorner  = Vertices[indexV(ix, iy, r_index)];
					Vertex* topLeftCorner     = Vertices[indexV(ix, iy + 1, r_index)];
					Vertex* topRightCorner    = Vertices[indexV(ix + 1, iy + 1, r_index)];
					Vertex* bottomRightCorner = Vertices[indexV(ix + 1, iy, r_index)];
					auto eleIndex = indexE(ix, iy, r_index);
					auto* rectangle = new RectangularPolygonalElement(eleIndex, bottomLeftCorner, topLeftCorner, topRightCorner, bottomRightCorner);
					//this->Elements.push_back(rectangle);
					this->Elements[eleIndex] = rectangle;
					rectangle->PhysicalPart = domain;
				}
			}
		}

		for (int l_index : l_indices)
		{
			for (BigNumber iy = 0; iy < ny_l; ++iy)
			{
				for (BigNumber ix = 0; ix < nx_l; ++ix)
				{
					Vertex* bottomLeftCorner  = Vertices[indexV(ix, iy, l_index)];
					Vertex* topLeftCorner     = Vertices[indexV(ix, iy + 1, l_index)];
					Vertex* topRightCorner    = Vertices[indexV(ix + 1, iy + 1, l_index)];
					Vertex* bottomRightCorner = Vertices[indexV(ix + 1, iy, l_index)];
					auto eleIndex = indexE(ix, iy, l_index);
					auto* rectangle = new RectangularPolygonalElement(eleIndex, bottomLeftCorner, topLeftCorner, topRightCorner, bottomRightCorner);
					//this->Elements.push_back(rectangle);
					this->Elements[eleIndex] = rectangle;
					rectangle->PhysicalPart = domain;
				}
			}
		}

		if (Elements.size()!=stripe_element_breaks.back())
			cout << "Elements size vs what should be: " << Elements.size() << "/" << stripe_element_breaks.back() << endl;

		//-------//
		// Faces //
		//-------//

		this->Faces.reserve(nx_l * (ny_l + 1) + ny_l * nx_l +
							  nx_r * (ny_r + 1) + ny_r * (nx_r + 1));
		BigNumber numberInterface = 0;

		for (int r_index : r_indices)
		{
			for (BigNumber ix = 0; ix < Nx_r; ++ix)
			{
				// South boundary
				auto* rectangle = dynamic_cast<RectangularPolygonalElement*>(this->Elements[indexE(ix, 0, r_index)]);
				auto* southBoundary = new CartesianEdge(numberInterface++, rectangle->BottomLeftCorner, rectangle->BottomRightCorner, rectangle, CartesianShapeOrientation::Horizontal);
				this->Faces.push_back(southBoundary);
				this->BoundaryFaces.push_back(southBoundary);
				rectangle->AddSouthFace(southBoundary);
				southBoundary->BoundaryPart = squareBottomBoundary;

				// North boundary
				rectangle = dynamic_cast<RectangularPolygonalElement*>(this->Elements[indexE(ix, ny_r - 1, r_index)]);
				auto* northBoundary = new CartesianEdge(numberInterface++, rectangle->TopLeftCorner, rectangle->TopRightCorner, rectangle, CartesianShapeOrientation::Horizontal);
				this->Faces.push_back(northBoundary);
				this->BoundaryFaces.push_back(northBoundary);
				rectangle->AddNorthFace(northBoundary);
				northBoundary->BoundaryPart = squareTopBoundary;
			}
		}

		for (int l_index : l_indices)
		{
			for (BigNumber ix = 0; ix < Nx_l; ++ix)
			{
				// South boundary
				auto* rectangle = dynamic_cast<RectangularPolygonalElement*>(this->Elements[indexE(ix, 0, l_index)]);
				auto* southBoundary = new CartesianEdge(numberInterface++, rectangle->BottomLeftCorner, rectangle->BottomRightCorner, rectangle, CartesianShapeOrientation::Horizontal);
				this->Faces.push_back(southBoundary);
				this->BoundaryFaces.push_back(southBoundary);
				rectangle->AddSouthFace(southBoundary);
				southBoundary->BoundaryPart = squareBottomBoundary;

				// North boundary
				rectangle = dynamic_cast<RectangularPolygonalElement*>(this->Elements[indexE(ix, ny_l - 1, l_index)]);
				auto* northBoundary = new CartesianEdge(numberInterface++, rectangle->TopLeftCorner, rectangle->TopRightCorner, rectangle, CartesianShapeOrientation::Horizontal);
				this->Faces.push_back(northBoundary);
				this->BoundaryFaces.push_back(northBoundary);
				rectangle->AddNorthFace(northBoundary);
				northBoundary->BoundaryPart = squareTopBoundary;
			}
		}

		// West boundary
		for (BigNumber iy = 0; iy < ny_r; ++iy)
		{
			auto* rectangle = dynamic_cast<RectangularPolygonalElement*>(this->Elements[indexE(0, iy, 0)]);
			auto* westBoundary = new CartesianEdge(numberInterface++, rectangle->BottomLeftCorner, rectangle->TopLeftCorner, rectangle, CartesianShapeOrientation::Vertical);
			this->Faces.push_back(westBoundary);
			this->BoundaryFaces.push_back(westBoundary);
			rectangle->AddWestFace(westBoundary);
			westBoundary->BoundaryPart = squareLeftBoundary;
		}

		// East boundary
		if (l_indices.back() == n_stripes - 1)
		{
			for (BigNumber iy = 0; iy < ny_l; ++iy)
			{
				auto* rectangle = dynamic_cast<RectangularPolygonalElement*>(this->Elements[indexE(nx_l - 1, iy, l_indices.back())]);
				auto* eastBoundary = new CartesianEdge(numberInterface++, rectangle->BottomRightCorner, rectangle->TopRightCorner, rectangle, CartesianShapeOrientation::Vertical);
				this->Faces.push_back(eastBoundary);
				this->BoundaryFaces.push_back(eastBoundary);
				rectangle->AddEastFace(eastBoundary);
				eastBoundary->BoundaryPart = squareRightBoundary;
			}
		}
		else
		{
			for (BigNumber iy = 0; iy < ny_r; ++iy)
			{
				auto* rectangle = dynamic_cast<RectangularPolygonalElement*>(this->Elements[indexE(nx_r - 1, iy, r_indices.back())]);
				auto* eastBoundary = new CartesianEdge(numberInterface++, rectangle->BottomRightCorner, rectangle->TopRightCorner, rectangle, CartesianShapeOrientation::Vertical);
				this->Faces.push_back(eastBoundary);
				this->BoundaryFaces.push_back(eastBoundary);
				rectangle->AddEastFace(eastBoundary);
				eastBoundary->BoundaryPart = squareRightBoundary;
			}
		}

		for (int r_index : r_indices)
		{
			for (BigNumber iy = 0; iy < ny_r; iy++)
			{
				for (BigNumber ix = 0; ix < nx_r; ix++)
				{
					auto* element = dynamic_cast<RectangularPolygonalElement*>(this->Elements[indexE(ix, iy, r_index)]);
					RectangularPolygonalElement* eastNeighbour;
					// If not on east boundary create east face
					if (!(ix == nx_r - 1 && r_index == n_stripes - 1))
					{
						if (ix == nx_r - 1)
						{
							// fetch cell on l-stripe to the east
							eastNeighbour = dynamic_cast<RectangularPolygonalElement*>(this->Elements[indexE(0, iy / ky, r_index + 1)]);
						}
						else
						{
							eastNeighbour = dynamic_cast<RectangularPolygonalElement*>(this->Elements[indexE(ix + 1, iy, r_index)]);
						}

						//auto* interface = new CartesianEdge(numberInterface++, eastNeighbour->BottomLeftCorner, eastNeighbour->TopLeftCorner, element, eastNeighbour, CartesianShapeOrientation::Vertical);
						auto* interface = new CartesianEdge(numberInterface++, element->BottomRightCorner, element->TopRightCorner, element, eastNeighbour, CartesianShapeOrientation::Vertical);
						this->Faces.push_back(interface);
						this->InteriorFaces.push_back(interface);
						element->AddEastFace(interface);
						eastNeighbour->AddWestFace(interface);
					}
					if (iy != ny_r - 1)
					{
						// North
						auto* northNeighbour = dynamic_cast<RectangularPolygonalElement*>(this->Elements[indexE(ix, iy + 1, r_index)]);
						auto* interface = new CartesianEdge(numberInterface++, northNeighbour->BottomLeftCorner, northNeighbour->BottomRightCorner, element, northNeighbour, CartesianShapeOrientation::Horizontal);
						this->Faces.push_back(interface);
						this->InteriorFaces.push_back(interface);
						element->AddNorthFace(interface);
						northNeighbour->AddSouthFace(interface);
					}
				}

				// If not on west boundary create west boundary of the stripe
				if (r_index != 0)
				{
					auto* element = dynamic_cast<RectangularPolygonalElement*>(this->Elements[indexE(0, iy, r_index)]);
					auto* westNeighbour = dynamic_cast<RectangularPolygonalElement*>(this->Elements[indexE(nx_l - 1, iy / ky, r_index - 1)]);
					auto* interface = new CartesianEdge(numberInterface++, element->BottomLeftCorner, element->TopLeftCorner, westNeighbour, element, CartesianShapeOrientation::Vertical);
					this->Faces.push_back(interface);
					this->InteriorFaces.push_back(interface);
					element->AddWestFace(interface);
					westNeighbour->AddEastFace(interface);
				}
			}
		}

		for (int l_index : l_indices)
		{
			for (BigNumber iy = 0; iy < ny_l; iy++)
			{
				for (BigNumber ix = 0; ix < nx_l; ix++)
				{
					auto* element = dynamic_cast<RectangularPolygonalElement*>(this->Elements[indexE(ix, iy, l_index)]);

					if (ix != nx_l - 1)
					{
						// East
						auto* eastNeighbour = dynamic_cast<RectangularPolygonalElement*>(this->Elements[indexE(ix + 1, iy, l_index)]);
						auto* interface = new CartesianEdge(numberInterface++, eastNeighbour->BottomLeftCorner, eastNeighbour->TopLeftCorner, element, eastNeighbour, CartesianShapeOrientation::Vertical);
						this->Faces.push_back(interface);
						this->InteriorFaces.push_back(interface);
						element->AddEastFace(interface);
						eastNeighbour->AddWestFace(interface);
					}

					if (iy != ny_l - 1)
					{
						// North
						auto* northNeighbour = dynamic_cast<RectangularPolygonalElement*>(this->Elements[indexE(ix, iy + 1, l_index)]);
						auto* interface = new CartesianEdge(numberInterface++, northNeighbour->BottomLeftCorner, northNeighbour->BottomRightCorner, element, northNeighbour, CartesianShapeOrientation::Horizontal);
						this->Faces.push_back(interface);
						this->InteriorFaces.push_back(interface);
						element->AddNorthFace(interface);
						northNeighbour->AddSouthFace(interface);
					}
				}
			}
		}
	}

private:
	static bool isHighResolutionStripe(const int index)
	{
		return index % 2 == 0;
	}

	inline BigNumber indexV(BigNumber x, BigNumber y, int stripe)
	{
		if (stripe % 2 == 0)
		{
			// r-stripe
			//BigNumber offset = stripe == 0 ? 0 : stripe_vertex_breaks[stripe - 1];
			BigNumber offset = stripe_vertex_breaks[stripe];
			auto xxxx = offset + y * (Nx_r + 1) + x;;
			return offset + y * (Nx_r + 1) + x;
		}

		// l-stripe

		if (x == 0)
		{
			// return rightmost vertex of western r-stripe
			return indexV(Nx_r, ky * y, stripe - 1);
		}

		BigNumber offset = stripe_vertex_breaks[stripe];

		if (stripe != n_stripes - 1)
		{
			if (x == Nx_l)
			{
				// return leftmost vertex of eastern r-stripe
				return indexV(0, ky * y, stripe + 1);
			}

			return offset + y * (Nx_l - 1) + x - 1;
		}

		return offset + y * Nx_l + x - 1;
	}

	inline BigNumber indexE(BigNumber x, BigNumber y, int stripe)
	{
		//BigNumber offset = stripe == 0 ? 0 : stripe_element_breaks[stripe - 1];
		BigNumber offset = stripe_element_breaks[stripe];

		if (stripe % 2 == 0)
		{
			// r-stripe
			return offset + y * Nx_r + x;
		}

		// l-stripe
		return offset + y * Nx_l + x;
	}

	inline bool isFinalStripe(const int index) const
	{
		return index == n_stripes - 1;
	}

public:
	string Description() override
	{
		return "CartesianEmilMesh: " + to_string(this->Nx_l + Nx_r) + " x " + to_string(this->Ny_l + Ny_r);
	}

	string FileNamePart() override
	{
		string geo = this->With4Quadrants ? "square4quadrants" : "square";
		return geo + "-inhouse-cartpoly-n" + to_string(this->Nx_l + Nx_r);
	}

	string GeometryDescription() override
	{
		return this->With4Quadrants ? "Square with 4 quadrants" : "Square";
	}

	double H() override
	{
		return (1.0 - this->r_length) / this->Nx_l;
	}

	double Regularity() override
	{
		return min((double)this->Nx_l, (double)this->Ny_l) / max((double)this->Nx_l, (double)this->Ny_l);
	}

	double AverageH() override
	{
		return H();
	}

	void CoarsenMesh(H_CoarsStgy elemCoarseningStgy, FaceCoarseningStrategy faceCoarseningStgy, FaceCollapsing bdryFaceCollapsing, double coarseningFactor) override
	{
		if (elemCoarseningStgy == H_CoarsStgy::StandardCoarsening)
		{
			if (faceCoarseningStgy == FaceCoarseningStrategy::InterfaceCollapsing)
				StandardCoarsening();
			else
				Utils::FatalError("Unmanaged face coarsening strategy");
		}
		else
			//PolyhedralMesh<2>::CoarsenMesh(elemCoarseningStgy, faceCoarseningStgy, bdryFaceCollapsing, coarseningFactor);
			Utils::FatalError("Unmanaged face coarsening strategy");
	}

	void StandardCoarsening()
	{
		BigNumber nx_l = this->Nx_l;
		BigNumber ny_l = this->Ny_l;
		BigNumber nx_r = this->Nx_r;
		BigNumber ny_r = this->Ny_r;

		if (nx_l % 2 != 0 || ny_l % 2 != 0)
		{
			cout << "Error: impossible to build coarse mesh. Nx and Ny must be even: Nx = " << nx_l << ", Ny = " << ny_l << "." << endl;
			return;
		}

		auto* coarseMesh = new Square_CartesianPolyStripeMesh(nx_l / 2, ny_l / 2, this->n_stripes, this->With4Quadrants, false);
		this->InitializeCoarsening(coarseMesh);
		coarseMesh->ComesFrom.CS = H_CoarsStgy::StandardCoarsening;
		//coarseMesh->ComesFrom.nFineElementsByCoarseElement = 4;
		//coarseMesh->ComesFrom.nFineFacesAddedByCoarseElement = 4;
		//coarseMesh->ComesFrom.nFineFacesByKeptCoarseFace = 2;
		coarseMesh->Build();

		for (int r_index : r_indices)
		{
			for (BigNumber iy = 0; iy < ny_r; ++iy)
			{
				for (BigNumber ix = 0; ix < nx_r; ++ix)
				{
					auto* fineElement = dynamic_cast<RectangularPolygonalElement*>(this->Elements[indexE(ix, iy, r_index)]);
					auto* coarseElement = dynamic_cast<RectangularPolygonalElement*>(coarseMesh->Elements[coarseMesh->indexE(ix / 2, iy / 2, r_index)]);

					coarseElement->FinerElements.push_back(fineElement);
					fineElement->CoarserElement = coarseElement;

					if (iy % 2 == 0)
					{
						fineElement->NorthFaces[0]->IsRemovedOnCoarserGrid = true;
						coarseElement->FinerFacesRemoved.push_back(fineElement->NorthFaces[0]);
						coarseElement->SouthFaces[0]->FinerFaces.push_back(fineElement->SouthFaces[0]);
						fineElement->SouthFaces[0]->CoarseFace = coarseElement->SouthFaces[0];
					}
					else if (iy == ny_r - 1)
					{
						coarseElement->NorthFaces[0]->FinerFaces.push_back(fineElement->NorthFaces[0]);
						fineElement->NorthFaces[0]->CoarseFace = coarseElement->NorthFaces[0];
					}

					if (ix % 2 == 0)
					{
						fineElement->EastFaces[0]->IsRemovedOnCoarserGrid = true;
						coarseElement->FinerFacesRemoved.push_back(fineElement->EastFaces[0]);
						coarseElement->WestFaces[0]->FinerFaces.push_back(fineElement->WestFaces[0]);
						fineElement->WestFaces[0]->CoarseFace = coarseElement->WestFaces[0];
					}
					else if (ix == nx_r - 1)
					{
						coarseElement->EastFaces[0]->FinerFaces.push_back(fineElement->EastFaces[0]);
						fineElement->EastFaces[0]->CoarseFace = coarseElement->EastFaces[0];
					}
				}
			}
		}

		for (int l_index : l_indices)
		{
			for (BigNumber iy = 0; iy < ny_l; ++iy)
			{
				for (BigNumber ix = 0; ix < nx_l; ++ix)
				{
					auto* fineElement = dynamic_cast<RectangularPolygonalElement*>(this->Elements[indexE(ix, iy, l_index)]);
					auto* coarseElement = dynamic_cast<RectangularPolygonalElement*>(coarseMesh->Elements[coarseMesh->indexE(ix / 2, iy / 2, l_index)]);

					coarseElement->FinerElements.push_back(fineElement);
					fineElement->CoarserElement = coarseElement;

					if (iy % 2 == 0)
					{
						fineElement->NorthFaces[0]->IsRemovedOnCoarserGrid = true;
						coarseElement->FinerFacesRemoved.push_back(fineElement->NorthFaces[0]);
						coarseElement->SouthFaces[0]->FinerFaces.push_back(fineElement->SouthFaces[0]);
						fineElement->SouthFaces[0]->CoarseFace = coarseElement->SouthFaces[0];
					}
					else if (iy == ny_l - 1)
					{
						coarseElement->NorthFaces[0]->FinerFaces.push_back(fineElement->NorthFaces[0]);
						fineElement->NorthFaces[0]->CoarseFace = coarseElement->NorthFaces[0];
					}

					if (ix % 2 == 0 && ix > 0)
					{
						fineElement->EastFaces[0]->IsRemovedOnCoarserGrid = true;
						coarseElement->FinerFacesRemoved.push_back(fineElement->EastFaces[0]);
						coarseElement->WestFaces[0]->FinerFaces.push_back(fineElement->WestFaces[0]);
						fineElement->WestFaces[0]->CoarseFace = coarseElement->WestFaces[0];
					}
				}
			}

			if (isFinalStripe(l_index))
			{
				// account for east domain boundary
				BigNumber ix = nx_l - 1;

				for (BigNumber iy = 0; iy < ny_l; ++iy)
				{
					auto* fineElement = dynamic_cast<RectangularPolygonalElement*>(this->Elements[indexE(ix, iy, l_index)]);
					auto* coarseElement = dynamic_cast<RectangularPolygonalElement*>(coarseMesh->Elements[coarseMesh->indexE(ix / 2, iy / 2, l_index)]);

					coarseElement->EastFaces[0]->FinerFaces.push_back(fineElement->EastFaces[0]);
					fineElement->EastFaces[0]->CoarseFace = coarseElement->EastFaces[0];
				}
			}
		}

		this->FinalizeCoarsening();
	}

	Mesh<2>* Copy() override
	{
		return new Square_CartesianPolyStripeMesh(this->Nx_l, this->Ny_l, this->With4Quadrants);
	}
};
