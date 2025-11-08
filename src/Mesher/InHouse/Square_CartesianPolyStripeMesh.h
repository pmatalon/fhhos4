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
	BigNumber Nx_l; //left half stat
	BigNumber Ny_l; //left half stat

	BigNumber Nx_r; //right half stat
	BigNumber Ny_r; //right half stat
	bool With4Quadrants;

	//double midline = 0.9;
	int ky = 16; // factor relating right half stats
	//vector<double> breaks = {0.1, 0.9};
	int n_breaks = 1;
	double r_length = 0.1;
	vector<double> breaks = {0.1};
	int n_stripes = breaks.size() + 1;
	//vector<int> r_indices = {0, 2};
	vector<int> r_indices = {0};
	vector<int> l_indices = {1};

	vector<BigNumber> stripe_vertex_breaks = {0};
	vector<BigNumber> stripe_element_breaks = {0};

	Square_CartesianPolyStripeMesh(BigNumber nx, BigNumber ny, bool with4Quadrants = false, bool buildMesh = true) : PolyhedralMesh()
	{
		// nx = ny falls down to square elements
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

		// todo
		double hx_l = (1 - r_length) / nx_l;
		double hy_l = 1.0 / ny_l;
		double hx_r = r_length / nx_r;
		double hy_r = 1.0 / ny_r;

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

		this->Vertices.reserve(stripe_vertex_breaks.back());

		cout << "vertices before adding them: " << Vertices.size() << endl;

		for (int r_index : r_indices)
		{
			for (BigNumber iy = 0; iy < ny_r + 1; ++iy)
			{
				for (BigNumber ix = 0; ix < nx_r + 1; ++ix)
				{
					auto* vertex = new Vertex(indexV(ix, iy, r_index), 0.0 + ix * hx_r, iy * hy_r);
					this->Vertices.push_back(vertex);
				}
			}
		}

		for (int l_index : l_indices)
		{
			for (BigNumber iy = 0; iy < ny_l + 1; ++iy)
			{
				for (BigNumber ix = 1; ix < nx_l; ++ix)
				{
					auto* vertex = new Vertex(indexV(ix, iy, l_index), r_length + ix * hx_l, iy * hy_l);
					this->Vertices.push_back(vertex);
				}

				if (l_index == n_stripes - 1)
				{
					// include vertices on left boundary
					BigNumber ix = nx_l;
					auto* vertex = new Vertex(indexV(ix, iy, l_index), r_length + ix * hx_l, iy * hy_l);
					this->Vertices.push_back(vertex);
				}
			}
		}

		if (Vertices.size()!=stripe_vertex_breaks.back())
			cout << "Vertices size vs what should be: " << Vertices.size() << "/" << stripe_vertex_breaks.back() << endl;


		//----------//
		// Elements //
		//----------//

		//this->Elements.reserve( r_indices.size() * nx_r * ny_r + l_indices.size() * nx_l * ny_l);
		this->Elements.reserve( stripe_element_breaks.back());

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
					auto* rectangle = new RectangularPolygonalElement(indexE(ix, iy, r_index), bottomLeftCorner, topLeftCorner, topRightCorner, bottomRightCorner);
					this->Elements.push_back(rectangle);
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
					auto* rectangle = new RectangularPolygonalElement(indexE(ix, iy, l_index), bottomLeftCorner, topLeftCorner, topRightCorner, bottomRightCorner);
					this->Elements.push_back(rectangle);
					rectangle->PhysicalPart = domain;
				}
			}
		}

		if (Elements.size()!=stripe_element_breaks.back())
			cout << "Elements size vs what should be: " << Elements.size() << "/" << stripe_element_breaks.size() << endl;

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

	/*
	inline BigNumber fetch_rIndexE(BigNumber x, BigNumber y, int stripe=0)
	{
		if (x == -1)
			return indexE(Nx_l - 1, y / ky, stripe - 1);
		else if (x == Nx_r)
			return indexE(0, y / ky, stripe + 1);

		return indexE(x, y, stripe);
	}

	inline BigNumber indexV_l(BigNumber x, BigNumber y)
	{
		//return y * (Nx_l + Nx_r + 1) + x;
		//return y * (Nx_l + 1) + x;
		return y * Nx_l + x;
	}

	inline BigNumber indexV_r(BigNumber x, BigNumber y)
	{
		//if (x == 0 && ((y % k) == 0))
		//{
	    //		// vertex coincides with
		//}
		//return y * (Nx_l + Nx_r + 1) + Nx_l + x;
		//return (Nx_l + 1) * (Ny_l + 1) + y * (Nx_r + 1) + x;
		return Nx_l * (Ny_l + 1) + y * (Nx_r + 1) + x;
	}

	inline BigNumber indexE_l(BigNumber x, BigNumber y)
	{
		return y * Nx_l + x;
	}

	inline BigNumber indexE_r(BigNumber x, BigNumber y)
	{
		return Nx_l * Ny_l + y * Nx_r + x;
	}
	*/

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
				//StandardCoarsening();
				Utils::FatalError("Unmanaged face coarsening strategy");
			else
				Utils::FatalError("Unmanaged face coarsening strategy");
		}
		else
			PolyhedralMesh<2>::CoarsenMesh(elemCoarseningStgy, faceCoarseningStgy, bdryFaceCollapsing, coarseningFactor);
	}
	/*
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

		auto* coarseMesh = new Square_CartesianPolyStripeMesh(nx_l / 2, ny_l / 2, this->With4Quadrants, false);
		this->InitializeCoarsening(coarseMesh);
		coarseMesh->ComesFrom.CS = H_CoarsStgy::StandardCoarsening;
		//coarseMesh->ComesFrom.nFineElementsByCoarseElement = 4;
		//coarseMesh->ComesFrom.nFineFacesAddedByCoarseElement = 4;
		//coarseMesh->ComesFrom.nFineFacesByKeptCoarseFace = 2;
		coarseMesh->Build();

		// left side excluding midline
		for (BigNumber iy = 0; iy < ny_l; ++iy)
		{
			for (BigNumber ix = 0; ix < nx_l - 2; ++ix)
			{
				RectangularElement* fineElement = dynamic_cast<RectangularElement*>(this->Elements[indexE_l(ix, iy)]);
				RectangularElement* coarseElement = dynamic_cast<RectangularElement*>(coarseMesh->Elements[coarseMesh->indexE_l(ix / 2, iy / 2)]);

				coarseElement->FinerElements.push_back(fineElement);
				fineElement->CoarserElement = coarseElement;
				if (iy % 2 == 0 && !fineElement->NorthFace->IsDomainBoundary)
				{
					fineElement->NorthFace->IsRemovedOnCoarserGrid = true;
					coarseElement->FinerFacesRemoved.push_back(fineElement->NorthFace);
					coarseElement->SouthFace->FinerFaces.push_back(fineElement->SouthFace);
					fineElement->SouthFace->CoarseFace = coarseElement->SouthFace;
				}
				if (iy == ny_l - 1)
				{
					coarseElement->NorthFace->FinerFaces.push_back(fineElement->NorthFace);
					fineElement->NorthFace->CoarseFace = coarseElement->NorthFace;
				}

				if (ix % 2 == 0 && !fineElement->EastFace->IsDomainBoundary)
				{
					fineElement->EastFace->IsRemovedOnCoarserGrid = true;
					coarseElement->FinerFacesRemoved.push_back(fineElement->EastFace);
					coarseElement->WestFace->FinerFaces.push_back(fineElement->WestFace);
					fineElement->WestFace->CoarseFace = coarseElement->WestFace;
				}
			}
		}

		// midline
		for (BigNumber iy = 0; iy < ny_l; ++iy)
		{
			auto* coarseElement = dynamic_cast<RectangularPolygonalElement*>(coarseMesh->Elements[coarseMesh->indexE_l((nx_l - 1) / 2, iy / 2)]);

			{
				auto* fineElement = dynamic_cast<RectangularElement*>(this->Elements[indexE_l(nx_l - 2, iy)]);

				coarseElement->FinerElements.push_back(fineElement);
				fineElement->CoarserElement = coarseElement;
				//if (iy % 2 == 0 && !fineElement->NorthFace->IsDomainBoundary)
				if (iy % 2 == 0)
				{
					fineElement->NorthFace->IsRemovedOnCoarserGrid = true;
					coarseElement->FinerFacesRemoved.push_back(fineElement->NorthFace);
					coarseElement->SouthFaces[0]->FinerFaces.push_back(fineElement->SouthFace);
					fineElement->SouthFace->CoarseFace = coarseElement->SouthFaces[0];
				}
				if (iy == ny_l - 1)
				{
					coarseElement->NorthFaces[0]->FinerFaces.push_back(fineElement->NorthFace);
					fineElement->NorthFace->CoarseFace = coarseElement->NorthFaces[0];
				}

				fineElement->EastFace->IsRemovedOnCoarserGrid = true;
				coarseElement->FinerFacesRemoved.push_back(fineElement->EastFace);
			}
			{
				auto* fineElement = dynamic_cast<RectangularPolygonalElement*>(this->Elements[indexE_l(nx_l - 1, iy)]);

				coarseElement->FinerElements.push_back(fineElement);
				fineElement->CoarserElement = coarseElement;
				//if (iy % 2 == 0 && !fineElement->NorthFace->IsDomainBoundary)
				if (iy % 2 == 0)
				{
					// there is only one south/north/west face
					fineElement->NorthFaces[0]->IsRemovedOnCoarserGrid = true;
					coarseElement->FinerFacesRemoved.push_back(fineElement->NorthFaces[0]);
					coarseElement->SouthFaces[0]->FinerFaces.push_back(fineElement->SouthFaces[0]);
					fineElement->SouthFaces[0]->CoarseFace = coarseElement->SouthFaces[0];
				}
				if (iy == ny_l - 1)
				{
					coarseElement->NorthFaces[0]->FinerFaces.push_back(fineElement->NorthFaces[0]);
					fineElement->NorthFaces[0]->CoarseFace = coarseElement->NorthFaces[0];
				}
			}
		}

		// right side
		for (BigNumber iy = 0; iy < ny_r; ++iy)
		{
			for (BigNumber ix = 0; ix < nx_r; ++ix)
			{
				RectangularElement* fineElement = dynamic_cast<RectangularElement*>(this->Elements[indexE_r(ix, iy)]);
				RectangularElement* coarseElement = dynamic_cast<RectangularElement*>(coarseMesh->Elements[coarseMesh->indexE_r(ix / 2, iy / 2)]);

				coarseElement->FinerElements.push_back(fineElement);
				fineElement->CoarserElement = coarseElement;
				if (iy % 2 == 0 && !fineElement->NorthFace->IsDomainBoundary)
				{
					fineElement->NorthFace->IsRemovedOnCoarserGrid = true;
					coarseElement->FinerFacesRemoved.push_back(fineElement->NorthFace);
					coarseElement->SouthFace->FinerFaces.push_back(fineElement->SouthFace);
					fineElement->SouthFace->CoarseFace = coarseElement->SouthFace;
				}
				if (iy == ny_l - 1)
				{
					coarseElement->NorthFace->FinerFaces.push_back(fineElement->NorthFace);
					fineElement->NorthFace->CoarseFace = coarseElement->NorthFace;
				}

				if (ix % 2 == 0 && !fineElement->EastFace->IsDomainBoundary)
				{
					fineElement->EastFace->IsRemovedOnCoarserGrid = true;
					coarseElement->FinerFacesRemoved.push_back(fineElement->EastFace);
					coarseElement->WestFace->FinerFaces.push_back(fineElement->WestFace);
					fineElement->WestFace->CoarseFace = coarseElement->WestFace;
				}
				if (ix == nx_l - 1)
				{
					coarseElement->EastFace->FinerFaces.push_back(fineElement->EastFace);
					fineElement->EastFace->CoarseFace = coarseElement->EastFace;
				}
			}
		}

		this->FinalizeCoarsening();
	}
	*/

	Mesh<2>* Copy() override
	{
		return new Square_CartesianPolyStripeMesh(this->Nx_l, this->Ny_l, this->With4Quadrants);
	}
};
