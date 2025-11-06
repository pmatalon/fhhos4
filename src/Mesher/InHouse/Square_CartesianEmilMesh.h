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

class Square_CartesianEmilMesh : public PolyhedralMesh<2>
{
public:
	BigNumber Nx_l; //left half stat
	BigNumber Ny_l; //left half stat
	int k = 16; // factor relating right half stats

	BigNumber Nx_r; //right half stat
	BigNumber Ny_r; //right half stat
	bool With4Quadrants;


	Square_CartesianEmilMesh(BigNumber nx, BigNumber ny, bool with4Quadrants = false, bool buildMesh = true) : PolyhedralMesh()
	{
		// nx = ny falls down to square elements
		this->Nx_l = nx;
		this->Ny_l = ny;
		this->Nx_r = k*nx;
		this->Ny_r = k*ny;
		this->With4Quadrants = with4Quadrants;

		if (with4Quadrants && (nx % 2 == 1 || ny % 2 == 1))
			Utils::FatalError("Building the mesh for a square with 4 quadrants requires the number of subdivisions in each direction to be even.");

		if (buildMesh)
			Build();
	}

	void Build() {
		BigNumber nx_l = this->Nx_l;
		BigNumber ny_l = this->Ny_l;

		BigNumber nx_r = this->Nx_r;
		BigNumber ny_r = this->Ny_r;

		double hx_l = 0.5 / nx_l;
		double hy_l = 0.5 / ny_l;
		double hx_r = 0.5 / nx_r;
		double hy_r = 0.5 / ny_r;

		// Physical parts
		PhysicalGroup<2>* domain = nullptr;
		PhysicalGroup<2>* quadrantBottomLeft = nullptr;
		PhysicalGroup<2>* quadrantBottomRight = nullptr;
		PhysicalGroup<2>* quadrantTopRight = nullptr;
		PhysicalGroup<2>* quadrantTopLeft = nullptr;
		if (this->With4Quadrants)
		{
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

		this->Vertices.reserve(nx_l * (ny_l + 1) + (nx_r + 1) * (ny_r + 1));
		// left half vertices
		for (BigNumber iy = 0; iy < ny_l + 1; ++iy)
		{
			for (BigNumber ix = 0; ix < nx_l; ++ix)
			{
				auto* vertex = new MeshVertex<2>(indexV_l(ix, iy), ix * hx_l, iy * hy_l);
				this->Vertices.push_back(vertex);
			}
		}
		// right half vertices
		for (BigNumber iy = 0; iy < ny_r + 1; ++iy)
		{
			for (BigNumber ix = 0; ix < nx_r + 1; ++ix)
			{
				auto* vertex = new MeshVertex<2>(indexV_r(ix, iy), 0.5 + ix * hx_r, iy * hy_r);
				this->Vertices.push_back(vertex);
			}
		}

		//----------//
		// Elements //
		//----------//

		this->Elements.reserve(nx_l * ny_l + nx_r * ny_r);
		// left half elements
		for (BigNumber iy = 0; iy < ny_l; ++iy)
		{
			for (BigNumber ix = 0; ix < nx_l - 1; ++ix)
			{
				BigNumber number = indexE_l(ix, iy);
				auto* bottomLeftCorner  = static_cast<MeshVertex<2>*>(Vertices[indexV_l(ix,     iy    )]);
				auto* topLeftCorner     = static_cast<MeshVertex<2>*>(Vertices[indexV_l(ix,     iy + 1)]);
				auto* topRightCorner    = static_cast<MeshVertex<2>*>(Vertices[indexV_l(ix + 1, iy + 1)]);
				auto* bottomRightCorner = static_cast<MeshVertex<2>*>(Vertices[indexV_l(ix + 1, iy    )]);
				auto* rectangle = new RectangularElement(number, bottomLeftCorner, topLeftCorner, topRightCorner, bottomRightCorner);
				this->Elements.push_back(rectangle);

				bottomLeftCorner->Elements.push_back(rectangle);
				topLeftCorner->Elements.push_back(rectangle);
				topRightCorner->Elements.push_back(rectangle);
				bottomRightCorner->Elements.push_back(rectangle);

				rectangle->PhysicalPart = domain;
			}

			// elements with boundary on midline
			{
				BigNumber ix = nx_l - 1;
				BigNumber number = indexE_l(ix, iy);
				Vertex* bottomLeftCorner  = Vertices[indexV_l(ix,     iy    )];
				Vertex* topLeftCorner     = Vertices[indexV_l(ix,     iy + 1)];
				Vertex* topRightCorner    = Vertices[indexV_r(0, k * (iy + 1))];
				Vertex* bottomRightCorner = Vertices[indexV_r(0, k * iy    )];
				auto* rectangle = new RectangularPolygonalElement(number, bottomLeftCorner, topLeftCorner, topRightCorner, bottomRightCorner);
				this->Elements.push_back(rectangle);

				rectangle->PhysicalPart = domain;
			}
		}

		// right half elements
		for (BigNumber iy = 0; iy < ny_r; ++iy)
		{
			for (BigNumber ix = 0; ix < nx_r; ++ix)
			{
				BigNumber number = indexE_r(ix, iy);
				auto* bottomLeftCorner  = static_cast<MeshVertex<2>*>(Vertices[indexV_r(ix,     iy    )]);
				auto* topLeftCorner     = static_cast<MeshVertex<2>*>(Vertices[indexV_r(ix,     iy + 1)]);
				auto* topRightCorner    = static_cast<MeshVertex<2>*>(Vertices[indexV_r(ix + 1, iy + 1)]);
				auto* bottomRightCorner = static_cast<MeshVertex<2>*>(Vertices[indexV_r(ix + 1, iy    )]);
				auto* rectangle = new RectangularElement(number, bottomLeftCorner, topLeftCorner, topRightCorner, bottomRightCorner);
				this->Elements.push_back(rectangle);

				bottomLeftCorner->Elements.push_back(rectangle);
				topLeftCorner->Elements.push_back(rectangle);
				topRightCorner->Elements.push_back(rectangle);
				bottomRightCorner->Elements.push_back(rectangle);

				rectangle->PhysicalPart = domain;
			}
		}

		//-------//
		// Faces //
		//-------//

		this->Faces.reserve(nx_l * (ny_l + 1) + ny_l * nx_l +
							  nx_r * (ny_r + 1) + ny_r * (nx_r + 1));
		BigNumber numberInterface = 0;

		// left half horizontal boundary
		for (BigNumber ix = 0; ix < nx_l - 1; ++ix)
		{
			// South boundary
			auto* rectangle = dynamic_cast<RectangularElement*>(this->Elements[indexE_l(ix, 0)]);
			CartesianEdge* southBoundary = new CartesianEdge(numberInterface++, rectangle->BottomLeftCorner, rectangle->BottomRightCorner, rectangle, CartesianShapeOrientation::Horizontal);
			this->Faces.push_back(southBoundary);
			this->BoundaryFaces.push_back(southBoundary);
			rectangle->SetSouthInterface(southBoundary);
			static_cast<MeshVertex<2>*>(rectangle->BottomLeftCorner)->Faces.push_back(southBoundary);
			static_cast<MeshVertex<2>*>(rectangle->BottomRightCorner)->Faces.push_back(southBoundary);
			southBoundary->BoundaryPart = squareBottomBoundary;

			// North boundary
			rectangle = dynamic_cast<RectangularElement*>(this->Elements[indexE_l(ix, ny_l - 1)]);
			CartesianEdge* northBoundary = new CartesianEdge(numberInterface++, rectangle->TopLeftCorner, rectangle->TopRightCorner, rectangle, CartesianShapeOrientation::Horizontal);
			this->Faces.push_back(northBoundary);
			this->BoundaryFaces.push_back(northBoundary);
			rectangle->SetNorthInterface(northBoundary);
			static_cast<MeshVertex<2>*>(rectangle->TopLeftCorner)->Faces.push_back(northBoundary);
			static_cast<MeshVertex<2>*>(rectangle->TopRightCorner)->Faces.push_back(northBoundary);
			northBoundary->BoundaryPart = squareTopBoundary;
		}

		{
			// South boundary
			RectangularPolygonalElement* rectangle = dynamic_cast<RectangularPolygonalElement*>(this->Elements[indexE_l(nx_l - 1, 0)]);
			CartesianEdge* southBoundary = new CartesianEdge(numberInterface++, rectangle->BottomLeftCorner, rectangle->BottomRightCorner, rectangle, CartesianShapeOrientation::Horizontal);
			this->Faces.push_back(southBoundary);
			this->BoundaryFaces.push_back(southBoundary);
			rectangle->AddSouthFace(southBoundary);
			southBoundary->BoundaryPart = squareBottomBoundary;

			// North boundary
			rectangle = dynamic_cast<RectangularPolygonalElement*>(this->Elements[indexE_l(nx_l - 1, ny_l - 1)]);
			CartesianEdge* northBoundary = new CartesianEdge(numberInterface++, rectangle->TopLeftCorner, rectangle->TopRightCorner, rectangle, CartesianShapeOrientation::Horizontal);
			this->Faces.push_back(northBoundary);
			this->BoundaryFaces.push_back(northBoundary);
			rectangle->AddNorthFace(northBoundary);
			northBoundary->BoundaryPart = squareTopBoundary;
		}

		// right half horizontal boundary
		for (BigNumber ix = 0; ix < nx_r; ++ix)
		{
			// South boundary
			auto* rectangle = dynamic_cast<RectangularElement*>(this->Elements[indexE_r(ix, 0)]);
			CartesianEdge* southBoundary = new CartesianEdge(numberInterface++, rectangle->BottomLeftCorner, rectangle->BottomRightCorner, rectangle, CartesianShapeOrientation::Horizontal);
			this->Faces.push_back(southBoundary);
			this->BoundaryFaces.push_back(southBoundary);
			rectangle->SetSouthInterface(southBoundary);
			static_cast<MeshVertex<2>*>(rectangle->BottomLeftCorner)->Faces.push_back(southBoundary);
			static_cast<MeshVertex<2>*>(rectangle->BottomRightCorner)->Faces.push_back(southBoundary);
			southBoundary->BoundaryPart = squareBottomBoundary;

			// North boundary
			rectangle = dynamic_cast<RectangularElement*>(this->Elements[indexE_r(ix, ny_r - 1)]);
			CartesianEdge* northBoundary = new CartesianEdge(numberInterface++, rectangle->TopLeftCorner, rectangle->TopRightCorner, rectangle, CartesianShapeOrientation::Horizontal);
			this->Faces.push_back(northBoundary);
			this->BoundaryFaces.push_back(northBoundary);
			rectangle->SetNorthInterface(northBoundary);
			static_cast<MeshVertex<2>*>(rectangle->TopLeftCorner)->Faces.push_back(northBoundary);
			static_cast<MeshVertex<2>*>(rectangle->TopRightCorner)->Faces.push_back(northBoundary);
			northBoundary->BoundaryPart = squareTopBoundary;
		}

		// left half vertical boundary
		for (BigNumber iy = 0; iy < ny_l; ++iy)
		{
			// West boundary
			RectangularElement* rectangle = dynamic_cast<RectangularElement*>(this->Elements[indexE_l(0, iy)]);
			CartesianEdge* westBoundary = new CartesianEdge(numberInterface++, rectangle->BottomLeftCorner, rectangle->TopLeftCorner, rectangle, CartesianShapeOrientation::Vertical);
			this->Faces.push_back(westBoundary);
			this->BoundaryFaces.push_back(westBoundary);
			rectangle->SetWestInterface(westBoundary);
			static_cast<MeshVertex<2>*>(rectangle->BottomLeftCorner)->Faces.push_back(westBoundary);
			static_cast<MeshVertex<2>*>(rectangle->TopLeftCorner)->Faces.push_back(westBoundary);
			westBoundary->BoundaryPart = squareLeftBoundary;
		}

		// right half vertical boundary
		for (BigNumber iy = 0; iy < ny_r; ++iy)
		{
			// East boundary
			auto* rectangle = dynamic_cast<RectangularElement*>(this->Elements[indexE_r(nx_r - 1, iy)]);
			auto* eastBoundary = new CartesianEdge(numberInterface++, rectangle->BottomRightCorner, rectangle->TopRightCorner, rectangle, CartesianShapeOrientation::Vertical);
			this->Faces.push_back(eastBoundary);
			this->BoundaryFaces.push_back(eastBoundary);
			rectangle->SetEastInterface(eastBoundary);
			static_cast<MeshVertex<2>*>(rectangle->BottomRightCorner)->Faces.push_back(eastBoundary);
			static_cast<MeshVertex<2>*>(rectangle->TopRightCorner)->Faces.push_back(eastBoundary);
			eastBoundary->BoundaryPart = squareRightBoundary;
		}

		// left side internal faces
		for (BigNumber iy = 0; iy < ny_l; iy++)
		{
			for (BigNumber ix = 0; ix < nx_l - 1; ix++)
			{
				RectangularElement* element = dynamic_cast<RectangularElement*>(this->Elements[indexE_l(ix, iy)]);

				// East
				if (ix != nx_l - 2)
				{
					RectangularElement* eastNeighbour = dynamic_cast<RectangularElement*>(this->Elements[indexE_l(ix + 1, iy)]);
					CartesianEdge* interface = new CartesianEdge(numberInterface++, eastNeighbour->BottomLeftCorner, eastNeighbour->TopLeftCorner, element, eastNeighbour, CartesianShapeOrientation::Vertical);
					this->Faces.push_back(interface);
					this->InteriorFaces.push_back(interface);
					element->SetEastInterface(interface);
					eastNeighbour->SetWestInterface(interface);
					static_cast<MeshVertex<2>*>(eastNeighbour->BottomLeftCorner)->Faces.push_back(interface);
					static_cast<MeshVertex<2>*>(eastNeighbour->TopLeftCorner)->Faces.push_back(interface);
				}
				else
				{
					RectangularPolygonalElement* eastNeighbour = dynamic_cast<RectangularPolygonalElement*>(this->Elements[indexE_l(ix + 1, iy)]);
					CartesianEdge* interface = new CartesianEdge(numberInterface++, eastNeighbour->BottomLeftCorner, eastNeighbour->TopLeftCorner, element, eastNeighbour, CartesianShapeOrientation::Vertical);
					this->Faces.push_back(interface);
					this->InteriorFaces.push_back(interface);
					element->SetEastInterface(interface);
					eastNeighbour->AddWestFace(interface);
					// why not link face to vertices here?
				}

				// North
				if (iy != ny_l - 1)
				{
					RectangularElement* northNeighbour = dynamic_cast<RectangularElement*>(this->Elements[indexE_l(ix, iy + 1)]);
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

		// midline vertical faces
		for (BigNumber iy = 0; iy < ny_r; iy++)
		{
			RectangularElement* element = dynamic_cast<RectangularElement*>(this->Elements[indexE_r(0, iy)]);
			RectangularPolygonalElement* westNeighbour = dynamic_cast<RectangularPolygonalElement*>(this->Elements[indexE_l(nx_l - 1, iy / k)]);
			CartesianEdge* interface = new CartesianEdge(numberInterface++, element->BottomLeftCorner, element->TopLeftCorner, westNeighbour, element, CartesianShapeOrientation::Vertical);
			this->Faces.push_back(interface);
			this->InteriorFaces.push_back(interface);
			element->SetWestInterface(interface);
			westNeighbour->AddEastFace(interface);
			static_cast<MeshVertex<2>*>(element->BottomLeftCorner)->Faces.push_back(interface);
			static_cast<MeshVertex<2>*>(element->TopLeftCorner)->Faces.push_back(interface);
		}

		// midline horizontal faces
		for (BigNumber iy = 0; iy < ny_l - 1; iy++)
		{
			RectangularPolygonalElement* element = dynamic_cast<RectangularPolygonalElement*>(this->Elements[indexE_l(nx_l - 1, iy)]);
			RectangularPolygonalElement* northNeighbour = dynamic_cast<RectangularPolygonalElement*>(this->Elements[indexE_l(nx_l - 1, iy + 1)]);
			CartesianEdge* interface = new CartesianEdge(numberInterface++, northNeighbour->BottomLeftCorner, northNeighbour->BottomRightCorner, element, northNeighbour, CartesianShapeOrientation::Horizontal);
			this->Faces.push_back(interface);
			this->InteriorFaces.push_back(interface);
			element->AddNorthFace(interface);
			northNeighbour->AddSouthFace(interface);
		}

		// right side internal faces
		for (BigNumber iy = 0; iy < ny_r; iy++)
		{
			for (BigNumber ix = 0; ix < nx_r; ix++)
			{
				RectangularElement* element = dynamic_cast<RectangularElement*>(this->Elements[indexE_r(ix, iy)]);
				if (ix != nx_r - 1)
				{
					// East
					RectangularElement* eastNeighbour = dynamic_cast<RectangularElement*>(this->Elements[indexE_r(ix + 1, iy)]);
					CartesianEdge* interface = new CartesianEdge(numberInterface++, eastNeighbour->BottomLeftCorner, eastNeighbour->TopLeftCorner, element, eastNeighbour, CartesianShapeOrientation::Vertical);
					this->Faces.push_back(interface);
					this->InteriorFaces.push_back(interface);
					element->SetEastInterface(interface);
					eastNeighbour->SetWestInterface(interface);
					static_cast<MeshVertex<2>*>(eastNeighbour->BottomLeftCorner)->Faces.push_back(interface);
					static_cast<MeshVertex<2>*>(eastNeighbour->TopLeftCorner)->Faces.push_back(interface);
				}
				if (iy != ny_r - 1)
				{
					// North
					RectangularElement* northNeighbour = dynamic_cast<RectangularElement*>(this->Elements[indexE_r(ix, iy + 1)]);
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
	//inline BigNumber indexV(BigNumber x, BigNumber y)
	//{
//		return y * (Nx + 1) + x;
	//}

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
		return 0.5 / this->Nx_l;
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
			if (faceCoarseningStgy == FaceCoarseningStrategy::None)
				//CoarsenByAgglomerationAndKeepFineFaces();
				Utils::FatalError("Not imp yet");
			else
				Utils::FatalError("Unmanaged face coarsening strategy");
		}
		else
			PolyhedralMesh<2>::CoarsenMesh(elemCoarseningStgy, faceCoarseningStgy, bdryFaceCollapsing, coarseningFactor);
	}


	/*
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
		coarseMesh->ComesFrom.CS = H_CoarsStgy::StandardCoarsening;
		//coarseMesh->ComesFrom.nFineElementsByCoarseElement = 4;
		//coarseMesh->ComesFrom.nFineFacesAddedByCoarseElement = 4;
		//coarseMesh->ComesFrom.nFineFacesByKeptCoarseFace = 1;

		BigNumber faceNumber = 0;
		for (BigNumber i = 0; i < ny / 2; ++i)
		{
			for (BigNumber j = 0; j < nx / 2; ++j)
			{
				RectangularPolygonalElement* bottomLeftElement = dynamic_cast<RectangularPolygonalElement*>(this->Elements[2 * i * nx + 2 * j]);
				RectangularPolygonalElement* bottomRightElement = dynamic_cast<RectangularPolygonalElement*>(this->Elements[2 * i * nx + 2 * j + 1]);
				RectangularPolygonalElement* topLeftElement = dynamic_cast<RectangularPolygonalElement*>(this->Elements[(2 * i + 1) * nx + 2 * j]);
				RectangularPolygonalElement* topRightElement = dynamic_cast<RectangularPolygonalElement*>(this->Elements[(2 * i + 1) * nx + 2 * j + 1]);

				// Coarse element
				RectangularPolygonalElement* coarseElement = new RectangularPolygonalElement(i*nx / 2 + j, bottomLeftElement->BottomLeftCorner, topLeftElement->TopLeftCorner, topRightElement->TopRightCorner, bottomRightElement->BottomRightCorner);
				coarseMesh->Elements.push_back(coarseElement);

				coarseElement->PhysicalPart = bottomLeftElement->PhysicalPart;

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
						f->CoarseFace = newFace;
						newFace->FinerFaces.push_back(f);
						coarseElement->AddWestFace(newFace);
						coarseMesh->AddFace(newFace);
					}
					for (Face<2>* f : topLeftElement->WestFaces)
					{
						Face<2>* newFace = f->CreateSameGeometricFace(faceNumber++, coarseElement);
						f->CoarseFace = newFace;
						newFace->FinerFaces.push_back(f);
						coarseElement->AddWestFace(newFace);
						coarseMesh->AddFace(newFace);
					}
				}
				else
				{
					RectangularPolygonalElement* leftNeighbour = dynamic_cast<RectangularPolygonalElement*>(coarseMesh->Elements[i * nx / 2 + j - 1]);
					coarseElement->SetWestFacesFromNeighbour(leftNeighbour);
				}

				// North faces
				for (Face<2>* f : topLeftElement->NorthFaces)
				{
					Face<2>* newFace = f->CreateSameGeometricFace(faceNumber++, coarseElement);
					f->CoarseFace = newFace;
					newFace->FinerFaces.push_back(f);
					coarseElement->AddNorthFace(newFace);
					coarseMesh->AddFace(newFace);
				}
				for (Face<2>* f : topRightElement->NorthFaces)
				{
					Face<2>* newFace = f->CreateSameGeometricFace(faceNumber++, coarseElement);
					f->CoarseFace = newFace;
					newFace->FinerFaces.push_back(f);
					coarseElement->AddNorthFace(newFace);
					coarseMesh->AddFace(newFace);
				}

				// East faces
				for (Face<2>* f : topRightElement->EastFaces)
				{
					Face<2>* newFace = f->CreateSameGeometricFace(faceNumber++, coarseElement);
					f->CoarseFace = newFace;
					newFace->FinerFaces.push_back(f);
					coarseElement->AddEastFace(newFace);
					coarseMesh->AddFace(newFace);
				}
				for (Face<2>* f : bottomRightElement->EastFaces)
				{
					Face<2>* newFace = f->CreateSameGeometricFace(faceNumber++, coarseElement);
					f->CoarseFace = newFace;
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
						f->CoarseFace = newFace;
						newFace->FinerFaces.push_back(f);
						coarseElement->AddSouthFace(newFace);
						coarseMesh->AddFace(newFace);
					}
					for (Face<2>* f : bottomLeftElement->SouthFaces)
					{
						Face<2>* newFace = f->CreateSameGeometricFace(faceNumber++, coarseElement);
						f->CoarseFace = newFace;
						newFace->FinerFaces.push_back(f);
						coarseElement->AddSouthFace(newFace);
						coarseMesh->AddFace(newFace);
					}
				}
				else
				{
					RectangularPolygonalElement* bottomNeighbour = dynamic_cast<RectangularPolygonalElement*>(coarseMesh->Elements[(i - 1) * nx / 2 + j]);
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
		this->FinalizeCoarsening();
	}
	*/

	Mesh<2>* Copy() override
	{
		return new Square_CartesianEmilMesh(this->Nx_l, this->Ny_l, this->With4Quadrants);
	}
};
