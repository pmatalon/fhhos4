#pragma once
#include "../../Mesh/3D/TetrahedralMesh.h"
#include "Cube_CartesianMesh.h"
using namespace std;

class Cube_CartesianTetrahedralMesh : public TetrahedralMesh
{
private:
	Cube_CartesianMesh* _cartMesh = nullptr;
	map<ParallelepipedElement*, vector<TetrahedralElement*>> _tetrasInCube;

public:
	Cube_CartesianTetrahedralMesh(BigNumber n) : Cube_CartesianTetrahedralMesh(new Cube_CartesianMesh(n, n, n))
	{}

	Cube_CartesianTetrahedralMesh(Cube_CartesianMesh* cartMesh, bool buildMesh = true) : TetrahedralMesh()
	{
		_cartMesh = cartMesh;
		if (buildMesh)
			Build();
	}

private:
	void Build()
	{
		// Physical parts
		PhysicalGroup<3>* domain = nullptr;
		if (this->PhysicalParts.empty())
			this->PhysicalParts = CubeGeometry::PhysicalParts();
		domain = this->PhysicalParts[0];

		// Vertices
		for (Vertex* cartV : _cartMesh->Vertices)
		{
			MeshVertex<3>* v = new MeshVertex<3>(*cartV);
			this->Vertices.push_back(v);
			this->_verticesByNumber.insert({ v->Number, v });
		}

		// Elements
		BigNumber elemNumber = 0;
		BigNumber faceNumber = 0;
		for (Element<3>* e : _cartMesh->Elements)
		{
			ParallelepipedElement* cube = dynamic_cast<ParallelepipedElement*>(e);

			MeshVertex<3>* v_000 = _verticesByNumber.at(cube->BackLeftBottomCorner->Number);
			MeshVertex<3>* v_100 = _verticesByNumber.at(cube->FrontLeftBottomCorner->Number);
			MeshVertex<3>* v_110 = _verticesByNumber.at(cube->FrontRightBottomCorner->Number);
			MeshVertex<3>* v_111 = _verticesByNumber.at(cube->FrontRightTopCorner->Number);
			MeshVertex<3>* v_101 = _verticesByNumber.at(cube->FrontLeftTopCorner->Number);
			MeshVertex<3>* v_010 = _verticesByNumber.at(cube->BackRightBottomCorner->Number);
			MeshVertex<3>* v_011 = _verticesByNumber.at(cube->BackRightTopCorner->Number);
			MeshVertex<3>* v_001 = _verticesByNumber.at(cube->BackLeftTopCorner->Number);

			TetrahedralElement* tetra0 = this->CreateAndAddNewTetra(elemNumber++, v_010, v_110, v_111, v_100);
			TetrahedralElement* tetra1 = this->CreateAndAddNewTetra(elemNumber++, v_111, v_010, v_100, v_000);
			TetrahedralElement* tetra2 = this->CreateAndAddNewTetra(elemNumber++, v_000, v_100, v_111, v_101);
			TetrahedralElement* tetra3 = this->CreateAndAddNewTetra(elemNumber++, v_111, v_010, v_000, v_011);
			TetrahedralElement* tetra4 = this->CreateAndAddNewTetra(elemNumber++, v_101, v_111, v_000, v_011);
			TetrahedralElement* tetra5 = this->CreateAndAddNewTetra(elemNumber++, v_101, v_000, v_001, v_011);
			tetra0->PhysicalPart = domain;
			tetra1->PhysicalPart = domain;
			tetra2->PhysicalPart = domain;
			tetra3->PhysicalPart = domain;
			tetra4->PhysicalPart = domain;
			tetra5->PhysicalPart = domain;

			vector<TetrahedralElement*> tetrasInCube(6);
			tetrasInCube[0] = tetra0;
			tetrasInCube[1] = tetra1;
			tetrasInCube[2] = tetra2;
			tetrasInCube[3] = tetra3;
			tetrasInCube[4] = tetra4;
			tetrasInCube[5] = tetra5;

			this->_tetrasInCube.insert({ cube, tetrasInCube });

			this->CreateFaceIfNeeded(tetra0, faceNumber, v_010, v_111, v_110);
			this->CreateFaceIfNeeded(tetra0, faceNumber, v_010, v_110, v_100);
			this->CreateInteriorFace(tetra0, tetra1, faceNumber, v_010, v_100, v_111);
			this->CreateFaceIfNeeded(tetra0, faceNumber, v_100, v_110, v_111);
			this->CreateInteriorFace(tetra1, tetra3, faceNumber, v_111, v_010, v_000);
			this->CreateInteriorFace(tetra1, tetra2, faceNumber, v_111, v_000, v_100);
			this->CreateFaceIfNeeded(tetra1, faceNumber, v_000, v_010, v_100);
			this->CreateFaceIfNeeded(tetra2, faceNumber, v_000, v_100, v_101);
			this->CreateInteriorFace(tetra2, tetra4, faceNumber, v_000, v_101, v_111);
			this->CreateFaceIfNeeded(tetra2, faceNumber, v_101, v_100, v_111);
			this->CreateFaceIfNeeded(tetra3, faceNumber, v_111, v_010, v_011);
			this->CreateInteriorFace(tetra3, tetra4, faceNumber, v_111, v_011, v_000);
			this->CreateFaceIfNeeded(tetra3, faceNumber, v_011, v_010, v_000);
			this->CreateFaceIfNeeded(tetra4, faceNumber, v_101, v_111, v_011);
			this->CreateInteriorFace(tetra4, tetra5, faceNumber, v_101, v_011, v_000);
			this->CreateFaceIfNeeded(tetra5, faceNumber, v_101, v_001, v_000);
			this->CreateFaceIfNeeded(tetra5, faceNumber, v_101, v_011, v_001);
			this->CreateFaceIfNeeded(tetra5, faceNumber, v_011, v_000, v_001);
		}

		for (Face<3>* f : this->Faces)
		{
			if (f->IsDomainBoundary)
				this->BoundaryFaces.push_back(f);
		}
	}

public:
	virtual string Description() override
	{
		return "Cartesian tetrahedral";
	}
	virtual string FileNamePart() override
	{
		return "cart-tetra";
	}
	string GeometryDescription() override
	{
		return "Cube";
	}

	void CoarsenMesh(CoarseningStrategy elemCoarseningStgy, FaceCoarseningStrategy faceCoarseningStgy, double coarseningFactor) override
	{
		if (elemCoarseningStgy != CoarseningStrategy::StandardCoarsening)
			TetrahedralMesh::CoarsenMesh(elemCoarseningStgy, faceCoarseningStgy, coarseningFactor);

		if (faceCoarseningStgy != FaceCoarseningStrategy::InterfaceCollapsing)
			Utils::FatalError("Unmanaged face coarsening strategy");

		_cartMesh->CoarsenMesh(elemCoarseningStgy, faceCoarseningStgy, coarseningFactor);
		Cube_CartesianTetrahedralMesh* coarseMesh = new Cube_CartesianTetrahedralMesh(dynamic_cast<Cube_CartesianMesh*>(_cartMesh->CoarseMesh), false);
		this->InitializeCoarsening(coarseMesh);
		coarseMesh->ComesFrom.CS = CoarseningStrategy::StandardCoarsening;
		coarseMesh->Build();

		for (Element<3>* fc : _cartMesh->Elements)
		{
			ParallelepipedElement* fineCube = dynamic_cast<ParallelepipedElement*>(fc);
			vector<TetrahedralElement*> fineTetras = this->_tetrasInCube.at(fineCube);
			assert(fineTetras.size() == 6);

			vector<TetrahedralElement*> coarseTetras = coarseMesh->_tetrasInCube.at(dynamic_cast<ParallelepipedElement*>(fineCube->CoarserElement));
			assert(coarseTetras.size() == 6);

			for (TetrahedralElement* fineTetra : fineTetras)
			{
				for (TetrahedralElement* coarseTetra : coarseTetras)
				{
					if (coarseTetra->Contains(fineTetra->Center()))
					{
						fineTetra->CoarserElement = coarseTetra;
						coarseTetra->FinerElements.push_back(fineTetra);
						break;
					}
				}
				assert(fineTetra->CoarserElement);

			}
		}

		this->LinkFacesToCoarseFaces();

		this->FinalizeCoarsening();
	}

	~Cube_CartesianTetrahedralMesh()
	{
		if (!this->FineMesh)
			delete _cartMesh;
	}
};