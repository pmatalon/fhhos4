#pragma once
#include "TetrahedralMesh.h"
#include "CartesianGrid3D.h"
using namespace std;

class CartesianTetrahedralMesh : public TetrahedralMesh
{
private:
	CartesianGrid3D* _cartMesh = nullptr;
	map<Parallelepiped*, vector<Tetrahedron*>> _tetrasInCube;

public:
	CartesianTetrahedralMesh(BigNumber n) : CartesianTetrahedralMesh(new CartesianGrid3D(n, n, n))
	{}

	CartesianTetrahedralMesh(CartesianGrid3D* cartMesh) : TetrahedralMesh()
	{
		_cartMesh = cartMesh;

		for (Vertex* cartV : _cartMesh->Vertices)
		{
			MeshVertex<3>* v = new MeshVertex<3>(*cartV);
			this->Vertices.push_back(v);
			this->_verticesByNumber.insert({ v->Number, v });
		}

		BigNumber elemNumber = 0;
		BigNumber faceNumber = 0;
		for (Element<3>* e : _cartMesh->Elements)
		{
			Parallelepiped* cube = dynamic_cast<Parallelepiped*>(e);

			MeshVertex<3>* v_000 = _verticesByNumber.at(cube->BackLeftBottomCorner->Number);
			MeshVertex<3>* v_100 = _verticesByNumber.at(cube->FrontLeftBottomCorner->Number);
			MeshVertex<3>* v_110 = _verticesByNumber.at(cube->FrontRightBottomCorner->Number);
			MeshVertex<3>* v_111 = _verticesByNumber.at(cube->FrontRightTopCorner->Number);
			MeshVertex<3>* v_101 = _verticesByNumber.at(cube->FrontLeftTopCorner->Number);
			MeshVertex<3>* v_010 = _verticesByNumber.at(cube->BackRightBottomCorner->Number);
			MeshVertex<3>* v_011 = _verticesByNumber.at(cube->BackRightTopCorner->Number);
			MeshVertex<3>* v_001 = _verticesByNumber.at(cube->BackLeftTopCorner->Number);

			Tetrahedron* tetra0 = this->CreateAndAddNewTetra(elemNumber++, v_010, v_110, v_111, v_100);
			Tetrahedron* tetra1 = this->CreateAndAddNewTetra(elemNumber++, v_111, v_010, v_100, v_000);
			Tetrahedron* tetra2 = this->CreateAndAddNewTetra(elemNumber++, v_000, v_100, v_111, v_101);
			Tetrahedron* tetra3 = this->CreateAndAddNewTetra(elemNumber++, v_111, v_010, v_000, v_011);
			Tetrahedron* tetra4 = this->CreateAndAddNewTetra(elemNumber++, v_101, v_111, v_000, v_011);
			Tetrahedron* tetra5 = this->CreateAndAddNewTetra(elemNumber++, v_101, v_000, v_001, v_011);

			vector<Tetrahedron*> tetrasInCube(6);
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

	virtual string Description() override
	{
		return "Cartesian tetrahedral";
	}

	virtual string FileNamePart() override
	{
		return "carttetra";
	}

	void CoarsenMesh(CoarseningStrategy strategy) override
	{
		if (strategy != CoarseningStrategy::StandardCoarsening)
			TetrahedralMesh::CoarsenMesh(strategy);

		_cartMesh->CoarsenMesh(strategy);
		CartesianTetrahedralMesh* coarseMesh = new CartesianTetrahedralMesh(dynamic_cast<CartesianGrid3D*>(_cartMesh->CoarseMesh));
		this->CoarseMesh = coarseMesh;
		coarseMesh->FineMesh = this;
		coarseMesh->ComesFrom.CS = CoarseningStrategy::StandardCoarsening;

		for (Element<3>* fc : _cartMesh->Elements)
		{
			Parallelepiped* fineCube = dynamic_cast<Parallelepiped*>(fc);
			vector<Tetrahedron*> fineTetras = this->_tetrasInCube.at(fineCube);
			assert(fineTetras.size() == 6);

			vector<Tetrahedron*> coarseTetras = coarseMesh->_tetrasInCube.at(dynamic_cast<Parallelepiped*>(fineCube->CoarserElement));
			assert(coarseTetras.size() == 6);

			for (Tetrahedron* fineTetra : fineTetras)
			{
				for (Tetrahedron* coarseTetra : coarseTetras)
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

		this->CoarseMesh->SetDiffusionCoefficient(this->_diffusionPartition);
		this->CoarseMesh->SetBoundaryConditions(this->_boundaryConditions);
	}

	~CartesianTetrahedralMesh()
	{
		if (!this->FineMesh)
			delete _cartMesh;
	}
};