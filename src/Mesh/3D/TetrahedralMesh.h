#pragma once
#include "../PolyhedralMesh.h"
#include "Tetrahedron.h"
using namespace std;

class TetrahedralMesh : virtual public PolyhedralMesh<3>
{
protected:
	map<size_t, MeshVertex<3>*> _verticesByNumber;
	double _h = -1;
	double _regularity = 1;
public:
	TetrahedralMesh() : PolyhedralMesh()
	{}

	virtual string Description() override
	{
		return "Tetrahedral";
	}

	virtual string FileNamePart() override
	{
		return "tetrahedral";
	}

	double H() override
	{
		return _h;
	}

	double Regularity() override
	{
		return _regularity;
	}

	virtual void RefineMesh(CoarseningStrategy strategy) override
	{
		if (this->FineMesh)
			assert(false && "Mesh already refined!");

		if (strategy == CoarseningStrategy::BeyRefinement)
			RefineMeshByBeyMethod();
		else
			PolyhedralMesh<3>::RefineMesh(strategy);
	}

protected:
	void RefineMeshByBeyMethod()
	{
		cout << "Mesh refinement using Bey's method" << endl;

		TetrahedralMesh* fineMesh = new TetrahedralMesh();
		fineMesh->ComesFrom.CS = CoarseningStrategy::BeyRefinement;
		fineMesh->ComesFrom.nFineElementsByCoarseElement = 8;
		fineMesh->ComesFrom.nFineFacesAddedByCoarseElement = 8;
		fineMesh->ComesFrom.nFineFacesByKeptCoarseFace = 4;

		this->FineMesh = fineMesh;
		fineMesh->CoarseMesh = this;

		//----------------------------//
		//    Copy coarse vertices    //
		//----------------------------//

		for (Vertex* coarseV : this->Vertices)
		{
			MeshVertex<3>* fineV = new MeshVertex<3>(*coarseV);
			fineMesh->Vertices.push_back(fineV);
			fineMesh->_verticesByNumber.insert({ fineV->Number, fineV });
		}

		//-----------------------------------------------//
		//    Create new vertices, elements and faces    //
		//-----------------------------------------------//

		map<DomPoint, MeshVertex<3>*> newVertices;
		BigNumber vertexNumber = fineMesh->Vertices.size();
		BigNumber elemNumber = 0;
		BigNumber faceNumber = 0;
		for (Element<3>* e : this->Elements)
		{
			Tetrahedron* coarseTetra = dynamic_cast<Tetrahedron*>(e);

			MeshVertex<3>* V0 = fineMesh->_verticesByNumber.at(coarseTetra->V1()->Number);
			MeshVertex<3>* V1 = fineMesh->_verticesByNumber.at(coarseTetra->V2()->Number);
			MeshVertex<3>* V2 = fineMesh->_verticesByNumber.at(coarseTetra->V3()->Number);
			MeshVertex<3>* V3 = fineMesh->_verticesByNumber.at(coarseTetra->V4()->Number);

			DomPoint p01 = Middle<3>(V0, V1);
			DomPoint p02 = Middle<3>(V0, V2);
			DomPoint p03 = Middle<3>(V0, V3);
			DomPoint p12 = Middle<3>(V1, V2);
			DomPoint p13 = Middle<3>(V1, V3);
			DomPoint p23 = Middle<3>(V2, V3);
			MeshVertex<3>* V01 = ExistingNewVertex(newVertices, p01);
			MeshVertex<3>* V02 = ExistingNewVertex(newVertices, p02);
			MeshVertex<3>* V03 = ExistingNewVertex(newVertices, p03);
			MeshVertex<3>* V12 = ExistingNewVertex(newVertices, p12);
			MeshVertex<3>* V13 = ExistingNewVertex(newVertices, p13);
			MeshVertex<3>* V23 = ExistingNewVertex(newVertices, p23);
			if (!V01)
			{
				V01 = new MeshVertex<3>(vertexNumber++, p01);
				fineMesh->Vertices.push_back(V01);
				newVertices.insert({ *V01, V01 });
			}
			if (!V02)
			{
				V02 = new MeshVertex<3>(vertexNumber++, p02);
				fineMesh->Vertices.push_back(V02);
				newVertices.insert({ *V02, V02 });
			}
			if (!V03)
			{
				V03 = new MeshVertex<3>(vertexNumber++, p03);
				fineMesh->Vertices.push_back(V03);
				newVertices.insert({ *V03, V03 });
			}
			if (!V12)
			{
				V12 = new MeshVertex<3>(vertexNumber++, p12);
				fineMesh->Vertices.push_back(V12);
				newVertices.insert({ *V12, V12 });
			}
			if (!V13)
			{
				V13 = new MeshVertex<3>(vertexNumber++, p13);
				fineMesh->Vertices.push_back(V13);
				newVertices.insert({ *V13, V13 });
			}
			if (!V23)
			{
				V23 = new MeshVertex<3>(vertexNumber++, p23);
				fineMesh->Vertices.push_back(V23);
				newVertices.insert({ *V23, V23 });
			}

			// Corners of the tetrahedra
			Tetrahedron* cornerTetra1 = fineMesh->CreateAndAddNewTetra(elemNumber++, V0, V01, V02, V03, coarseTetra);
			Tetrahedron* cornerTetra2 = fineMesh->CreateAndAddNewTetra(elemNumber++, V01, V1, V12, V13, coarseTetra);
			Tetrahedron* cornerTetra3 = fineMesh->CreateAndAddNewTetra(elemNumber++, V02, V12, V2, V23, coarseTetra);
			Tetrahedron* cornerTetra4 = fineMesh->CreateAndAddNewTetra(elemNumber++, V03, V13, V23, V3, coarseTetra);

			// Octahedron
			Tetrahedron* octaTetra1 = fineMesh->CreateAndAddNewTetra(elemNumber++, V01, V02, V03, V13, coarseTetra);
			Tetrahedron* octaTetra2 = fineMesh->CreateAndAddNewTetra(elemNumber++, V01, V02, V12, V13, coarseTetra);
			Tetrahedron* octaTetra3 = fineMesh->CreateAndAddNewTetra(elemNumber++, V02, V03, V13, V23, coarseTetra);
			Tetrahedron* octaTetra4 = fineMesh->CreateAndAddNewTetra(elemNumber++, V02, V12, V13, V23, coarseTetra);

			//--------------------//
			//    Create faces    //
			//--------------------//

			// Exterior faces
			// cornerTetra1
			fineMesh->CreateFaceIfNeeded(cornerTetra1, faceNumber, V0, V01, V02);
			fineMesh->CreateFaceIfNeeded(cornerTetra1, faceNumber, V03, V0, V01);
			fineMesh->CreateFaceIfNeeded(cornerTetra1, faceNumber, V02, V03, V0);

			// cornerTetra2
			fineMesh->CreateFaceIfNeeded(cornerTetra2, faceNumber, V1, V12, V13);
			fineMesh->CreateFaceIfNeeded(cornerTetra2, faceNumber, V01, V1, V12);
			fineMesh->CreateFaceIfNeeded(cornerTetra2, faceNumber, V13, V01, V1);

			// cornerTetra3
			fineMesh->CreateFaceIfNeeded(cornerTetra3, faceNumber, V2, V23, V02);
			fineMesh->CreateFaceIfNeeded(cornerTetra3, faceNumber, V12, V2, V23);
			fineMesh->CreateFaceIfNeeded(cornerTetra3, faceNumber, V02, V12, V2);

			// cornerTetra4
			fineMesh->CreateFaceIfNeeded(cornerTetra4, faceNumber, V3, V03, V13);
			fineMesh->CreateFaceIfNeeded(cornerTetra4, faceNumber, V23, V3, V03);
			fineMesh->CreateFaceIfNeeded(cornerTetra4, faceNumber, V13, V23, V3);

			// octaTetra
			fineMesh->CreateFaceIfNeeded(octaTetra1, faceNumber, V03, V13, V01);
			fineMesh->CreateFaceIfNeeded(octaTetra2, faceNumber, V01, V02, V12);
			fineMesh->CreateFaceIfNeeded(octaTetra3, faceNumber, V23, V02, V03);
			fineMesh->CreateFaceIfNeeded(octaTetra4, faceNumber, V12, V13, V23);

			// Interior faces
			fineMesh->CreateInteriorFace(cornerTetra1, octaTetra1, faceNumber, V01, V02, V03);
			fineMesh->CreateInteriorFace(cornerTetra2, octaTetra2, faceNumber, V12, V13, V01);
			fineMesh->CreateInteriorFace(cornerTetra3, octaTetra4, faceNumber, V23, V02, V12);
			fineMesh->CreateInteriorFace(cornerTetra4, octaTetra3, faceNumber, V03, V13, V23);
			fineMesh->CreateInteriorFace(octaTetra1, octaTetra3, faceNumber, V02, V03, V13);
			fineMesh->CreateInteriorFace(octaTetra1, octaTetra2, faceNumber, V13, V01, V02);
			fineMesh->CreateInteriorFace(octaTetra2, octaTetra4, faceNumber, V02, V12, V13);
			fineMesh->CreateInteriorFace(octaTetra3, octaTetra4, faceNumber, V13, V23, V02);
		}
		
		newVertices.clear();

		for (Face<3>* f : fineMesh->Faces)
		{
			if (f->IsDomainBoundary)
				fineMesh->BoundaryFaces.push_back(f);
		}

		//------------------------------//
		//    Link coarse/fine faces    //
		//------------------------------//

		fineMesh->LinkFacesToCoarseFaces();
	}

	MeshVertex<3>* ExistingNewVertex(map<DomPoint, MeshVertex<3>*> &newVertices, DomPoint p)
	{
		map<DomPoint, MeshVertex<3>*>::iterator it = newVertices.find(p);
		if (it != newVertices.end())
			return it->second;
		return nullptr;
	}

	Tetrahedron* CreateAndAddNewTetra(BigNumber elemNumber, MeshVertex<3>* V1, MeshVertex<3>* V2, MeshVertex<3>* V3, MeshVertex<3>* V4, Tetrahedron* coarseTetra = nullptr)
	{
		Tetrahedron* tetra = new Tetrahedron(elemNumber, V1, V2, V3, V4);

		V1->Elements.push_back(tetra);
		V2->Elements.push_back(tetra);
		V3->Elements.push_back(tetra);
		V4->Elements.push_back(tetra);

		if (coarseTetra)
		{
			coarseTetra->FinerElements.push_back(tetra);
			tetra->CoarserElement = coarseTetra;
		}

		this->Elements.push_back(tetra);

		if (tetra->Diameter() > this->_h)
			this->_h = tetra->Diameter();

		if (tetra->Regularity() < this->_regularity)
			this->_regularity = tetra->Regularity();

		return tetra;
	}

	void CreateFaceIfNeeded(Tetrahedron* tetra, BigNumber& faceNumber, MeshVertex<3>* V1, MeshVertex<3>* V2, MeshVertex<3>* V3)
	{
		Face<3>* face = this->ExistingFaceWithVertices(vector<MeshVertex<3>*> {V1, V2, V3});
		if (!face)
		{
			face = new TriangularFace(faceNumber++, V1, V2, V3, tetra);
			face->IsDomainBoundary = true;
			this->Faces.push_back(face);

			V1->Faces.push_back(face);
			V2->Faces.push_back(face);
			V3->Faces.push_back(face);
		}
		else
		{
			face->Element2 = tetra;
			face->IsDomainBoundary = false;
			this->InteriorFaces.push_back(face);
		}
		tetra->AddFace(face);
	}

	void CreateInteriorFace(Tetrahedron* tetra1, Tetrahedron* tetra2, BigNumber& faceNumber, MeshVertex<3>* V1, MeshVertex<3>* V2, MeshVertex<3>* V3)
	{
		Face<3>* face = new TriangularFace(faceNumber++, V1, V2, V3, tetra1, tetra2);
		face->IsDomainBoundary = false;
		this->Faces.push_back(face);
		this->InteriorFaces.push_back(face);

		V1->Faces.push_back(face);
		V2->Faces.push_back(face);
		V3->Faces.push_back(face);

		tetra1->AddFace(face);
		tetra2->AddFace(face);
	}
};