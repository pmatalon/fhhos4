#pragma once
#include "../PolyhedralMesh.h"
#include "Tetrahedron.h"
using namespace std;

class TetrahedralMesh : virtual public PolyhedralMesh<3>
{
private:
	map<size_t, MeshVertex<3>*> _verticesByNumber;
	double _h = -1;
	double _regularity = 1;
public:
	TetrahedralMesh() : PolyhedralMesh()
	{}

public:
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

	virtual void RefineMesh(CoarseningStrategy strategy)
	{
		if (this->FineMesh)
			assert(false && "Mesh already refined!");

		if (strategy == CoarseningStrategy::BeyRefinement)
			RefineMeshByBeyMethod();
		else
			assert(false && "Unmanaged refinement strategy");
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
			MeshVertex<3>* V01 = ExistingNewVertex(fineMesh, p01);
			MeshVertex<3>* V02 = ExistingNewVertex(fineMesh, p02);
			MeshVertex<3>* V03 = ExistingNewVertex(fineMesh, p03);
			MeshVertex<3>* V12 = ExistingNewVertex(fineMesh, p12);
			MeshVertex<3>* V13 = ExistingNewVertex(fineMesh, p13);
			MeshVertex<3>* V23 = ExistingNewVertex(fineMesh, p23);
			if (!V01)
			{
				V01 = new MeshVertex<3>(vertexNumber++, p01);
				fineMesh->Vertices.push_back(V01);
			}
			if (!V02)
			{
				V02 = new MeshVertex<3>(vertexNumber++, p02);
				fineMesh->Vertices.push_back(V02);
			}
			if (!V03)
			{
				V03 = new MeshVertex<3>(vertexNumber++, p03);
				fineMesh->Vertices.push_back(V03);
			}
			if (!V12)
			{
				V12 = new MeshVertex<3>(vertexNumber++, p12);
				fineMesh->Vertices.push_back(V12);
			}
			if (!V13)
			{
				V13 = new MeshVertex<3>(vertexNumber++, p13);
				fineMesh->Vertices.push_back(V13);
			}
			if (!V23)
			{
				V23 = new MeshVertex<3>(vertexNumber++, p23);
				fineMesh->Vertices.push_back(V23);
			}

			// Corners of the tetrahedra
			Tetrahedron* cornerTetra1 = CreateAndAddNewTetra(coarseTetra, fineMesh, elemNumber++, V0, V01, V02, V03);
			Tetrahedron* cornerTetra2 = CreateAndAddNewTetra(coarseTetra, fineMesh, elemNumber++, V01, V1, V12, V13);
			Tetrahedron* cornerTetra3 = CreateAndAddNewTetra(coarseTetra, fineMesh, elemNumber++, V02, V12, V2, V23);
			Tetrahedron* cornerTetra4 = CreateAndAddNewTetra(coarseTetra, fineMesh, elemNumber++, V03, V13, V23, V3);

			// Octahedron
			Tetrahedron* octaTetra1 = CreateAndAddNewTetra(coarseTetra, fineMesh, elemNumber++, V01, V02, V03, V13);
			Tetrahedron* octaTetra2 = CreateAndAddNewTetra(coarseTetra, fineMesh, elemNumber++, V01, V02, V12, V13);
			Tetrahedron* octaTetra3 = CreateAndAddNewTetra(coarseTetra, fineMesh, elemNumber++, V02, V03, V13, V23);
			Tetrahedron* octaTetra4 = CreateAndAddNewTetra(coarseTetra, fineMesh, elemNumber++, V02, V12, V13, V23);

			//--------------------//
			//    Create faces    //
			//--------------------//

			// Exterior faces
			// cornerTetra1
			CreateExteriorFaceIfNeeded(cornerTetra1, fineMesh, faceNumber, V0, V01, V02);
			CreateExteriorFaceIfNeeded(cornerTetra1, fineMesh, faceNumber, V03, V0, V01);
			CreateExteriorFaceIfNeeded(cornerTetra1, fineMesh, faceNumber, V02, V03, V0);

			// cornerTetra2
			CreateExteriorFaceIfNeeded(cornerTetra2, fineMesh, faceNumber, V1, V12, V13);
			CreateExteriorFaceIfNeeded(cornerTetra2, fineMesh, faceNumber, V01, V1, V12);
			CreateExteriorFaceIfNeeded(cornerTetra2, fineMesh, faceNumber, V13, V01, V1);

			// cornerTetra3
			CreateExteriorFaceIfNeeded(cornerTetra3, fineMesh, faceNumber, V2, V23, V02);
			CreateExteriorFaceIfNeeded(cornerTetra3, fineMesh, faceNumber, V12, V2, V23);
			CreateExteriorFaceIfNeeded(cornerTetra3, fineMesh, faceNumber, V02, V12, V2);

			// cornerTetra4
			CreateExteriorFaceIfNeeded(cornerTetra4, fineMesh, faceNumber, V3, V03, V13);
			CreateExteriorFaceIfNeeded(cornerTetra4, fineMesh, faceNumber, V23, V3, V03);
			CreateExteriorFaceIfNeeded(cornerTetra4, fineMesh, faceNumber, V13, V23, V3);

			// octaTetra
			CreateExteriorFaceIfNeeded(octaTetra1, fineMesh, faceNumber, V03, V13, V01);
			CreateExteriorFaceIfNeeded(octaTetra2, fineMesh, faceNumber, V01, V02, V12);
			CreateExteriorFaceIfNeeded(octaTetra3, fineMesh, faceNumber, V23, V02, V03);
			CreateExteriorFaceIfNeeded(octaTetra4, fineMesh, faceNumber, V12, V13, V23);

			// Interior faces
			CreateInteriorFace(cornerTetra1, octaTetra1, fineMesh, faceNumber, V01, V02, V03);
			CreateInteriorFace(cornerTetra2, octaTetra2, fineMesh, faceNumber, V12, V13, V01);
			CreateInteriorFace(cornerTetra3, octaTetra4, fineMesh, faceNumber, V23, V02, V12);
			CreateInteriorFace(cornerTetra4, octaTetra3, fineMesh, faceNumber, V03, V13, V23);
			CreateInteriorFace(octaTetra1, octaTetra3, fineMesh, faceNumber, V02, V03, V13);
			CreateInteriorFace(octaTetra1, octaTetra2, fineMesh, faceNumber, V13, V01, V02);
			CreateInteriorFace(octaTetra2, octaTetra4, fineMesh, faceNumber, V02, V12, V13);
			CreateInteriorFace(octaTetra3, octaTetra4, fineMesh, faceNumber, V13, V23, V02);
		}

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

	MeshVertex<3>* ExistingNewVertex(Mesh<3>* fineMesh, DomPoint p)
	{
		for (BigNumber i = this->Vertices.size(); i < fineMesh->Vertices.size(); i++)
		{
			MeshVertex<3>* newVertex = static_cast<MeshVertex<3>*>(fineMesh->Vertices[i]);
			if (*newVertex == p)
				return newVertex;
		}
		return nullptr;
	}

	Tetrahedron* CreateAndAddNewTetra(Tetrahedron* coarseTetra, TetrahedralMesh* fineMesh, BigNumber elemNumber, MeshVertex<3>* V1, MeshVertex<3>* V2, MeshVertex<3>* V3, MeshVertex<3>* V4)
	{
		Tetrahedron* fineTetra = new Tetrahedron(elemNumber, V1, V2, V3, V4);

		V1->Elements.push_back(fineTetra);
		V2->Elements.push_back(fineTetra);
		V3->Elements.push_back(fineTetra);
		V4->Elements.push_back(fineTetra);

		coarseTetra->FinerElements.push_back(fineTetra);
		fineTetra->CoarserElement = coarseTetra;

		fineMesh->Elements.push_back(fineTetra);

		if (fineTetra->Diameter() > fineMesh->_h)
			fineMesh->_h = fineTetra->Diameter();

		if (fineTetra->Regularity() < fineMesh->_regularity)
			fineMesh->_regularity = fineTetra->Regularity();

		return fineTetra;
	}

	void CreateExteriorFaceIfNeeded(Tetrahedron* fineTetra, Mesh<3>* fineMesh, BigNumber& faceNumber, MeshVertex<3>* V1, MeshVertex<3>* V2, MeshVertex<3>* V3)
	{
		Face<3>* fineFace = fineMesh->ExistingFaceWithVertices(vector<MeshVertex<3>*> {V1, V2, V3});
		if (!fineFace)
		{
			fineFace = new TriangularFace(faceNumber++, V1, V2, V3, fineTetra);
			fineFace->IsDomainBoundary = true;
			fineMesh->Faces.push_back(fineFace);

			V1->Faces.push_back(fineFace);
			V2->Faces.push_back(fineFace);
			V3->Faces.push_back(fineFace);
		}
		else
		{
			fineFace->Element2 = fineTetra;
			fineFace->IsDomainBoundary = false;
			fineMesh->InteriorFaces.push_back(fineFace);
		}
		fineTetra->AddFace(fineFace);
	}

	void CreateInteriorFace(Tetrahedron* fineTetra1, Tetrahedron* fineTetra2, Mesh<3>* fineMesh, BigNumber& faceNumber, MeshVertex<3>* V1, MeshVertex<3>* V2, MeshVertex<3>* V3)
	{
		Face<3>* fineFace = new TriangularFace(faceNumber++, V1, V2, V3, fineTetra1, fineTetra2);
		fineFace->IsDomainBoundary = false;
		fineMesh->Faces.push_back(fineFace);
		fineMesh->InteriorFaces.push_back(fineFace);

		V1->Faces.push_back(fineFace);
		V2->Faces.push_back(fineFace);
		V3->Faces.push_back(fineFace);

		fineTetra1->AddFace(fineFace);
		fineTetra2->AddFace(fineFace);
	}
};