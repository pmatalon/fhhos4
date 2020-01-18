#pragma once
#include <gmsh.h>
#include "PolyhedralMesh.h"
#include "2D/Quadrilateral.h"
#include "3D/Tetrahedron.h"
using namespace std;

enum GMSHElementTypes
{
	GMSH_Segment = 1,
	GMSH_Triangle = 2,
	GMSH_Quadrilateral = 3,
	GMSH_Tetrahedron = 4,
	GMSH_Hexahedron = 8
};

enum GMSHFaceTypes
{
	GMSH_TriangleFace = 3,
	GMSH_QuadrilateralFace = 4
};

template <int Dim>
class MeshVertex : public Vertex
{
public:
	vector<Element<Dim>*> Elements;
	vector<Face<Dim>*> Faces;

	MeshVertex(BigNumber number, double x) : Vertex(number, x) {}
	MeshVertex(BigNumber number, double x, double y) : Vertex(number, x, y) {}
	MeshVertex(BigNumber number, double x, double y, double z) : Vertex(number, x, y, z) {}
};

template <int Dim>
class GMSHMesh : public PolyhedralMesh<Dim>
{
protected:
	string _description = "GMSH file";
	string _fileNamePart = "gmsh-file";
private:
	map<size_t, BigNumber> _elementExternalNumbers;
	map<size_t, BigNumber> _vertexExternalNumbers;
	double _h = -1;
public:
	GMSHMesh(string mshFile, int initialRefinements = 0) : PolyhedralMesh<Dim>()
	{
		if (!Utils::FileExists(mshFile))
		{
			if (Utils::FileExists(Mesh<Dim>::MeshDirectory + mshFile))
				mshFile = Mesh<Dim>::MeshDirectory + mshFile;
			else
			{
				cout << Utils::BeginRed << "File not found: " << mshFile << Utils::EndColor;
				exit(EXIT_FAILURE);
			}
		}

		gmsh::initialize();
		//gmsh::option::setNumber("General.Terminal", 1);
		//gmsh::option::setNumber("General.Verbosity", 99);

		gmsh::open(mshFile);

		for (int i = 0; i < initialRefinements; i++)
			gmsh::model::mesh::refine();

		Build();
	}

private:
	GMSHMesh() : PolyhedralMesh<Dim>()
	{
		Build();
	}

	// Dim-specific functions
	Element<Dim>* CreateElement(int elemType, size_t elementTag, const vector<size_t>& elementNodes, size_t& elemNodeIndex, BigNumber elemNumber) { return nullptr; }
	void CreateFaces(int elemType, BigNumber& faceNumber) { }

	void Build()
	{
		//----------//
		// Vertices //
		//----------//

		vector<size_t> nodeTags;
		vector<double> coord;
		vector<double> parametricCoord;
		bool returnParametricCoord = false;
		bool includeBoundary = true;
		gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord, Dim, -1, includeBoundary, returnParametricCoord);

		// If .geo file, the mesh hasn't been generated
		if (nodeTags.empty())
		{
			gmsh::model::mesh::generate(Dim);
			gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord, Dim, -1, includeBoundary, returnParametricCoord);
			if (nodeTags.empty())
			{
				cout << Utils::BeginRed << "Error: the mesh generation failed. Try commenting the line SetFactory(\"OpenCASCADE\"); (if there is one)." << Utils::EndColor;
				exit(EXIT_FAILURE);
			}
		}

		this->Vertices.reserve(nodeTags.size());
		BigNumber vertexNumber = 0;
		for (size_t i = 0; i < nodeTags.size(); i++)
		{
			if (_vertexExternalNumbers.find(nodeTags[i]) == _vertexExternalNumbers.end())
			{
				MeshVertex<Dim>* v = new MeshVertex<Dim>(vertexNumber, coord[3 * i], coord[3 * i + 1], coord[3 * i + 2]);
				_vertexExternalNumbers.insert({ nodeTags[i], vertexNumber });
				this->Vertices.push_back(v);
				vertexNumber++;
			}
		}

		//----------//
		// Elements //
		//----------//

		vector<int> elementTypes;
		vector<vector<size_t>> elementTags;
		vector<vector<size_t>> elemNodeTags;
		gmsh::model::mesh::getElements(elementTypes, elementTags, elemNodeTags, Dim, -1);

		BigNumber elemNumber = 0;
		for (int i = 0; i < elementTypes.size(); i++)
		{
			int elemType = elementTypes[i];
			vector<size_t> elements = elementTags[i];
			vector<size_t> elementNodes = elemNodeTags[i];

			size_t elemNodeIndex = 0;
			for (size_t j = 0; j < elements.size(); j++)
			{
				Element<Dim>* e = CreateElement(elemType, elements[j], elementNodes, elemNodeIndex, elemNumber);

				for (Vertex* v : e->Shape()->Vertices())
				{
					MeshVertex<Dim>* mv = (MeshVertex<Dim>*)v;
					mv->Elements.push_back(e);
				}

				this->Elements.push_back(e);
				elemNumber++;

				if (e->Diameter() > this->_h)
					this->_h = e->Diameter();
			}
		}

		//-----------//
		//   Faces   //
		//-----------//

		BigNumber faceNumber = 0;
		for (int i = 0; i < elementTypes.size(); i++)
		{
			int elemType = elementTypes[i];
			CreateFaces(elemType, faceNumber);
		}
	}

private:
	inline Vertex* GetVertex(BigNumber nodeTag)
	{
		BigNumber internalNumber = _vertexExternalNumbers.at(nodeTag);
		return this->Vertices[internalNumber];
	}

	inline Element<Dim>* GetElement(BigNumber elementTag)
	{
		BigNumber internalNumber = _elementExternalNumbers.at(elementTag);
		return this->Elements[internalNumber];
	}

public:
	~GMSHMesh()
	{
		if (!this->CoarseMesh)
			gmsh::finalize();
	}

	string Description() override
	{
		return _description;
	}

	string FileNamePart() override
	{
		return _fileNamePart;
	}

	double H() override
	{
		return _h;
	}

	void CoarsenMesh(CoarseningStrategy strategy) override
	{
		if (strategy != CoarseningStrategy::StructuredRefinement)
			assert(false && "Unmanaged coarsening strategy");
	}

	virtual void RefineMesh()
	{
		if (this->FineMesh)
			assert(false && "Mesh already refined!");

		cout << "Mesh refinement" << endl;

		string coarse_mesh_tmp_file = "./temporary_coarse.msh";
		string fine_mesh_tmp_file = "./temporary_fine.msh";

		// Save the current mesh in a temporary file, because the refinement will delete it
		gmsh::write(coarse_mesh_tmp_file);

		// Mesh refinement
		gmsh::model::mesh::refine();

		// Building our own mesh objects from the GMSH ones
		GMSHMesh<Dim>* fineMesh = new GMSHMesh<Dim>();
		fineMesh->ComesFrom.CS = CoarseningStrategy::StructuredRefinement;
		fineMesh->_description = this->_description;
		fineMesh->_fileNamePart = this->_fileNamePart;

		this->FineMesh = fineMesh;
		this->FineMesh->CoarseMesh = this;

		// Save the fine mesh in a temporary file because we're going to reload the coarse one
		gmsh::write(fine_mesh_tmp_file);

		// Reloading the now coarse mesh to be able to link its objects with the finer ones
		gmsh::open(coarse_mesh_tmp_file);

		// Linking fine/coarse elements
		for (Element<Dim>* fine : fineMesh->Elements)
		{
			Element<Dim>* coarse = LocateElementThatEmbeds(fine);
			coarse->FinerElements.push_back(fine);
			fine->CoarserElement = coarse;
		}

		// Remove coarse temporary file and reload the fine mesh
		remove(coarse_mesh_tmp_file.c_str());
		gmsh::open(fine_mesh_tmp_file);
		remove(fine_mesh_tmp_file.c_str());

		// Continue linking
		for (Face<Dim>* fineFace : fineMesh->Faces)
		{
			if (fineFace->IsDomainBoundary)
			{
				fineFace->IsRemovedOnCoarserGrid = false;
				for (Face<Dim>* coarseFace : fineFace->Element1->CoarserElement->Faces)
				{
					if (coarseFace->IsDomainBoundary && coarseFace->Contains(fineFace->Center()))
					{
						fineFace->CoarseFace = coarseFace;
						coarseFace->FinerFaces.push_back(fineFace);
						break;
					}
				}
				if (!fineFace->CoarseFace)
					assert(false && "A coarse face should have been found.");
			}
			else
			{
				Element<Dim>* coarseElement1 = fineFace->Element1->CoarserElement;
				Element<Dim>* coarseElement2 = fineFace->Element2->CoarserElement;
				fineFace->IsRemovedOnCoarserGrid = coarseElement1 == coarseElement2;
				if (!fineFace->IsRemovedOnCoarserGrid)
				{
					// TODO: change InterfaceWith so it returns a vector of faces (local refinement)
					Face<Dim>* coarseFace = coarseElement1->InterfaceWith(coarseElement2);
					assert(coarseFace != nullptr);
					fineFace->CoarseFace = coarseFace;
					coarseFace->FinerFaces.push_back(fineFace);
				}
				else
					coarseElement1->FinerFacesRemoved.push_back(fineFace);
			}
		}
	}

	GMSHMesh<Dim>* RefineNTimes(int nRefinements)
	{
		GMSHMesh<Dim>* mesh = this;
		for (int i = 0; i < nRefinements; i++)
		{
			mesh->RefineMesh();
			mesh = static_cast<GMSHMesh<Dim>*>(mesh->FineMesh);
		}
		return mesh;
	}

	GMSHMesh<Dim>* RefineUntilNElements(BigNumber nElements)
	{
		GMSHMesh<Dim>* mesh = this;
		while (mesh->Elements.size() < nElements)
		{
			mesh->RefineMesh();
			mesh = static_cast<GMSHMesh<Dim>*>(mesh->FineMesh);
		}
		return mesh;
	}

private:
	Element<Dim>* LocateElementThatEmbeds(Element<Dim>* finerElement)
	{
		size_t coarseElementTag;
		int coarseElementType;
		DomPoint fineCenter = finerElement->Center();
		vector<size_t> coarseElementNodes;
		double u, v, w;
		bool strictResearch = true;

		vector<std::string> names;
		gmsh::model::list(names);

		gmsh::model::mesh::getElementByCoordinates(fineCenter.X, fineCenter.Y, fineCenter.Z, coarseElementTag,
			coarseElementType, coarseElementNodes, u, v, w, // useless parameters for us
			Dim, strictResearch);

		return GetElement(coarseElementTag);
	}

};

//-------------//
// 2D elements //
//-------------//

template <>
Element<2>* GMSHMesh<2>::CreateElement(int elemType, size_t elementTag, const vector<size_t>& elementNodes, size_t& elemNodeIndex, BigNumber elemNumber)
{ 
	Element<2>* e = nullptr;
	if (elemType == GMSH_Quadrilateral)
	{
		size_t nodeTag1 = elementNodes[elemNodeIndex];
		size_t nodeTag2 = elementNodes[elemNodeIndex + 1];
		size_t nodeTag3 = elementNodes[elemNodeIndex + 2];
		size_t nodeTag4 = elementNodes[elemNodeIndex + 3];

		e = new Quadrilateral(elemNumber, GetVertex(nodeTag1), GetVertex(nodeTag2), GetVertex(nodeTag3), GetVertex(nodeTag4));
		_elementExternalNumbers.insert({ elementTag, elemNumber });
		elemNodeIndex += 4;
	}
	else if (elemType == GMSH_Triangle)
	{
		size_t nodeTag1 = elementNodes[elemNodeIndex];
		size_t nodeTag2 = elementNodes[elemNodeIndex + 1];
		size_t nodeTag3 = elementNodes[elemNodeIndex + 2];

		e = new Triangle(elemNumber, GetVertex(nodeTag1), GetVertex(nodeTag2), GetVertex(nodeTag3));
		_elementExternalNumbers.insert({ elementTag, elemNumber });
		elemNodeIndex += 3;
	}
	else
		assert(false && "GMSH element type not managed.");

	return e;
}

//-------------//
// 3D elements //
//-------------//

template <>
Element<3>* GMSHMesh<3>::CreateElement(int elemType, size_t elementTag, const vector<size_t>& elementNodes, size_t& elemNodeIndex, BigNumber elemNumber)
{
	Element<3>* e = nullptr;
	if (elemType == GMSH_Tetrahedron)
	{
		size_t nodeTag1 = elementNodes[elemNodeIndex];
		size_t nodeTag2 = elementNodes[elemNodeIndex + 1];
		size_t nodeTag3 = elementNodes[elemNodeIndex + 2];
		size_t nodeTag4 = elementNodes[elemNodeIndex + 3];

		e = new Tetrahedron(elemNumber, GetVertex(nodeTag1), GetVertex(nodeTag2), GetVertex(nodeTag3), GetVertex(nodeTag4));
		_elementExternalNumbers.insert({ elementTag, elemNumber });
		elemNodeIndex += 4;
	}
	else
		assert(false && "GMSH element type not managed.");
	return e;
}

//--------------//
//   2D faces   //
//--------------//

template <>
void GMSHMesh<2>::CreateFaces(int elemType, BigNumber& faceNumber)
{
	vector<size_t> edgeNodes;
	bool onlyPrimaryNodes = true;
	gmsh::model::mesh::getElementEdgeNodes(elemType, edgeNodes, -1, onlyPrimaryNodes);

	for (size_t j = 0; j < edgeNodes.size(); j += 2)
	{
		MeshVertex<2>* v1 = (MeshVertex<2>*)GetVertex(edgeNodes[j]);
		MeshVertex<2>* v2 = (MeshVertex<2>*)GetVertex(edgeNodes[j + 1]);

		bool edgeAlreadyExists = false;
		for (Face<2>* f : v1->Faces)
		{
			if (f->HasVertex(v2))
			{
				edgeAlreadyExists = true;
				break;
			}
		}
		if (edgeAlreadyExists)
			continue;

		vector<Element<2>*> neighbours;
		for (Element<2>* e1 : v1->Elements)
		{
			for (Element<2>* e2 : v2->Elements)
			{
				if (e1 == e2)
				{
					neighbours.push_back(e1);
					break;
				}
			}
		}

		Edge* edge;
		if (neighbours.size() == 1)
		{
			edge = new Edge(faceNumber++, v1, v2, neighbours[0]);
			neighbours[0]->AddFace(edge);
			this->BoundaryFaces.push_back(edge);
		}
		else if (neighbours.size() == 2)
		{
			edge = new Edge(faceNumber++, v1, v2, neighbours[0], neighbours[1]);
			neighbours[0]->AddFace(edge);
			neighbours[1]->AddFace(edge);
			this->InteriorFaces.push_back(edge);
		}
		else
			assert(false);
		this->Faces.push_back(edge);
		v1->Faces.push_back(edge);
		v2->Faces.push_back(edge);
	}
}

//--------------//
//   3D faces   //
//--------------//

template <>
void GMSHMesh<3>::CreateFaces(int elemType, BigNumber& faceNumber)
{
	assert(elemType == GMSHElementTypes::GMSH_Tetrahedron);

	int faceType = GMSHFaceTypes::GMSH_TriangleFace;
	vector<size_t> faceNodes;
	bool onlyPrimaryNodes = true;
	gmsh::model::mesh::getElementFaceNodes(elemType, faceType, faceNodes, -1, onlyPrimaryNodes);

	for (size_t j = 0; j < faceNodes.size(); j += 3)
	{
		MeshVertex<3>* v1 = (MeshVertex<3>*)GetVertex(faceNodes[j]);
		MeshVertex<3>* v2 = (MeshVertex<3>*)GetVertex(faceNodes[j + 1]);
		MeshVertex<3>* v3 = (MeshVertex<3>*)GetVertex(faceNodes[j + 2]);

		bool faceAlreadyExists = false;
		for (Face<3>* f : v1->Faces)
		{
			if (f->HasVertex(v2) && f->HasVertex(v3))
			{
				faceAlreadyExists = true;
				break;
			}
		}
		if (faceAlreadyExists)
			continue;

		vector<Element<3>*> neighbours;
		for (Element<3>* f1 : v1->Elements)
		{
			for (Element<3>* f2 : v2->Elements)
			{
				if (f1 == f2)
				{
					for (Element<3>* f3 : v3->Elements)
					{
						if (f2 == f3)
						{
							neighbours.push_back(f1);
							break;
						}
					}
					break;
				}
			}
		}

		TriangularFace* face;
		if (neighbours.size() == 1)
		{
			face = new TriangularFace(faceNumber++, v1, v2, v3, neighbours[0]);
			neighbours[0]->AddFace(face);
			this->BoundaryFaces.push_back(face);
		}
		else if (neighbours.size() == 2)
		{
			face = new TriangularFace(faceNumber++, v1, v2, v3, neighbours[0], neighbours[1]);
			neighbours[0]->AddFace(face);
			neighbours[1]->AddFace(face);
			this->InteriorFaces.push_back(face);
		}
		else
			assert(false);
		this->Faces.push_back(face);
		v1->Faces.push_back(face);
		v2->Faces.push_back(face);
		v3->Faces.push_back(face);
	}
}