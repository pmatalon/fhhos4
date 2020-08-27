#pragma once
#include <gmsh.h>
#include "PolyhedralMesh.h"
#include "2D/QuadrilateralElement.h"
#include "3D/TetrahedralElement.h"
#include "2D/Square_TriangularMesh.h"
using namespace std;

enum GMSHElementTypes
{
	GMSH_Segment = 1,
	GMSH_Triangle = 2,
	GMSH_Quadrilateral = 3,
	GMSH_Tetrahedron = 4,
	GMSH_Hexahedron = 5,
	GMSH_Prism = 6,
	GMSH_Pyramid = 7
};

enum GMSHFaceTypes
{
	GMSH_TriangleFace = 3,
	GMSH_QuadrilateralFace = 4
};

template <int Dim>
class GMSHMesh : virtual public PolyhedralMesh<Dim>
{
protected:
	string _gmshFilePath;
	string _description = "GMSH file";
	string _fileNamePart = "gmsh-file";
	string _geometryDescription = "unknown";

	map<size_t, Element<Dim>*> _elementExternalNumbers;
	map<size_t, MeshVertex<Dim>*> _vertexExternalNumbers;
	double _h = -1;
	double _regularity = 1;
	BigNumber _N;
public:
	GMSHMesh(string mshFile, string description, string fileNamePart, BigNumber n = 0) 
		: PolyhedralMesh<Dim>()
	{
		_description = description;
		_fileNamePart = fileNamePart;
		_geometryDescription = Utils::FileNameWithoutExtension(mshFile);

		_N = n;

		mshFile = SearchFile(mshFile);
		_gmshFilePath = mshFile;

		bool tmpFileCreated = false;
		if (n != 0)
		{
			mshFile = SetN(mshFile, n);
			tmpFileCreated = true;
		}

		gmsh::initialize();
		//gmsh::option::setNumber("General.Terminal", 1);
		//gmsh::option::setNumber("General.Verbosity", 99);

		if (n == 0)
			cout << "Opening file " << mshFile << endl;

		gmsh::open(mshFile);

		Build();

		if (tmpFileCreated)
			remove(mshFile.c_str());
	}

	GMSHMesh(string mshFile, string description, string fileNamePart, string geometryDescription, BigNumber N = 0) :
		GMSHMesh(mshFile, description, fileNamePart, N)
	{
		_geometryDescription = geometryDescription;
	}

	GMSHMesh(string mshFile, BigNumber N = 0) : 
		GMSHMesh(mshFile, "GMSH file", "gmsh-file", N)
	{}
protected:
	GMSHMesh(string description, string fileNamePart) : PolyhedralMesh<Dim>()
	{
		_description = description;
		_fileNamePart = fileNamePart;
		Build();
	}

	// Dim-specific functions
	virtual Element<Dim>* CreateElement(int elemType, size_t elementTag, const vector<size_t>& elementNodes, size_t& elemNodeIndex, BigNumber elemNumber) { return nullptr; }
	void CreateFaces(int elemType, BigNumber& faceNumber) { }
	virtual Face<Dim>* GetBoundaryFaceFromGMSHNodes(int faceType, const vector<size_t>& faceNodes, size_t& faceNodeIndex) { return nullptr; }

public:
	BigNumber N()
	{
		return _N > 0 ? _N : (1/_h - 1);
	}

	static int GetDimension(string mshFile)
	{
		mshFile = SearchFile(mshFile);
		gmsh::initialize();
		gmsh::open(mshFile);
		int d = gmsh::model::getDimension();
		gmsh::finalize();
		return d;
	}

private:
	static string SearchFile(string mshFile)
	{
		if (!Utils::FileExists(mshFile))
		{
			if (Utils::FileExists(Mesh<Dim>::MeshDirectory + mshFile))
				mshFile = Mesh<Dim>::MeshDirectory + mshFile;
			else if (Utils::FileExists(Mesh<Dim>::MeshDirectory + "2D/" + mshFile))
				mshFile = Mesh<Dim>::MeshDirectory + "2D/" + mshFile;
			else if (Utils::FileExists(Mesh<Dim>::MeshDirectory + "3D/" + mshFile))
				mshFile = Mesh<Dim>::MeshDirectory + "3D/" + mshFile;
			else
			{
				Utils::Error("File not found: " + mshFile);
				Utils::Error("File not found: " + Mesh<Dim>::MeshDirectory + mshFile);
				Utils::FatalError("File not found: " + Mesh<Dim>::MeshDirectory + "2D/" + mshFile);
				Utils::FatalError("File not found: " + Mesh<Dim>::MeshDirectory + "3D/" + mshFile);
			}
		}
		return mshFile;
	}

	static string SetN(string geoFile, BigNumber n)
	{
		assert(n > 1);
		string geoFileWithN = Utils::Directory(geoFile) + "/" + Utils::FileNameWithoutExtension(geoFile) + "_tmp_set_N.geo";
		ofstream gmshScriptWithN(geoFileWithN);

		cout << "Opening file " << geoFile << endl;
		ifstream gmshScript(geoFile);

		cout << "Setting N = " << n << endl;
		if (!gmshScriptWithN.is_open())
		{
			gmshScript.close();
			Utils::FatalError("Unable to create temporary file " + geoFileWithN);
		}

		bool NFound = false;
		string line;
		while (getline(gmshScript, line))
		{
			if (line.find("N =") == 0 || line.find("N=") == 0)
			{
				NFound = true;
				gmshScriptWithN << "N = " << n << ";\r" << endl;
			}
			else
				gmshScriptWithN << line << endl;
		}

		gmshScript.close();
		gmshScriptWithN.close();

		if (!NFound)
		{
			remove(geoFileWithN.c_str());
			Utils::FatalError("No variable N declared in the file " + geoFile);
		}

		return geoFileWithN;
	}

	void Build()
	{
		//-----------------//
		// Physical groups //
		//-----------------//

		vector<PhysicalGroup> physicalGroups = GetPhysicalGroups(Dim);
		vector<PhysicalGroup> boundaries = GetPhysicalGroups(Dim-1);
		
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
			gmsh::model::mesh::generate(Dim); // mesh generation
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
				_vertexExternalNumbers.insert({ nodeTags[i], v });
				this->Vertices.push_back(v);
				vertexNumber++;
			}
		}

		//----------//
		// Elements //
		//----------//

		vector<int> elementTypes;
		BigNumber elemNumber = 0;

		// Get entities (surfaces in 2D, volumes in 3D). For simple geometries, entities corresponds to physical groups.
		gmsh::vectorpair entitiesDimTags;
		gmsh::model::getEntities(entitiesDimTags, Dim);
		for (int k = 0; k < entitiesDimTags.size(); k++)
		{
			int entityTag = entitiesDimTags[k].second;

			vector<int> physicalTags;
			gmsh::model::getPhysicalGroupsForEntity(Dim, entityTag, physicalTags);
			if (!physicalGroups.empty() && physicalTags.empty())
				Utils::FatalError("Entity " + to_string(entityTag) + " has no physical group. Check GMSH file.");
			if (physicalGroups.empty() && !physicalTags.empty())
				Utils::FatalError("Entity " + to_string(entityTag) + " has a physical group although no physical group has been defined. This should never happen.");
			if (physicalTags.size() > 1)
				Utils::FatalError("Entity " + to_string(entityTag) + " must have only one physical group (" + to_string(physicalTags.size()) + " found). Check GMSH file.");
			
			int physicalGroupId = physicalGroups.empty() ? 0 : physicalTags[0];


			// Get elements in this entity
			elementTypes.clear();
			vector<vector<size_t>> elementTags;
			vector<vector<size_t>> elemNodeTags;
			gmsh::model::mesh::getElements(elementTypes, elementTags, elemNodeTags, Dim, entityTag);

			for (int i = 0; i < elementTypes.size(); i++)
			{
				int elemType = elementTypes[i];
				vector<size_t> elements = elementTags[i];
				vector<size_t> elementNodes = elemNodeTags[i];

				size_t elemNodeIndex = 0;
				for (size_t j = 0; j < elements.size(); j++)
				{
					Element<Dim>* e = CreateElement(elemType, elements[j], elementNodes, elemNodeIndex, elemNumber);

					e->PhysicalGroupId = physicalGroupId;

					for (Vertex* v : e->Vertices())
					{
						MeshVertex<Dim>* mv = (MeshVertex<Dim>*)v;
						mv->Elements.push_back(e);
					}

					this->Elements.push_back(e);
					elemNumber++;

					if (e->Diameter() > this->_h)
						this->_h = e->Diameter();

					if (e->Regularity() < this->_regularity)
						this->_regularity = e->Regularity();
				}
			}
		}

		//-----------//
		//   Faces   //
		//-----------//

		// Create all the faces
		BigNumber faceNumber = 0;
		for (int i = 0; i < elementTypes.size(); i++)
		{
			int elemType = elementTypes[i];
			CreateFaces(elemType, faceNumber);
		}

		// Affectation of the physical boundaries
		if (!boundaries.empty())
		{
			gmsh::vectorpair faceEntitiesDimTags;
			gmsh::model::getEntities(faceEntitiesDimTags, Dim - 1);

			vector<int> faceTypes;
			int nFaces = 0;
			for (int k = 0; k < faceEntitiesDimTags.size(); k++)
			{
				int entityTag = faceEntitiesDimTags[k].second;

				vector<int> physicalTags;
				gmsh::model::getPhysicalGroupsForEntity(Dim - 1, entityTag, physicalTags);

				if (physicalTags.empty())
					continue;
				if (physicalTags.size() > 1)
					Utils::FatalError("Entity " + to_string(entityTag) + " must have only one physical group (" + to_string(physicalTags.size()) + " found). Check GMSH file.");

				int physicalBoundaryId = physicalTags[0];

				// Get faces in this entity
				faceTypes.clear();
				vector<vector<size_t>> faceTags;
				vector<vector<size_t>> faceNodeTags;
				gmsh::model::mesh::getElements(faceTypes, faceTags, faceNodeTags, Dim - 1, entityTag);

				for (int i = 0; i < faceTypes.size(); i++)
				{
					int faceType = faceTypes[i];
					vector<size_t> faces = faceTags[i];
					vector<size_t> faceNodes = faceNodeTags[i];

					size_t faceNodeIndex = 0;
					for (size_t j = 0; j < faces.size(); j++)
					{
						Face<Dim>* f = GetBoundaryFaceFromGMSHNodes(faceType, faceNodes, faceNodeIndex);
						f->PhysicalBoundaryId = physicalBoundaryId;
					}
				}
			}
		}

	}

protected:
	inline MeshVertex<Dim>* GetVertexFromGMSHTag(BigNumber nodeTag)
	{
		return _vertexExternalNumbers.at(nodeTag);
	}

	inline Element<Dim>* GetElementFromGMSHTag(BigNumber elementTag)
	{
		return _elementExternalNumbers.at(elementTag);
	}

	inline Face<Dim>* GetFaceFromGMSHTag(BigNumber faceTag)
	{
		return _faceExternalNumbers.at(faceTag);
	}

	static vector<PhysicalGroup> GetPhysicalGroups(int dim)
	{
		gmsh::vectorpair phyGroupsDimTags;
		gmsh::model::getPhysicalGroups(phyGroupsDimTags, dim);
		vector<PhysicalGroup> physicalGroups;
		for (size_t i = 0; i < phyGroupsDimTags.size(); i++)
		{
			PhysicalGroup phyGroup(phyGroupsDimTags[i].second);
			gmsh::model::getPhysicalName(dim, phyGroup.Id, phyGroup.Name);
			physicalGroups.push_back(phyGroup);
		}
		return physicalGroups;
	}

public:
	virtual ~GMSHMesh()
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
	string GeometryDescription() override
	{
		return _geometryDescription;
	}

	double H() override
	{
		return _h;
	}

	double Regularity() override
	{
		return _regularity;
	}

	void CoarsenMesh(CoarseningStrategy strategy) override
	{
		if (strategy == CoarseningStrategy::GMSHSplittingRefinement || strategy == CoarseningStrategy::BeyRefinement)
			return;
		else if (strategy == CoarseningStrategy::IndependentRemeshing)
			IndependentRemesh();
		else
			PolyhedralMesh<Dim>::CoarsenMesh(strategy);

		this->CoarseMesh->SetDiffusionField(this->_diffusionField);
		this->CoarseMesh->SetBoundaryConditions(this->_boundaryConditions);
	}

	virtual void RefineMesh(CoarseningStrategy strategy)
	{
		if (this->FineMesh)
			assert(false && "Mesh already refined!");

		if (strategy == CoarseningStrategy::GMSHSplittingRefinement)
			RefineMeshBySplitting();
		else
			PolyhedralMesh<Dim>::RefineMesh(strategy);
	}

	virtual void RefineMeshBySplitting()
	{
		cout << "Mesh refinement by splitting" << endl;

		string coarse_mesh_tmp_file = "./temporary_coarse.msh";
		string fine_mesh_tmp_file = "./temporary_fine.msh";

		// Save the current mesh in a temporary file, because the refinement will delete it
		gmsh::write(coarse_mesh_tmp_file);

		// Mesh refinement by splitting
		gmsh::model::mesh::refine();

		// Building our own mesh objects from the GMSH ones
		GMSHMesh<Dim>* fineMesh = CreateNewGMSHMesh();
		fineMesh->ComesFrom.CS = CoarseningStrategy::GMSHSplittingRefinement;

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
		fineMesh->LinkFacesToCoarseFaces();
	}

	virtual GMSHMesh<Dim>* CreateNewGMSHMesh()
	{
		GMSHMesh<Dim>* mesh = new GMSHMesh<Dim>(this->_description, this->_fileNamePart);
		mesh->_geometryDescription = this->_geometryDescription;
		return mesh;
	}

private:
	void IndependentRemesh()
	{
		cout << "Coarsening by independent remeshing" << endl;

		GMSHMesh<Dim>* fineMesh = this;

		GMSHMesh<Dim>* coarseMesh = new GMSHMesh<Dim>(_gmshFilePath, _description, _fileNamePart, fineMesh->N() / 2);
		coarseMesh->ComesFrom.CS = CoarseningStrategy::IndependentRemeshing;

		fineMesh->CoarseMesh = coarseMesh;
		coarseMesh->FineMesh = fineMesh;

		for (Element<Dim>* fine : fineMesh->Elements)
		{
			Element<Dim>* coarse = coarseMesh->LocateElementThatEmbeds(fine);
			coarse->FinerElements.push_back(fine);
			fine->CoarserElement = coarse;
		}
	}


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

		return GetElementFromGMSHTag(coarseElementTag);
	}

protected:
	void RenumberLike(Mesh<Dim>* myMesh)
	{
		// Renumbering the elements like myMesh
		assert(this->Elements.size() == myMesh->Elements.size());
		for (auto gmshElem : this->Elements)
		{
			bool matchFound = false;
			for (auto myElem : myMesh->Elements)
			{
				if (myElem->HasSameVertices(gmshElem, true))
				{
					gmshElem->Number = myElem->Number;
					matchFound = true;
					break;
				}
			}
			assert(matchFound);
		}

		// Renumbering the faces like myMesh
		assert(this->Faces.size() == myMesh->Faces.size());
		for (auto gmshFace : this->Faces)
		{
			bool matchFound = false;
			for (auto myFace : myMesh->Faces)
			{
				if (myFace->HasSameVertices(gmshFace, true))
				{
					gmshFace->Number = myFace->Number;
					matchFound = true;
					break;
				}
			}
			assert(matchFound);
		}

		// Reordering the faces in the lists so that we still get the same numbers as myMesh when HHO renumbers the faces
		// in the order it finds them on the lists
		/*vector<Face<Dim>*> facesCopy = this->Faces;
		this->Faces.clear();
		for (auto myFace : myMesh->Faces)
		{
			for (auto gmshFace : facesCopy)
			{
				if (gmshFace->Number == myFace->Number)
				{
					this->Faces.push_back(gmshFace);
					break;
				}
			}
		}
		assert(facesCopy.size() == this->Faces.size());*/ // not mandatory to retrieve performance

		Reorder(this->InteriorFaces, myMesh->InteriorFaces);
		//std::srand(unsigned(std::time(0)));
		//std::random_shuffle(this->InteriorFaces.begin(), this->InteriorFaces.end());

		//Reorder(this->BoundaryFaces, myMesh->BoundaryFaces); // not mandatory to retrieve performance

		if (this->CoarseMesh)
		{
			myMesh->CoarsenMesh(CoarseningStrategy::StandardCoarsening);
			dynamic_cast<GMSHMesh<Dim>*>(this->CoarseMesh)->RenumberLike(myMesh->CoarseMesh);
		}
	}

	void Reorder(vector<Face<Dim>*>& listToReorder, vector<Face<Dim>*> orderedList)
	{
		assert(listToReorder.size() == orderedList.size());
		listToReorder.clear();
		for (auto myFace : orderedList)
		{
			for (auto gmshFace : this->Faces)
			{
				if (gmshFace->Number == myFace->Number)
				{
					listToReorder.push_back(gmshFace);
					break;
				}
			}
		}
	}

public:
	virtual void RenumberLikeMe()
	{
		if (this->_fileNamePart.compare("gmsh-tri") == 0)
		{
			BigNumber n = sqrt(this->Elements.size() / 2);
			Square_TriangularMesh* myMesh = new Square_TriangularMesh(n, n);

			RenumberLike(dynamic_cast<Mesh<Dim>*>(myMesh));
			delete myMesh;
		}
		else
			assert(false);
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

		e = new QuadrilateralElement(elemNumber, GetVertexFromGMSHTag(nodeTag1), GetVertexFromGMSHTag(nodeTag2), GetVertexFromGMSHTag(nodeTag3), GetVertexFromGMSHTag(nodeTag4));
		_elementExternalNumbers.insert({ elementTag, e });
		elemNodeIndex += 4;
	}
	else if (elemType == GMSH_Triangle)
	{
		size_t nodeTag1 = elementNodes[elemNodeIndex];
		size_t nodeTag2 = elementNodes[elemNodeIndex + 1];
		size_t nodeTag3 = elementNodes[elemNodeIndex + 2];

		e = new TriangularElement(elemNumber, GetVertexFromGMSHTag(nodeTag1), GetVertexFromGMSHTag(nodeTag2), GetVertexFromGMSHTag(nodeTag3));
		_elementExternalNumbers.insert({ elementTag, e });
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

		e = new TetrahedralElement(elemNumber, GetVertexFromGMSHTag(nodeTag1), GetVertexFromGMSHTag(nodeTag2), GetVertexFromGMSHTag(nodeTag3), GetVertexFromGMSHTag(nodeTag4));
		_elementExternalNumbers.insert({ elementTag, e });
		elemNodeIndex += 4;
	}
	else if (elemType == GMSH_Hexahedron && this->FileNamePart().compare("gmsh-cart") == 0)
	{
		Vertex* v1 = GetVertexFromGMSHTag(elementNodes[elemNodeIndex]);     // (0, 1, 1)
		Vertex* v2 = GetVertexFromGMSHTag(elementNodes[elemNodeIndex + 1]); // (0, 1, 0)
		Vertex* v3 = GetVertexFromGMSHTag(elementNodes[elemNodeIndex + 2]); // (1, 1, 0)
		Vertex* v4 = GetVertexFromGMSHTag(elementNodes[elemNodeIndex + 3]); // (1, 1, 1)
		Vertex* v5 = GetVertexFromGMSHTag(elementNodes[elemNodeIndex + 4]); // (0, 0, 1)
		Vertex* v6 = GetVertexFromGMSHTag(elementNodes[elemNodeIndex + 5]); // (0, 0, 0)
		Vertex* v7 = GetVertexFromGMSHTag(elementNodes[elemNodeIndex + 6]); // (1, 0, 0)
		Vertex* v8 = GetVertexFromGMSHTag(elementNodes[elemNodeIndex + 7]); // (1, 0, 1)

		Vertex* backLeftBottomCorner = v6;
		Vertex* frontLeftBottomCorner = v7;
		Vertex* backRightBottomCorner = v2;
		Vertex* backLeftTopCorner = v5;
		Vertex* frontLeftTopCorner = v8;
		Vertex* backRightTopCorner = v1;
		Vertex* frontRightBottomCorner = v3;
		Vertex* frontRightTopCorner = v4;

		e = new ParallelepipedElement(elemNumber, backLeftBottomCorner, frontLeftBottomCorner, backRightBottomCorner, backLeftTopCorner, frontLeftTopCorner, backRightTopCorner, frontRightBottomCorner, frontRightTopCorner);
		_elementExternalNumbers.insert({ elementTag, e });
		elemNodeIndex += 8;
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
		MeshVertex<2>* v1 = (MeshVertex<2>*)GetVertexFromGMSHTag(edgeNodes[j]);
		MeshVertex<2>* v2 = (MeshVertex<2>*)GetVertexFromGMSHTag(edgeNodes[j + 1]);

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
	int faceType = -1;
	int nFaceVertices = -1;
	if (elemType == GMSHElementTypes::GMSH_Tetrahedron)
	{
		faceType = GMSHFaceTypes::GMSH_TriangleFace;
		nFaceVertices = 3;
	}
	else if (elemType == GMSHElementTypes::GMSH_Hexahedron)
	{
		faceType = GMSHFaceTypes::GMSH_QuadrilateralFace;
		nFaceVertices = 4;
	}
	else
		assert(false);

	vector<size_t> faceNodes;
	bool onlyPrimaryNodes = true;
	gmsh::model::mesh::getElementFaceNodes(elemType, faceType, faceNodes, -1, onlyPrimaryNodes);

	for (size_t j = 0; j < faceNodes.size(); j += nFaceVertices)
	{
		vector<MeshVertex<3>*> vertices(nFaceVertices);
		for (int k = 0; k < nFaceVertices; k++)
		{
			MeshVertex<3>* v = (MeshVertex<3>*)GetVertexFromGMSHTag(faceNodes[j+k]);
			vertices[k] = v;
		}

		// Check if a face with those vertices already exists
		bool faceAlreadyExists = false;
		for (Face<3>* f : vertices[0]->Faces)
		{
			bool thisFaceHasAllVertices = true;
			for (int k = 1; k < nFaceVertices; k++)
			{
				if (!f->HasVertex(vertices[k]))
				{
					thisFaceHasAllVertices = false;
					break;
				}
			}
			if (thisFaceHasAllVertices)
			{
				faceAlreadyExists = true;
				break;
			}
		}
		if (faceAlreadyExists)
			continue;

		// Find the neighbouring elements this future face is the interface of
		vector<Element<3>*> neighbours;
		for (Element<3>* e : vertices[0]->Elements)
		{
			bool eHasAllVertices = true;
			for (int k = 1; k < nFaceVertices; k++)
			{
				bool eHasVertexk = false;
				for (Element<3>* ek : vertices[k]->Elements)
				{
					if (ek == e)
					{
						eHasVertexk = true;
						break;
					}
				}
				if (!eHasVertexk)
				{
					eHasAllVertices = false;
					break;
				}
			}
			if (eHasAllVertices)
				neighbours.push_back(e);
		}


		// Creation of the face
		Face<3>* face;
		if (faceType == GMSHFaceTypes::GMSH_TriangleFace)
		{
			if (neighbours.size() == 1)
				face = new TriangularFace(faceNumber++, vertices[0], vertices[1], vertices[2], neighbours[0]);
			else if (neighbours.size() == 2)
				face = new TriangularFace(faceNumber++, vertices[0], vertices[1], vertices[2], neighbours[0], neighbours[1]);
		}
		else if (faceType == GMSHFaceTypes::GMSH_QuadrilateralFace && this->FileNamePart().compare("gmsh-cart") == 0)
		{
			Vertex* v1 = vertices[0]; // (0, 1, 1)
			Vertex* v2 = vertices[1]; // (1, 1, 1)
			Vertex* v3 = vertices[2]; // (1, 1, 0)
			Vertex* v4 = vertices[3]; // (0, 1, 0)

			CartesianShapeOrientation orientation = CartesianShapeOrientation::None;
			Vertex* origin;
			Vertex* vertex1;
			Vertex* vertex2;

			DimVector<3> v12 = Vect<3>(v1, v2);
			DimVector<3> v13 = Vect<3>(v1, v3);
			//DimVector<3> normal = Vect<3>(v2, v1).cross(Vect<3>(v3, v1));
			DimVector<3> unitX; unitX << 1, 0, 0;
			DimVector<3> unitY; unitY << 0, 1, 0;
			DimVector<3> unitZ; unitZ << 0, 0, 1;
			if (abs(v12.dot(unitX)) < 1e-14 && abs(v13.dot(unitX)) < 1e-14)
			{
				orientation = CartesianShapeOrientation::InYOZ;
			}
			else if (abs(v12.dot(unitY)) < 1e-14 && abs(v13.dot(unitY)) < 1e-14)
			{
				orientation = CartesianShapeOrientation::InXOZ;
				origin = v4;
				vertex1 = v3;
				vertex2 = v1;
			}
			else if (abs(v12.dot(unitZ)) < 1e-14 && abs(v13.dot(unitZ)) < 1e-14)
			{
				orientation = CartesianShapeOrientation::InXOY;
			}
			else
				assert(false);

			DimVector<3> OV1; OV1 << v1->X, v1->Y, v1->Z;
			DimVector<3> OV2; OV2 << v2->X, v2->Y, v2->Z;
			DimVector<3> OV3; OV3 << v3->X, v3->Y, v3->Z;
			DimVector<3> OV4; OV4 << v4->X, v4->Y, v4->Z;

			double minNorm = min({OV1.norm(), OV2.norm(), OV3.norm(), OV4.norm()});
			if (minNorm == OV1.norm())
			{
				origin = v1;
				vertex1 = v4;
				vertex2 = v2;
			}
			else if (minNorm == OV2.norm())
			{
				origin = v2;
				vertex1 = v3;
				vertex2 = v1;
			}
			else if (minNorm == OV3.norm())
			{
				origin = v3;
				vertex1 = v2;
				vertex2 = v4;
			}
			else if (minNorm == OV4.norm())
			{
				origin = v4;
				vertex1 = v1;
				vertex2 = v3;
			}
			else
				assert(false);

			if (neighbours.size() == 1)
				face = new RectangularFace(faceNumber++, origin, vertex1, vertex2, neighbours[0], orientation);
			else if (neighbours.size() == 2)
				face = new RectangularFace(faceNumber++, origin, vertex1, vertex2, neighbours[0], neighbours[1], orientation);
		}
		else
			assert(false);

		
		if (neighbours.size() == 1)
		{
			neighbours[0]->AddFace(face);
			this->BoundaryFaces.push_back(face);
		}
		else if (neighbours.size() == 2)
		{
			neighbours[0]->AddFace(face);
			neighbours[1]->AddFace(face);
			this->InteriorFaces.push_back(face);
		}
		else
			assert(false);

		this->Faces.push_back(face);
		for (MeshVertex<3>* v : vertices)
			v->Faces.push_back(face);
	}
}




template<>
Face<2>* GMSHMesh<2>::GetBoundaryFaceFromGMSHNodes(int faceType, const vector<size_t>& faceNodes, size_t& faceNodeIndex)
{
	assert(faceType == GMSHElementTypes::GMSH_Segment);

	MeshVertex<2>* v1 = (MeshVertex<2>*)GetVertexFromGMSHTag(faceNodes[faceNodeIndex]);
	MeshVertex<2>* v2 = (MeshVertex<2>*)GetVertexFromGMSHTag(faceNodes[faceNodeIndex + 1]);
	faceNodeIndex += 2;

	for (Face<2>* f : v1->Faces)
	{
		if (f->IsDomainBoundary && f->HasVertex(v2))
			return f;
	}
	assert(false && "Face not found");
}