#pragma once
#include <gmsh.h>
#include "../../Mesh/PolyhedralMesh.h"
#include "../InHouse/Square_TriangularMesh.h"
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
	GMSH_SegmentFace,
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
private:
	string _mshFilePath;
	bool _mshFileIsTmp = false;
public:
	static bool GMSHLogEnabled;
	static bool UseCache;

	GMSHMesh(string gmshFile, string description, string fileNamePart, BigNumber n = 0, bool buildMesh = true) 
		: PolyhedralMesh<Dim>()
	{
		_description = description;
		_fileNamePart = fileNamePart;
		_geometryDescription = FileSystem::FileNameWithoutExtension(gmshFile);

		_N = n;

		gmshFile = SearchFile(gmshFile);
		_gmshFilePath = gmshFile;

		// Initialize GMSH
		gmsh::initialize();

		// Enable the logs from GMSH to be printed in the console
		if (GMSHLogEnabled)
		{
			gmsh::option::setNumber("General.Terminal", 1);
			gmsh::option::setNumber("General.Verbosity", 99);
		}

		bool tmpGeoFileCreated = false;
		string mshFile;
		string extension = FileSystem::Extension(gmshFile);
		if (extension.compare("geo") == 0)
		{
			string mshFileInCache = Mesh<Dim>::CacheDirectory + _geometryDescription + (n == 0 ? "" : "_N" + to_string(n)) + ".msh";

			if (UseCache && FileSystem::FileExists(mshFileInCache))
			{
				cout << "Loading file from cache: " << mshFileInCache << endl;
				gmsh::open(mshFileInCache);
			}
			else
			{
				if (n != 0)
				{
					string tmpGeoFile = SetN(gmshFile, n);
					tmpGeoFileCreated = true;
					cout << "Loading file: " << tmpGeoFile << endl;
					gmsh::open(tmpGeoFile);
					remove(tmpGeoFile.c_str());
				}
				else
				{
					cout << "Loading file: " << gmshFile << endl;
					gmsh::open(gmshFile);
				}

				// Generate the mesh
				cout << "Generating the mesh by GMSH..." << endl;
				gmsh::model::mesh::generate(Dim);

				if (UseCache)
				{
					// Write the mesh in cache
					cout << "Saving the mesh in cache" << endl;
					FileSystem::CreateDirectoryIfNotExist(Mesh<Dim>::CacheDirectory);
					gmsh::write(mshFileInCache);
					if (!FileSystem::FileExists(mshFileInCache))
						Utils::Error("Failed to write the file " + mshFileInCache);
				}
			}

			if (UseCache)
				_mshFilePath = mshFileInCache;
		}
		else if (extension.compare("msh") == 0)
		{
			cout << "Loading file: " << gmshFile << endl;
			gmsh::open(gmshFile);
		}
		else
			Utils::FatalError("Unknown file extension: '" + extension + "'. Only .geo and .msh allowed.");

		if (buildMesh)
			Build();
	}

	GMSHMesh(string geoFile, string description, string fileNamePart, string geometryDescription, BigNumber N = 0) :
		GMSHMesh(geoFile, description, fileNamePart, N)
	{
		_geometryDescription = geometryDescription;
	}

	GMSHMesh(string geoFile, BigNumber N = 0) :
		GMSHMesh(geoFile, "GMSH file", "gmsh-file", N)
	{}
protected:
	GMSHMesh(string description, string fileNamePart) : PolyhedralMesh<Dim>()
	{
		_description = description;
		_fileNamePart = fileNamePart;
	}

	// Dim-specific functions
	virtual Element<Dim>* CreateElement(int elemType, const vector<size_t>& elementNodes, size_t start, size_t elemIndex) { return nullptr; }
	void CreateFaces(int elemType, BigNumber& faceNumber) { }
	virtual Face<Dim>* GetBoundaryFaceFromGMSHNodes(int faceType, const vector<size_t>& faceNodes, size_t& faceNodeIndex) { return nullptr; }

public:
	BigNumber N()
	{
		return _N > 0 ? _N : (1/_h - 1);
	}

	static int GetDimension(string geoFile)
	{
		geoFile = SearchFile(geoFile);
		gmsh::initialize();
		gmsh::open(geoFile);
		int d = gmsh::model::getDimension();
		gmsh::finalize();
		return d;
	}

	static void CloseGMSH()
	{
		gmsh::finalize();
	}

private:
	static string SearchFile(string geoFile)
	{
		if (!FileSystem::HasExtension(geoFile))
			geoFile += ".geo";
		if (!FileSystem::FileExists(geoFile))
		{
			if (FileSystem::FileExists(Mesh<Dim>::MeshDirectory + geoFile))
				geoFile = Mesh<Dim>::MeshDirectory + geoFile;
			else if (FileSystem::FileExists(Mesh<Dim>::MeshDirectory + "2D/" + geoFile))
				geoFile = Mesh<Dim>::MeshDirectory + "2D/" + geoFile;
			else if (FileSystem::FileExists(Mesh<Dim>::MeshDirectory + "3D/" + geoFile))
				geoFile = Mesh<Dim>::MeshDirectory + "3D/" + geoFile;
			else
			{
				Utils::Error("File not found: " + geoFile);
				Utils::Error("File not found: " + Mesh<Dim>::MeshDirectory + geoFile);
				Utils::FatalError("File not found: " + Mesh<Dim>::MeshDirectory + "2D/" + geoFile);
				Utils::FatalError("File not found: " + Mesh<Dim>::MeshDirectory + "3D/" + geoFile);
			}
		}
		return geoFile;
	}

	static string SetN(string geoFile, BigNumber n)
	{
		assert(n > 1);
		string geoFileWithN = FileSystem::Directory(geoFile) + "/" + FileSystem::FileNameWithoutExtension(geoFile) + "_N" + to_string(n) + ".tmp.geo";
		ofstream gmshScriptWithN(geoFileWithN);

		cout << "Opening file: " << geoFile << endl;
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
			Utils::Warning("No variable N declared in the file " + geoFile + ". Using the file as is.");
		}

		return geoFileWithN;
	}

	void Build()
	{
		//-----------------//
		// Physical groups //
		//-----------------//

		if (this->PhysicalParts.empty())
			this->PhysicalParts = GetPhysicalGroups();
		if (this->BoundaryParts.empty())
			this->BoundaryParts = GetBoundaryGroups();
		
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
			Utils::FatalError("The mesh isn't readable.");

		cout << "Building internal mesh objects..." << endl;

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

		// Get entities (surfaces in 2D, volumes in 3D). For simple geometries, entities corresponds to physical groups.
		gmsh::vectorpair entitiesDimTags;
		gmsh::model::getEntities(entitiesDimTags, Dim);

		if (entitiesDimTags.size() > 1) // if = 1, then the memory allocation is done when processing elements
		{
			// Get the total number of elements to allocate the memory

			//gmsh::model::mesh::getElementTypes(elementTypes, Dim);
			vector<int> elementTypes;
			vector<vector<size_t>> elementTags;
			vector<vector<size_t>> elemNodeTags;
			gmsh::model::mesh::getElements(elementTypes, elementTags, elemNodeTags, Dim);
			for (size_t i = 0; i < elementTypes.size(); i++)
			{
				int elemType = elementTypes[i];
				vector<size_t> elements = elementTags[i];
				vector<size_t> elementNodes = elemNodeTags[i];

				AllocateContiguousMemoryForElements(elemType, elements.size());
			}
		}


		vector<int> elementTypes;
		elementTypes.clear();

		map<int, size_t> elementsAlreadyInMemory; // the key (int) is the GMSH elemType

		for (int k = 0; k < entitiesDimTags.size(); k++)
		{
			int entityTag = entitiesDimTags[k].second;

			vector<int> physicalTags;
			gmsh::model::getPhysicalGroupsForEntity(Dim, entityTag, physicalTags);
			if (!this->PhysicalParts.empty() && physicalTags.empty())
				Utils::FatalError("Entity " + to_string(entityTag) + " has no physical group. Check GMSH file.");
			if (this->PhysicalParts.empty() && !physicalTags.empty())
				Utils::FatalError("Entity " + to_string(entityTag) + " has a physical group although no physical group has been defined. This should never happen.");
			if (physicalTags.size() > 1)
				Utils::FatalError("Entity " + to_string(entityTag) + " must have only one physical group (" + to_string(physicalTags.size()) + " found). Check GMSH file.");
			
			PhysicalGroup<Dim>* physicalPart = nullptr;
			for (PhysicalGroup<Dim>* pp : this->PhysicalParts)
			{
				if (pp->Id == physicalTags[0])
				{
					physicalPart = pp;
					break;
				}
			}
			assert(physicalPart);

			// Get elements in this entity
			elementTypes.clear();
			vector<vector<size_t>> elementTags;
			vector<vector<size_t>> elemNodeTags;
			gmsh::model::mesh::getElements(elementTypes, elementTags, elemNodeTags, Dim, entityTag);

			struct ChunkResult
			{
				vector<Element<Dim>*> Elements;
				double h = -1;
				double regularity = 2;
			};

			for (size_t i = 0; i < elementTypes.size(); i++)
			{
				int elemType = elementTypes[i];
				vector<size_t> elements = elementTags[i];
				vector<size_t> elementNodes = elemNodeTags[i];

				if (entitiesDimTags.size() == 1) // if > 1, then already allocated above
					AllocateContiguousMemoryForElements(elemType, elements.size());

				if (elementsAlreadyInMemory.find(elemType) == elementsAlreadyInMemory.end())
					elementsAlreadyInMemory.insert({ elemType, 0 });
				size_t nElementsAlreadyInMemory = elementsAlreadyInMemory[elemType];

				NumberParallelLoop<ChunkResult> parallelLoop(elements.size());
				parallelLoop.Execute([this, elemType, nElementsAlreadyInMemory, &elements, &elementNodes, &physicalPart](BigNumber j, ParallelChunk<ChunkResult>* chunk)
					{
						size_t elementTag = elements[j];
						Element<Dim>* e = CreateElement(elemType, elementNodes, nElementsAlreadyInMemory, j);

						e->Id = elementTag;
						e->PhysicalPart = physicalPart;

						for (Vertex* v : e->Vertices())
						{
							MeshVertex<Dim>* mv = (MeshVertex<Dim>*)v;
							mv->Mutex.lock();
							mv->Elements.push_back(e);
							mv->Mutex.unlock();
						}

						chunk->Results.Elements.push_back(e);

						chunk->Results.h          = max(chunk->Results.h,          e->Diameter());
						chunk->Results.regularity = min(chunk->Results.regularity, e->Regularity());
					});

				parallelLoop.AggregateChunkResults([this](ChunkResult& chunk)
					{
						for (Element<Dim>* e : chunk.Elements)
						{
							this->AddElement(e, false);
							_elementExternalNumbers.insert({ e->Id, e });
						}

						_h          = max(_h,          chunk.h);
						_regularity = min(_regularity, chunk.regularity);
					});

				elementsAlreadyInMemory[elemType] += elements.size();
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
		if (!this->BoundaryParts.empty())
		{
			gmsh::vectorpair faceEntitiesDimTags;
			gmsh::model::getEntities(faceEntitiesDimTags, Dim - 1);

			vector<int> faceTypes;
			for (int k = 0; k < faceEntitiesDimTags.size(); k++)
			{
				int entityTag = faceEntitiesDimTags[k].second;

				vector<int> physicalTags;
				gmsh::model::getPhysicalGroupsForEntity(Dim - 1, entityTag, physicalTags);

				if (physicalTags.empty())
					continue;
				if (physicalTags.size() > 1)
					Utils::FatalError("Entity " + to_string(entityTag) + " must have only one physical group (" + to_string(physicalTags.size()) + " found). Check GMSH file.");

				BoundaryGroup* boundaryPart = nullptr;
				for (BoundaryGroup* bp : this->BoundaryParts)
				{
					if (bp->Id == physicalTags[0])
					{
						boundaryPart = bp;
						break;
					}
				}
				assert(boundaryPart);

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
						f->BoundaryPart = boundaryPart;
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

	vector<PhysicalGroup<Dim>*> GetPhysicalGroups()
	{
		gmsh::vectorpair phyGroupsDimTags;
		gmsh::model::getPhysicalGroups(phyGroupsDimTags, Dim);
		vector<PhysicalGroup<Dim>*> physicalGroups;
		for (size_t i = 0; i < phyGroupsDimTags.size(); i++)
		{
			PhysicalGroup<Dim>* phyGroup = new PhysicalGroup<Dim>(phyGroupsDimTags[i].second);
			gmsh::model::getPhysicalName(Dim, phyGroup->Id, phyGroup->Name);
			physicalGroups.push_back(phyGroup);
		}
		return physicalGroups;
	}

	vector<BoundaryGroup*> GetBoundaryGroups()
	{
		gmsh::vectorpair phyGroupsDimTags;
		gmsh::model::getPhysicalGroups(phyGroupsDimTags, Dim-1);
		vector<BoundaryGroup*> physicalGroups;
		for (size_t i = 0; i < phyGroupsDimTags.size(); i++)
		{
			BoundaryGroup* phyGroup = new BoundaryGroup(phyGroupsDimTags[i].second);
			gmsh::model::getPhysicalName(Dim-1, phyGroup->Id, phyGroup->Name);
			physicalGroups.push_back(phyGroup);
		}
		return physicalGroups;
	}

	void AllocateContiguousMemoryForElements(int elemType, BigNumber numberOfElements)
	{
		this->Elements.reserve(numberOfElements);

		if (elemType == GMSH_Quadrilateral)
			this->_quadrilateralElements = vector<QuadrilateralElement>(numberOfElements);
		else if (elemType == GMSH_Triangle)
			this->_triangularElements = vector<TriangularElement>(numberOfElements);
		else if (elemType == GMSH_Tetrahedron)
			this->_tetrahedralElements = vector<TetrahedralElement>(numberOfElements);
		else if (elemType == GMSH_Hexahedron)
			this->_parallelepipedElements = vector<ParallelepipedElement>(numberOfElements);
		else
			Utils::FatalError("GMSH element type not managed.");
	}

	void AllocateContiguousMemoryForFaces(GMSHFaceTypes faceType, BigNumber numberOfFaces)
	{
		this->Faces.reserve(numberOfFaces);

		if (Dim == 2)
			this->_edgeFaces.reserve(numberOfFaces);
		else if (faceType == GMSHFaceTypes::GMSH_TriangleFace)
			this->_triangularFaces.reserve(numberOfFaces);
		else
			Utils::FatalError("GMSH element type not managed.");
	}

public:
	virtual ~GMSHMesh()
	{}

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
		if (Utils::IsRefinementStrategy(strategy))
			return;
		else if (strategy == CoarseningStrategy::IndependentRemeshing)
			IndependentRemesh();
		else
			PolyhedralMesh<Dim>::CoarsenMesh(strategy);
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
		this->InitializeRefinement(fineMesh);
		fineMesh->ComesFrom.CS = CoarseningStrategy::GMSHSplittingRefinement;
		if (dynamic_cast<TetrahedralElement*>(this->Elements[0]))
		{
			fineMesh->ComesFrom.nFineElementsByCoarseElement = 8;
			fineMesh->ComesFrom.nFineFacesAddedByCoarseElement = 8;
			fineMesh->ComesFrom.nFineFacesByKeptCoarseFace = 4;
		}
		fineMesh->Build();

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

		if (fineMesh->ComesFrom.HasDetails())
		{
			bool warnings = false;
			for (Element<Dim>* coarse : this->Elements)
			{
				if (coarse->FinerElements.size() != fineMesh->ComesFrom.nFineElementsByCoarseElement)
				{
					Utils::Warning("After refinement, " + to_string(coarse->FinerElements.size()) + " fine tetrahedra have been found for the tetrahedron of Tag " + to_string(coarse->Id) + ", which is weird.");
					warnings = true;
				}
			}
			if (warnings)
			{
				string msg = "These warnings usually happen when the refined mesh has some highly distorted elements (confirm by enabling GMSH logging with the option -gmsh-log).\n";
				msg += "As a consequence, the link between fine and coarse elements might be impossible to retrieve, causing a fatal error later.\n";
				msg += "To avoid this problem, try starting the refinement process from a finer mesh by using the option -coarse-n.";
				Utils::Error(msg);
			}
		}

		// Remove coarse temporary file and reload the fine mesh
		remove(coarse_mesh_tmp_file.c_str());
		gmsh::open(fine_mesh_tmp_file);
		
		fineMesh->_mshFilePath = fine_mesh_tmp_file;
		fineMesh->_mshFileIsTmp = true;

		// Continue linking
		fineMesh->LinkFacesToCoarseFaces();

		this->FinalizeRefinement();
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
		if (this->N() == 1)
			return;

		cout << "Coarsening by independent remeshing" << endl;

		GMSHMesh<Dim>* fineMesh = this;

		GMSHMesh<Dim>* coarseMesh = new GMSHMesh<Dim>(_gmshFilePath, _description, _fileNamePart, fineMesh->N() / 2, false);
		this->InitializeCoarsening(coarseMesh);
		coarseMesh->ComesFrom.CS = CoarseningStrategy::IndependentRemeshing;
		coarseMesh->Build();

		for (Element<Dim>* fine : fineMesh->Elements)
		{
			Element<Dim>* coarse = coarseMesh->LocateElementThatEmbeds(fine);
			if (!coarse->IsInSamePhysicalPartAs(fine))
				Utils::Warning("Coarse and fine elements aren't located in the same physical part. Something is wrong...");
			coarse->FinerElements.push_back(fine);
			fine->CoarserElement = coarse;
		}

		CloseGMSH();

		/*for (Element<Dim>* coarse : coarseMesh->Elements)
		{
			if (coarse->FinerElements.empty())
			{
				Utils::Warning("This coarse element does not have any fine element.");
				for (Element<Dim>* neighbour : coarse->NeighboursInSamePhysicalPart())
				{
					for (Element<Dim>* fe : neighbour->FinerElements)
					{
						if (!PolyhedralMesh<Dim>::Intersection(coarse, fe).empty())
						{
							coarse->FinerElements.push_back(fe);
							break;
						}
					}
				}
				assert(!coarse->FinerElements.empty());
			}
		}*/

		this->FinalizeCoarsening();
	}


	Element<Dim>* LocateElementThatEmbeds(Element<Dim>* finerElement)
	{
		size_t coarseElementTag = LocateGMSHElementContaining(finerElement->Center());
		return GetElementFromGMSHTag(coarseElementTag);
	}

	size_t LocateGMSHElementContaining(const DomPoint& p)
	{
		size_t elementTag;
		int elementType;
		vector<size_t> elementNodes;
		double u, v, w;
		bool strictResearch = true;

		// If this function fails, check that a Physical Surface (if 2D) is defined in the file (if Physical Lines exist, then a Physical Surface must also exist!).
		gmsh::model::mesh::getElementByCoordinates(p.X, p.Y, p.Z, elementTag,
			elementType, elementNodes, u, v, w, // useless parameters for us
			Dim, strictResearch);

		return elementTag;
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

	void Reorder(vector<Face<Dim>*>& listToReorder, const vector<Face<Dim>*>& orderedList)
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

	void ExportSolutionToGMSH(FunctionalBasis<Dim>* basis, const Vector &solution, const string& outputFilePathPrefix) override
	{
		assert(!_mshFilePath.empty());
		gmsh::initialize();
		gmsh::open(_mshFilePath);
		if (_mshFileIsTmp)
			remove(_mshFilePath.c_str());

		int viewId = gmsh::view::add("potential");

		vector<std::string> modelNames;
		gmsh::model::list(modelNames);
		string modelName = modelNames[modelNames.size() - 1];

		vector<size_t> elementTags;
		vector<vector<double>> values;
		elementTags.reserve(this->Elements.size());
		values.reserve(this->Elements.size());
		for (Element<Dim>* e : this->Elements)
		{
			elementTags.push_back(e->Id);
			double value = e->EvalApproximateSolution(basis, solution, e->Center());
			values.push_back({ value });
		}
		gmsh::view::addModelData(viewId, 0, modelName, "ElementData", elementTags, values);

		string meshFilePath = outputFilePathPrefix + ".msh";
		string dataFilePath = outputFilePathPrefix + ".pos";

		gmsh::write(meshFilePath);
		gmsh::view::write(viewId, dataFilePath);

		cout << "Solution exported for GMSH to " << dataFilePath << endl;
		gmsh::finalize();
	}

};

//-------------//
// 2D elements //
//-------------//

template <>
Element<2>* GMSHMesh<2>::CreateElement(int elemType, const vector<size_t>& elementNodes, size_t start, size_t elemIndex)
{ 
	Element<2>* e = nullptr;
	if (elemType == GMSH_Quadrilateral)
	{
		size_t elemNodeIndex = 4 * elemIndex;
		size_t nodeTag1 = elementNodes[elemNodeIndex];
		size_t nodeTag2 = elementNodes[elemNodeIndex + 1];
		size_t nodeTag3 = elementNodes[elemNodeIndex + 2];
		size_t nodeTag4 = elementNodes[elemNodeIndex + 3];

		_quadrilateralElements[start + elemIndex] = QuadrilateralElement(0, GetVertexFromGMSHTag(nodeTag1), GetVertexFromGMSHTag(nodeTag2), GetVertexFromGMSHTag(nodeTag3), GetVertexFromGMSHTag(nodeTag4));
		e = &_quadrilateralElements[start + elemIndex];
	}
	else if (elemType == GMSH_Triangle)
	{
		size_t elemNodeIndex = 3 * elemIndex;
		size_t nodeTag1 = elementNodes[elemNodeIndex];
		size_t nodeTag2 = elementNodes[elemNodeIndex + 1];
		size_t nodeTag3 = elementNodes[elemNodeIndex + 2];

		_triangularElements[start + elemIndex] = TriangularElement(0, GetVertexFromGMSHTag(nodeTag1), GetVertexFromGMSHTag(nodeTag2), GetVertexFromGMSHTag(nodeTag3));
		e = &_triangularElements[start + elemIndex];
	}
	else
		assert(false && "GMSH element type not managed.");
	return e;
}

//-------------//
// 3D elements //
//-------------//

template <>
Element<3>* GMSHMesh<3>::CreateElement(int elemType, const vector<size_t>& elementNodes, size_t start, size_t elemIndex)
{
	Element<3>* e = nullptr;
	if (elemType == GMSH_Tetrahedron)
	{
		size_t elemNodeIndex = 4 * elemIndex;
		size_t nodeTag1 = elementNodes[elemNodeIndex];
		size_t nodeTag2 = elementNodes[elemNodeIndex + 1];
		size_t nodeTag3 = elementNodes[elemNodeIndex + 2];
		size_t nodeTag4 = elementNodes[elemNodeIndex + 3];

		_tetrahedralElements[start + elemIndex] = TetrahedralElement(0, GetVertexFromGMSHTag(nodeTag1), GetVertexFromGMSHTag(nodeTag2), GetVertexFromGMSHTag(nodeTag3), GetVertexFromGMSHTag(nodeTag4));
		e = &_tetrahedralElements[start + elemIndex];
	}
	else if (elemType == GMSH_Hexahedron && this->FileNamePart().compare("gmsh-cart") == 0)
	{
		size_t elemNodeIndex = 8 * elemIndex;
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

		_parallelepipedElements[start + elemIndex] = ParallelepipedElement(0, backLeftBottomCorner, frontLeftBottomCorner, backRightBottomCorner, backLeftTopCorner, frontLeftTopCorner, backRightTopCorner, frontRightBottomCorner, frontRightTopCorner);
		e = &_parallelepipedElements[start + elemIndex];
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

	AllocateContiguousMemoryForFaces(GMSHFaceTypes::GMSH_SegmentFace, edgeNodes.size() / 2);

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
			this->_edgeFaces.emplace_back(faceNumber++, v1, v2, neighbours[0]);
			edge = &this->_edgeFaces.back();
			neighbours[0]->AddFace(edge);
			this->BoundaryFaces.push_back(edge);
		}
		else if (neighbours.size() == 2)
		{
			this->_edgeFaces.emplace_back(faceNumber++, v1, v2, neighbours[0], neighbours[1]);
			edge = &this->_edgeFaces.back();
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
	GMSHFaceTypes faceType;
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

	AllocateContiguousMemoryForFaces(faceType, faceNodes.size() / nFaceVertices);

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
				if (!vertices[k]->IsVertexOf(e))
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
			{
				this->_triangularFaces.emplace_back(faceNumber++, vertices[0], vertices[1], vertices[2], neighbours[0]);
				face = &this->_triangularFaces.back();
			}
			else if (neighbours.size() == 2)
			{
				this->_triangularFaces.emplace_back(faceNumber++, vertices[0], vertices[1], vertices[2], neighbours[0], neighbours[1]);
				face = &this->_triangularFaces.back();
			}
			else
				assert(false);
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

template<>
Face<3>* GMSHMesh<3>::GetBoundaryFaceFromGMSHNodes(int faceType, const vector<size_t>& faceNodes, size_t& faceNodeIndex)
{
	assert(faceType == GMSHElementTypes::GMSH_Triangle);

	MeshVertex<3>* v1 = (MeshVertex<3>*)GetVertexFromGMSHTag(faceNodes[faceNodeIndex]);
	MeshVertex<3>* v2 = (MeshVertex<3>*)GetVertexFromGMSHTag(faceNodes[faceNodeIndex + 1]);
	MeshVertex<3>* v3 = (MeshVertex<3>*)GetVertexFromGMSHTag(faceNodes[faceNodeIndex + 2]);
	faceNodeIndex += 3;

	for (Face<3>* f : v1->Faces)
	{
		if (f->IsDomainBoundary && f->HasVertex(v2) && f->HasVertex(v3))
			return f;
	}
	assert(false && "Face not found");
}

template <>
bool GMSHMesh<1>::GMSHLogEnabled = false;
template <>
bool GMSHMesh<2>::GMSHLogEnabled = false;
template <>
bool GMSHMesh<3>::GMSHLogEnabled = false;

template <>
bool GMSHMesh<1>::UseCache = true;
template <>
bool GMSHMesh<2>::UseCache = true;
template <>
bool GMSHMesh<3>::UseCache = true;