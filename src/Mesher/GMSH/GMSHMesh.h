#pragma once
#include <gmsh.h>
#include "../../Mesh/PolyhedralMesh.h"
#ifdef ENABLE_2D
	#include "../InHouse/Square_TriangularMesh.h"
#endif
#ifdef ENABLE_3D
	#include "../../Mesh/3D/RectangularFace.h"
#endif
#include "../../FunctionalBasis/Monomials/MonomialBasis.h"
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
	TestCase<Dim>* _testCase = nullptr;
	string _gmshFilePath;
	string _description = "GMSH file";
	string _fileNamePart = "";
	string _geometryDescription = "unknown";

	map<size_t, Element<Dim>*> _elementExternalNumbers;
	map<size_t, MeshVertex<Dim>*> _vertexExternalNumbers;
	double _h = -1;
	double _average_h = -1;
	double _regularity = 1;
	BigNumber _N;
private:
	string _mshFilePath;
	bool _mshFileIsTmp = false;
public:
	static bool GMSHLogEnabled;
	static bool UseCache;

	GMSHMesh(TestCase<Dim>* testCase, string gmshFile, string description, string fileNamePart, BigNumber n = 0, bool buildMesh = true)
		: PolyhedralMesh<Dim>()
	{
		_testCase = testCase;
		_description = description;
		_fileNamePart = fileNamePart;
		_geometryDescription = FileSystem::FileNameWithoutExtension(gmshFile);

		_N = n;

		gmshFile = SearchFile(gmshFile);
		_gmshFilePath = gmshFile;

		// Initialize GMSH
		gmsh::initialize();

		ManageGMSHLog();

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

				string optimizationMethod;
				gmsh::model::mesh::optimize(optimizationMethod);

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

	GMSHMesh(TestCase<Dim>* testCase, string geoFile, string description, string fileNamePart, string geometryDescription, BigNumber N = 0) :
		GMSHMesh(testCase, geoFile, description, fileNamePart, N)
	{
		_geometryDescription = geometryDescription;
	}

	GMSHMesh(TestCase<Dim>* testCase, string geoFile, BigNumber N = 0) :
		GMSHMesh(testCase, geoFile, "GMSH file", "", N)
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
		assert(n > 0);
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

	static void ManageGMSHLog()
	{
		// Enable the logs from GMSH to be printed in the console
		if (GMSHLogEnabled)
		{
			gmsh::option::setNumber("General.Terminal", 1);
			gmsh::option::setNumber("General.Verbosity", 99);
		}
		else
			gmsh::option::setNumber("General.Terminal", 0);
	}

	void Build()
	{
		//-----------------//
		// Physical groups //
		//-----------------//

		const int DefaultPhyGroup = -1;

		if (this->PhysicalParts.empty())
		{
			this->PhysicalParts = GetPhysicalGroups();
			if (this->PhysicalParts.empty())
			{
				cout << "Default physical part created." << endl;
				PhysicalGroup<Dim>* defaultPhyPart = new PhysicalGroup<Dim>(DefaultPhyGroup, "domain");
				this->PhysicalParts.push_back(defaultPhyPart);
			}
		}

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

		//--------------------//
		// Geometric vertices //
		//--------------------//

		// Get vertices defining the geometry.
		// These vertices must be conserved at every level of coarsening.

		if (this->ComesFrom.CS == H_CoarsStgy::IndependentRemeshing)
			this->_geometricVertices.clear();

		gmsh::vectorpair entitiesDim0Tags;
		gmsh::model::getEntities(entitiesDim0Tags, 0);
		nodeTags = vector<size_t>();
		map<int, MeshVertex<Dim>*> pointTagVertex;
		for (int i = 0; i < entitiesDim0Tags.size(); i++)
		{
			// Geometric point
			int pointTag = entitiesDim0Tags[i].second;
			if (_testCase)
			{
				auto it = find(_testCase->GeometricPointExclusionList.begin(), _testCase->GeometricPointExclusionList.end(), pointTag);
				if (it != _testCase->GeometricPointExclusionList.end())
					continue; // point excluded
			}

			// Find associated vertex
			gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord, 0, pointTag, true, false);
			if (nodeTags.size() == 0)
			{
				// TODO - I'm commmenting this warning because it always occurs and it's annoying.
				// Find out why and uncomment.
				
				//Utils::Warning("No mesh vertex has been found for the geometric point " + to_string(pointTag) + ". If it corresponds to the center of a circle, you can ignore this warning.");
				continue;
			}
			else if (nodeTags.size() > 1)
				Utils::Warning("Multiple mesh vertices have been found for the geometric point " + to_string(pointTag) + ". Using the first one.");
			
			MeshVertex<Dim>* v = GetVertexFromGMSHTag(nodeTags[0]);
			if (!v)
				Utils::Warning("The GMSH node of tag " + to_string(nodeTags[0]) + " corresponding to the geometric point " + to_string(pointTag) + " has not been found in the node list given by GMSH.");
			this->_geometricVertices.insert(v);
			pointTagVertex.insert({ pointTag , v });
		}
		cout << this->_geometricVertices.size() << " geometric points to conserve found." << endl;

		//------------------------//
		//   Re-entrant corners   // 
		//------------------------//

		if (Utils::ProgramArgs.Solver.MG.ReEntrantCornerManagement != ReEntrantCornerMgmt::Disabled && this->_reEntrantCorners.empty() && _testCase)
		{
			for (auto it = _testCase->ReEntrantGeometricPoints.begin(); it != _testCase->ReEntrantGeometricPoints.end(); it++)
			{
				PhysicalGroup<Dim>* phyPart = this->GetPhysicalGroup(it->first);
				if (!phyPart)
					Utils::FatalError("Unknown PhysicalPart '" + it->first + "'. Check ReEntrantGeometricPoints in the test case.");

				auto itRC = this->_reEntrantCorners.find(phyPart);
				if (itRC == this->_reEntrantCorners.end())
					this->_reEntrantCorners.insert({ phyPart, {} });

				for (int pointTag : it->second)
				{
					MeshVertex<Dim>* v = pointTagVertex.at(pointTag);
					this->_reEntrantCorners.at(phyPart).push_back(v);
				}
			}

			if (_testCase->Code().compare("magnetism") == 0)
			{
				double smallCircleR = 0.4;
				double mediumCircleR = 0.5;
				double bigCircleR = 0.8;
				for (Vertex* v : this->Vertices)
				{
					if (abs(v->X*v->X + v->Y*v->Y - smallCircleR * smallCircleR) < 1e-8)
						this->_reEntrantCorners.at(this->PhysicalParts[1]).push_back(v); // Middle
					else if (abs(v->X*v->X + v->Y*v->Y - bigCircleR * bigCircleR) < 1e-8)
						this->_reEntrantCorners.at(this->PhysicalParts[0]).push_back(v); // Exterior
					else if (abs(v->X*v->X + v->Y*v->Y - mediumCircleR * mediumCircleR) < 1e-8 && abs(v->X) < 0.2 && v->Y > 0)
						this->_reEntrantCorners.at(this->PhysicalParts[0]).push_back(v); // Exterior
				}
			}
		}

		//----------//
		// Elements //
		//----------//

		cout << "Building elements..." << endl;

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


		set<int> meshElementTypes;

		map<int, size_t> elementsAlreadyInMemory; // the key (int) is the GMSH elemType
		_average_h = 0;

		for (int k = 0; k < entitiesDimTags.size(); k++)
		{
			int entityTag = entitiesDimTags[k].second;

			// Get physical group of this entity
			vector<int> physicalTags;
			gmsh::model::getPhysicalGroupsForEntity(Dim, entityTag, physicalTags);

			PhysicalGroup<Dim>* physicalPart = nullptr;
			// If no physical group, affect the default one, created above
			if (physicalTags.empty())
			{
				physicalPart = this->GetPhysicalGroup(DefaultPhyGroup);
				if (!physicalPart)
					Utils::FatalError("Entity " + to_string(entityTag) + " has no physical part. Check GMSH file.");
			}
			else
			{
				if (physicalTags.size() > 1)
					Utils::FatalError("Entity " + to_string(entityTag) + " must have only one physical part (" + to_string(physicalTags.size()) + " found). Check GMSH file.");

				physicalPart = this->GetPhysicalGroup(physicalTags[0]);
				if (!physicalPart)
				{
					// Id not found, search by name
					string phyName;
					gmsh::model::getPhysicalName(Dim, physicalTags[0], phyName);
					physicalPart = this->GetPhysicalGroup(phyName);
					if (physicalPart)
						Utils::FatalError("Entity " + to_string(entityTag) + " has an unknown physical part (id=" + to_string(physicalTags[0]) + ", name=\"" + phyName + "\"). However it has the same name of a known one. If the GMSH file has changed, there may be discrepencies with the meshes stored in cache. Try emptying the cache for this geometry or use the option -no-cache.");
					else
						Utils::FatalError("Entity " + to_string(entityTag) + " has an unknown physical part (id=" + to_string(physicalTags[0]) + ", name=\"" + phyName + "\"). Check GMSH file.");
				}
			}

			// Get elements in this entity
			vector<int> entityElementTypes;
			vector<vector<size_t>> elementTags;
			vector<vector<size_t>> elemNodeTags;
			gmsh::model::mesh::getElements(entityElementTypes, elementTags, elemNodeTags, Dim, entityTag);

			if (entityElementTypes.empty())
				Utils::FatalError("No element type returned by GMSH for entity " + to_string(entityTag) + "!");

			for (int elemType : entityElementTypes)
				meshElementTypes.insert(elemType);

			struct ChunkResult
			{
				vector<Element<Dim>*> Elements;
				double h = -1;
				double regularity = 2;
				double sum_h = 0;
			};

			for (size_t i = 0; i < entityElementTypes.size(); i++)
			{
				int elemType = entityElementTypes[i];
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
						chunk->Results.sum_h     += e->Diameter();
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
						_average_h += chunk.sum_h;
					});

				elementsAlreadyInMemory[elemType] += elements.size();
			}
		}

		_average_h /= this->Elements.size();

		//-----------//
		//   Faces   //
		//-----------//

		cout << "Building faces..." << endl;

		// Create all the faces
		BigNumber faceNumber = 0;
		for (int elemType : meshElementTypes)
			CreateFaces(elemType, faceNumber);

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

				BoundaryGroup* boundaryPart = this->GetBoundaryGroup(physicalTags[0]);
				if (!boundaryPart)
				{
					// Id not found, search by name
					string phyName;
					gmsh::model::getPhysicalName(Dim-1, physicalTags[0], phyName);
					boundaryPart = this->GetBoundaryGroup(phyName);
					if (boundaryPart)
						Utils::FatalError("Entity " + to_string(entityTag) + " has an unknown boundary part (id=" + to_string(physicalTags[0]) + ", name=\"" + phyName + "\"). However it has the same name of a known one. If the GMSH file has changed, there may be discrepencies with the meshes stored in cache. Try emptying the cache for this geometry or use the option -no-cache.");
					else
						Utils::FatalError("Entity " + to_string(entityTag) + " has an unknown boundary part (id=" + to_string(physicalTags[0]) + ", name=\"" + phyName + "\"). Check GMSH file.");
				}

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


		if (this->ComesFrom.CS == H_CoarsStgy::None)
		{
			for (auto it = this->_reEntrantCorners.begin(); it != this->_reEntrantCorners.end(); it++)
			{
				PhysicalGroup<Dim>* phyPart = it->first;
				for (Vertex* v : it->second)
				{
					MeshVertex<Dim>* rc = dynamic_cast<MeshVertex<Dim>*>(v);
					if (!rc)
						break;
					for (Element<Dim>* e : rc->Elements)
					{
						if (e->PhysicalPart == phyPart)
							e->HasReEntrantCorner = true;
					}
				}
			}
		}

	}

protected:
	inline MeshVertex<Dim>* GetVertexFromGMSHTag(BigNumber nodeTag)
	{
		auto it = _vertexExternalNumbers.find(nodeTag);
		if (it != _vertexExternalNumbers.end())
			return it->second;
		else
			return nullptr;
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

		if (physicalGroups.empty())
			cout << "No physical parts found." << endl;
		else
		{
			cout << physicalGroups.size() << " physical parts found: ";
			for (PhysicalGroup<Dim>* pp : physicalGroups)
			{
				cout << "(" << pp->Id << ") " << pp->Name;
				if (pp != physicalGroups.back())
					cout << ", ";
			}
			cout << "." << endl;
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

		if (physicalGroups.empty())
			cout << "No boundary parts found." << endl;
		else
		{
			cout << physicalGroups.size() << " boundary parts found: ";
			for (BoundaryGroup* bp : physicalGroups)
			{
				cout << "(" << bp->Id << ") " << bp->Name;
				if (bp != physicalGroups.back())
					cout << ", ";
			}
			cout << "." << endl;
		}

		return physicalGroups;
	}

	// Defined at the end of the file for each Dim
	void AllocateContiguousMemoryForElements(int elemType, BigNumber numberOfElements) { assert(false); }

	// Defined at the end of the file for each Dim
	void AllocateContiguousMemoryForFaces(GMSHFaceTypes faceType, BigNumber numberOfFaces) { assert(false); }

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

	double AverageH() override
	{
		return _average_h;
	}

	double Regularity() override
	{
		return _regularity;
	}

	void CoarsenMesh(H_CoarsStgy elemCoarseningStgy, FaceCoarseningStrategy faceCoarseningStgy, FaceCollapsing boundaryFaceCollapsing, double coarseningFactor) override
	{
		if (Utils::IsRefinementStrategy(elemCoarseningStgy))
			return;
		else if (elemCoarseningStgy == H_CoarsStgy::IndependentRemeshing)
			IndependentRemesh(round(coarseningFactor));
		else
			PolyhedralMesh<Dim>::CoarsenMesh(elemCoarseningStgy, faceCoarseningStgy, boundaryFaceCollapsing, coarseningFactor);
	}

	void RefineMesh(H_CoarsStgy strategy) override
	{
		if (this->FineMesh)
			assert(false && "Mesh already refined!");

		if (strategy == H_CoarsStgy::GMSHSplittingRefinement)
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
		fineMesh->ComesFrom.CS = H_CoarsStgy::GMSHSplittingRefinement;
#ifdef ENABLE_3D
		if (dynamic_cast<TetrahedralElement*>(this->Elements[0]))
		{
			fineMesh->ComesFrom.nFineElementsByCoarseElement = 8;
			fineMesh->ComesFrom.nFineFacesAddedByCoarseElement = 8;
			fineMesh->ComesFrom.nFineFacesByKeptCoarseFace = 4;
		}
#endif // ENABLE_3D
		fineMesh->Build();

		// Save the fine mesh in a temporary file because we're going to reload the coarse one
		gmsh::write(fine_mesh_tmp_file);

		// Reloading the now coarse mesh to be able to link its objects with the finer ones
		gmsh::open(coarse_mesh_tmp_file);

		// Linking fine/coarse elements
		for (Element<Dim>* fine : fineMesh->Elements)
		{
			Element<Dim>* coarse = LocateElementThatEmbedsMostOf(fine);
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
	void IndependentRemesh(int coarseningFactor)
	{
		if (this->N() == 1)
			return;

		cout << "Coarsening by independent remeshing" << endl;

		GMSHMesh<Dim>* fineMesh = this;

		GMSHMesh<Dim>* coarseMesh = new GMSHMesh<Dim>(_testCase, _gmshFilePath, _description, _fileNamePart, fineMesh->N() / coarseningFactor, false);
		this->InitializeCoarsening(coarseMesh);
		coarseMesh->ComesFrom.CS = H_CoarsStgy::IndependentRemeshing;
		coarseMesh->Build();

		for (Element<Dim>* fine : fineMesh->Elements)
		{
			Element<Dim>* coarse = coarseMesh->LocateElementThatEmbedsMostOf(fine);
			if (!coarse->IsInSamePhysicalPartAs(fine))
				Utils::Warning("Coarse and fine elements aren't located in the same physical part. Something is wrong...");
			coarse->FinerElements.push_back(fine);
			fine->CoarserElement = coarse;
			fine->IsFullyEmbeddedInCoarseElement = false;
		}

		CloseGMSH();

		int nCoarseWithoutFine = 0;
		for (Element<Dim>* coarse : coarseMesh->Elements)
		{
			//if (coarse->FinerElements.size() > 5 || coarse->FinerElements.size() < 2)
			/*if (coarse->Id % (coarseMesh->Elements.size()/13) == 0)
			{
				cout << "%---------- coarse " << endl;
				coarse->ExportToMatlab("b");
				for (auto fe : coarse->FinerElements)
				{
					cout << "% fine " << endl;
					fe->ExportToMatlab("m");
				}
			}*/

			if (coarse->FinerElements.empty())
			{
				nCoarseWithoutFine++;

				// Associate the closest fine element
				/*Element<Dim>* closestFine = nullptr;
				Element<Dim>* coarseToClosestFine = nullptr;
				double closestDistance = -1;
				for (Element<Dim>* n : coarse->VertexNeighbours())
				{
					for (Element<Dim>* fine : n->FinerElements)
					{
						double distance = Vect<Dim>(coarse->Center(), fine->Center()).norm();
						if (!closestFine || distance < closestDistance)
						{
							closestFine = fine;
							coarseToClosestFine = n;
							closestDistance = distance;
						}
					}
				}
				auto it = find(coarseToClosestFine->FinerElements.begin(), coarseToClosestFine->FinerElements.end(), closestFine);
				assert(it != coarseToClosestFine->FinerElements.end());
				coarseToClosestFine->FinerElements.erase(it);
				//assert(!coarseToClosestFine->FinerElements.empty());
				coarse->FinerElements.push_back(closestFine);
				closestFine->CoarserElement = coarse;*/
			}
		}
		if (nCoarseWithoutFine > 0)
			Utils::Warning(to_string(nCoarseWithoutFine) + " coarse elements do not have any fine element.");

		this->FinalizeCoarsening();
	}


	Element<Dim>* LocateElementThatEmbedsMostOf(Element<Dim>* fineElement)
	{
		DomPoint center = fineElement->Center();
		size_t coarseElementTag = LocateGMSHElementContaining(center, true);
		if (coarseElementTag > 0)
		{
			Element<Dim>* ce = GetElementFromGMSHTag(coarseElementTag);
			assert(ce->Contains(center));
			return ce;
		}
		else
		{
			Utils::Warning("No coarse element contains the barycenter of this fine element (Id=" + to_string(fineElement->Id) + ", Number=" + to_string(fineElement->Number) + ", IsOnBoundary=" + to_string(fineElement->IsOnBoundary()) + "). Getting the closest.");
			coarseElementTag = LocateGMSHElementContaining(center, false);
			Element<Dim>* ce = GetElementFromGMSHTag(coarseElementTag);

			if (!ce->IsInSamePhysicalPartAs(fineElement))
				Utils::FatalError("TODO: find closest element in the same PhysicalPart");

			Element<Dim>* closest = ce;
			double closestDistance = Vect<Dim>(center, ce->Center()).norm();
			for (Element<Dim>* e : ce->VertexNeighbours())
			{
				double distance = Vect<Dim>(center, e->Center()).norm();
				if (distance < closestDistance)
				{
					closest = e;
					closestDistance = distance;
				}
			}
			return closest;
		}
	}

	// if !strictSearch, then finds the closest element
	size_t LocateGMSHElementContaining(const DomPoint& p, bool strictSearch)
	{
		// outputs of the gmsh function
		size_t elementTag;
		int elementType;
		vector<size_t> elementNodes;
		double u, v, w;

		// If this function fails, check that a Physical Surface (if 2D) is defined in the file (if Physical Lines exist, then a Physical Surface must also exist!).
		try
		{
			gmsh::model::mesh::getElementByCoordinates(p.X, p.Y, p.Z, elementTag,
				elementType, elementNodes, u, v, w, // useless parameters for us
				Dim, strictSearch);
		}
		catch (...)
		{
			return 0;
		}

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
			myMesh->CoarsenMesh(H_CoarsStgy::StandardCoarsening, FaceCoarseningStrategy::InterfaceCollapsing, FaceCollapsing::OnlyCollinear, 2);
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
#ifdef ENABLE_2D
		if (this->_fileNamePart.compare("gmsh_tri") == 0)
		{
			BigNumber n = sqrt(this->Elements.size() / 2);
			Square_TriangularMesh* myMesh = new Square_TriangularMesh(n, n);

			RenumberLike(dynamic_cast<Mesh<Dim>*>(myMesh));
			delete myMesh;
		}
		else
			assert(false);
#endif // ENABLE_2D
	}

	void ExportToGMSH_Elements(FunctionalBasis<Dim>* basis, const Vector &coeffs, const string& outputFilePathPrefix, const string& viewName, double tolerance=1e-3, int maxRefinements=6, bool takeAbsoluteValue=false) override
	{
		assert(!_mshFilePath.empty());
		gmsh::initialize();
		ManageGMSHLog();

		gmsh::open(_mshFilePath);
		if (_mshFileIsTmp)
			remove(_mshFilePath.c_str());

		int viewId = gmsh::view::add(viewName);

		if (!takeAbsoluteValue && basis->BasisCode().compare(MonomialBasis<2>::Code()) == 0 && Utils::ProgramArgs.Discretization.OrthogonalizeElemBasesCode == 0)
		{
			// Refer to:
			// http://www.manpagez.com/info/gmsh/gmsh-2.4.0/gmsh_52.php
			// https://gitlab.onelab.info/gmsh/gmsh/-/blob/7cc02aa46e8fec38eac639f11952b7f68dea2c7c/tutorial/c++/x3.cpp
			// http://ww1.metal.agh.edu.pl/~banas/MMNT/gmsh.pdf (p. 76)

			// Example of quadrilateral element
			// Geometric mapping: (x, y, z) = sum_i a_i * phi_i(u, v, w)
			//     4 basis functions (phi_i)_i: phi_1 = u^0*v^0 (and w^0)
			//                                  phi_2 = u^1*v^0 (and w^0)
			//                                  phi_3 = u^0*v^1 (and w^0)
			//                                  phi_4 = u^1*v^1 (and w^0)
			//     So coefMapping = Id(4,4)
			//     and expMapping = [0 0 0]
			//						[1 0 0]
			//						[0 1 0]
			//						[1 1 0].
			//     The coefficients (a_i) are given separately (it is element-dependent).
			
			// Definition of the basis functions
			int basisSize = basis->Size();
			vector<double> coefMonomials(basisSize * basisSize); // Identity
			for (int i = 0; i < basisSize; i++)
				coefMonomials[i * basisSize + i] = 1; // because this is the monomial basis! If Legendre function, it would be the coefficients in front of its monomial decomposition.

			vector<double> expMonomials(basisSize * 3); // 3 = 3D
			auto localFunctions = basis->LocalFunctions();
			for (int i = 0; i < basisSize; i++)
			{
#ifdef ENABLE_2D
				if (Dim == 2)
				{
					TensorPolynomial2D* phi = dynamic_cast<TensorPolynomial2D*>(localFunctions[i]);
					expMonomials[i * 3 + 0] = phi->FuncX->GetDegree(); // exponent of u
					expMonomials[i * 3 + 1] = phi->FuncY->GetDegree(); // exponent of v
					expMonomials[i * 3 + 2] = 0;                       // exponent of w
				}
#endif // ENABLE_2D
#ifdef ENABLE_3D
				if (Dim == 3)
				{
					TensorPolynomial3D* phi = dynamic_cast<TensorPolynomial3D*>(localFunctions[i]);
					expMonomials[i * 3 + 0] = phi->FuncX->GetDegree(); // exponent of u
					expMonomials[i * 3 + 1] = phi->FuncY->GetDegree(); // exponent of v
					expMonomials[i * 3 + 2] = phi->FuncZ->GetDegree(); // exponent of w
				}
#endif // ENABLE_3D
			}

			// Geometric mapping and high-order data
#ifdef ENABLE_2D
			if (Dim == 2)
			{
				if (!this->_quadrilateralElements.empty())
				{
					GeometricMapping quadMapping = Quadrilateral::MappingInfo();
					gmsh::view::setInterpolationMatrices(viewId, "Quadrangle", basisSize, coefMonomials, expMonomials,
						quadMapping.NFunctions, quadMapping.Coeffs, quadMapping.Exponents);

					vector<double> data = GeoMappingAndHighOrderData(this->_quadrilateralElements, coeffs, quadMapping.NFunctions, basisSize);
					gmsh::view::addListData(viewId, "SQ", this->_quadrilateralElements.size(), data);
				}
				if (!this->_triangularElements.empty())
				{
					GeometricMapping triMapping = Triangle::MappingInfo();
					gmsh::view::setInterpolationMatrices(viewId, "Triangle", basisSize, coefMonomials, expMonomials,
						triMapping.NFunctions, triMapping.Coeffs, triMapping.Exponents);

					vector<double> data = GeoMappingAndHighOrderData(this->_triangularElements, coeffs, triMapping.NFunctions, basisSize);
					gmsh::view::addListData(viewId, "ST", this->_triangularElements.size(), data);
				}
			}
#endif // ENABLE_2D
#ifdef ENABLE_3D
			if (Dim == 3)
			{
				if (!this->_tetrahedralElements.empty())
				{
					GeometricMapping tetMapping = Tetrahedron::MappingInfo();
					gmsh::view::setInterpolationMatrices(viewId, "Tetrahedron", basisSize, coefMonomials, expMonomials,
						tetMapping.NFunctions, tetMapping.Coeffs, tetMapping.Exponents);

					vector<double> data = GeoMappingAndHighOrderData(this->_tetrahedralElements, coeffs, tetMapping.NFunctions, basisSize);
					gmsh::view::addListData(viewId, "SS", this->_tetrahedralElements.size(), data);
				}
			}
#endif // ENABLE_3D

			// Options for the high-order:
			// 
			// In order to visualize the high-order field, one must activate adaptive visualization
			gmsh::view::option::setNumber(viewId, "AdaptVisualizationGrid", 1);
			// set a visualization error threshold
			gmsh::view::option::setNumber(viewId, "TargetError", tolerance);
			// set a maximum subdivision level (GMSH does automatic mesh refinement to visualize the high-order field with the requested accuracy)
			gmsh::view::option::setNumber(viewId, "MaxRecursionLevel", maxRefinements);
		}
		else
		{
			if (basis->BasisCode().compare(MonomialBasis<Dim>::Code()) != 0 || Utils::ProgramArgs.Discretization.OrthogonalizeElemBasesCode > 0)
				Utils::Warning("A piecewise-constant function will be exported. For a high-order visualization, use non-orthogonalized monomial bases for the elements (-e-basis monomials -e-ogb 0).");

			// Set a single value in the element
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
				double value = e->EvalApproximateSolution(basis, coeffs, e->Center());
				if (takeAbsoluteValue)
					values.push_back({ abs(value) });
				else
					values.push_back({ value });
			}
			gmsh::view::addModelData(viewId, 0, modelName, "ElementData", elementTags, values);
		}


		string meshFilePath = outputFilePathPrefix + (outputFilePathPrefix.back() == '/' ? "" : "_") + viewName + ".msh";
		string dataFilePath = outputFilePathPrefix + (outputFilePathPrefix.back() == '/' ? "" : "_") + viewName + ".pos";

		gmsh::write(meshFilePath);
		gmsh::view::write(viewId, dataFilePath);

		cout << viewName << " exported for GMSH to " << dataFilePath << endl;

		// To open the GMSH window
		//gmsh::fltk::run();

		gmsh::finalize();
	}

	static void ExportToGMSH_Faces(PolyhedralMesh<Dim>* mesh, FunctionalBasis<Dim-1>* basis, const Vector& coeffs, const string& outputFilePathPrefix, const string& viewName, double tolerance = 1e-3, int maxRefinements = 6)
	{
		if (basis->BasisCode().compare(MonomialBasis<2>::Code()) != 0 || Utils::ProgramArgs.Discretization.OrthogonalizeFaceBasesCode > 0)
		{
			Utils::Error("For face visualization, use non-orthogonalized monomial bases (-f-basis monomials -f-ogb 0).");
			return;
		}
/*#ifdef ENABLE_2D
		if (Dim == 2 && mesh->_edgeFaces.empty())
		{
			Utils::Error("Face export to GMSH aborted, because mesh->_edgeFaces is empty. If it is a polygonal mesh, mesh->_edgeFaces is not filled (TODO).");
			return;
		}
#endif // ENABLE_2D*/

		gmsh::initialize();
		ManageGMSHLog();

		int viewId = gmsh::view::add(viewName);

		// Definition of the basis functions
		int basisSize = basis->Size();
		vector<double> coefMonomials(basisSize * basisSize); // Identity
		for (int i = 0; i < basisSize; i++)
			coefMonomials[i * basisSize + i] = 1; // because this is the monomial basis! If Legendre function, it would be the coefficients in front of its monomial decomposition.

		vector<double> expMonomials(basisSize * 3); // 3 = 3D
		auto localFunctions = basis->LocalFunctions();
		for (int i = 0; i < basisSize; i++)
		{
#ifdef ENABLE_2D
			if (Dim == 2)
			{
				BasisFunction<1>* phi = dynamic_cast<BasisFunction<1>*>(localFunctions[i]);
				expMonomials[i * 3 + 0] = phi->GetDegree(); // exponent of u
				expMonomials[i * 3 + 1] = 0;                // exponent of v
				expMonomials[i * 3 + 2] = 0;                // exponent of w
			}
#endif // ENABLE_2D
#ifdef ENABLE_3D
			if (Dim == 3)
			{
				TensorPolynomial2D* phi = dynamic_cast<TensorPolynomial2D*>(localFunctions[i]);
				expMonomials[i * 3 + 0] = phi->FuncX->GetDegree(); // exponent of u
				expMonomials[i * 3 + 1] = phi->FuncY->GetDegree(); // exponent of v
				expMonomials[i * 3 + 2] = 0;                       // exponent of w
			}
#endif // ENABLE_3D
		}

		// Geometric mapping and high-order data
#ifdef ENABLE_2D
		if (Dim == 2)
		{
			GeometricMapping segmentMapping = Segment::MappingInfo();
			gmsh::view::setInterpolationMatrices(viewId, "Line", basisSize, coefMonomials, expMonomials,
				segmentMapping.NFunctions, segmentMapping.Coeffs, segmentMapping.Exponents);

			if (!mesh->_edgeFaces.empty())
			{
				vector<double> data = GeoMappingAndHighOrderData(mesh->_edgeFaces, coeffs, segmentMapping.NFunctions, basisSize);
				gmsh::view::addListData(viewId, "SL", mesh->_edgeFaces.size(), data); // SL = Scalar Line
			}
			else
			{
				vector<Edge> edges;
				edges.reserve(mesh->Faces.size());
				for (Face<Dim>* f : mesh->Faces)
				{
					Edge* e = dynamic_cast<Edge*>(f);
					edges.push_back(*e);
				}
				vector<double> data = GeoMappingAndHighOrderData(edges, coeffs, segmentMapping.NFunctions, basisSize);
				gmsh::view::addListData(viewId, "SL", edges.size(), data); // SL = Scalar Line
			}
		}
#endif // ENABLE_2D
#ifdef ENABLE_3D
		if (Dim == 3)
		{
			if (!mesh->_triangularFaces.empty())
			{
				GeometricMapping triMapping = TriangleIn3D::MappingInfo();
				gmsh::view::setInterpolationMatrices(viewId, "Triangle", basisSize, coefMonomials, expMonomials,
					triMapping.NFunctions, triMapping.Coeffs, triMapping.Exponents);

				vector<double> data = GeoMappingAndHighOrderData(mesh->_triangularFaces, coeffs, triMapping.NFunctions, basisSize);
				gmsh::view::addListData(viewId, "ST", mesh->_triangularFaces.size(), data); // SL = Scalar Triangle
			}
			if (!mesh->_rectangularFaces.empty())
			{
				GeometricMapping cartMapping = CartesianShape<3, 2>::MappingInfo();
				gmsh::view::setInterpolationMatrices(viewId, "Quadrangle", basisSize, coefMonomials, expMonomials,
					cartMapping.NFunctions, cartMapping.Coeffs, cartMapping.Exponents);

				vector<double> data = GeoMappingAndHighOrderData(mesh->_rectangularFaces, coeffs, cartMapping.NFunctions, basisSize);
				gmsh::view::addListData(viewId, "SQ", mesh->_rectangularFaces.size(), data); // SQ = Scalar Quadrangle
			}
		}
#endif // ENABLE_3D

		// Options for the high-order:
		// 
		// In order to visualize the high-order field, one must activate adaptive visualization
		gmsh::view::option::setNumber(viewId, "AdaptVisualizationGrid", 1);
		// set a visualization error threshold
		gmsh::view::option::setNumber(viewId, "TargetError", tolerance);
		// set a maximum subdivision level (GMSH does automatic mesh refinement to visualize the high-order field with the requested accuracy)
		gmsh::view::option::setNumber(viewId, "MaxRecursionLevel", maxRefinements);


		string dataFilePath = outputFilePathPrefix + (outputFilePathPrefix.back() == '/' ? "" : "_") + viewName + ".pos";
		gmsh::view::write(viewId, dataFilePath);
		cout << viewName << " exported for GMSH to " << dataFilePath << endl;

		// To open the GMSH window
		//gmsh::fltk::run();

		gmsh::finalize();
	}

private:
	template<class ElementType>
	static vector<double> GeoMappingAndHighOrderData(const vector<ElementType>& elements, const Vector& coeffs, int nbMappingFunctions, int basisSize)
	{
		assert(elements.size() * basisSize == coeffs.rows());

		vector<double> data;
		data.reserve(elements.size() * (nbMappingFunctions * 3 + basisSize));
		for (const ElementType& e : elements)
		{
			// 1. Geometric mapping of the element w.r.t. the reference element
			// 
			// Instead of adding the coordinates of the vertices, we add the mapping coefficients.
			// The actual coordinates are retrieved using the mapping interpolation matrices given in gmsh::view::setInterpolationMatrices().
			vector<double> mappingCoeffs = e.Shape()->MappingCoefficients();
			std::copy(mappingCoeffs.begin(), mappingCoeffs.end(), std::back_inserter(data));

			// 2. High-order data
			// 
			// Add the coefficient of the linear combination in the basis functions.
			Vector elemCoeffs = coeffs.segment(e.Number * basisSize, basisSize);
			for (int i = 0; i < basisSize; i++)
				data.push_back(elemCoeffs[i]);
		}
		return data;
	}

public:
	/*void ExportAbsoluteValueToGMSH_Elements(FunctionalBasis<Dim>* basis, const Vector& coeffs, const string& outputFilePathPrefix, const string& viewName) override
	{
		assert(!_mshFilePath.empty());
		gmsh::initialize();
		ManageGMSHLog();

		gmsh::open(_mshFilePath);
		if (_mshFileIsTmp)
			remove(_mshFilePath.c_str());

		int viewId = gmsh::view::add(viewName);

		// Set a single value in the element
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
			double value = e->EvalApproximateSolution(basis, coeffs, e->Center());
			values.push_back({ abs(value) });
		}
		gmsh::view::addModelData(viewId, 0, modelName, "ElementData", elementTags, values);


		string meshFilePath = outputFilePathPrefix + (outputFilePathPrefix.back() == '/' ? "" : "_") + viewName + ".msh";
		string dataFilePath = outputFilePathPrefix + (outputFilePathPrefix.back() == '/' ? "" : "_") + viewName + ".pos";

		gmsh::write(meshFilePath);
		gmsh::view::write(viewId, dataFilePath);

		cout << viewName << " exported for GMSH to " << dataFilePath << endl;

		// To open the GMSH window
		//gmsh::fltk::run();

		gmsh::finalize();
	}*/

	void ExportToGMSH_Nodes(const Vector& nodeValues, const string& outputFilePathPrefix, const string& suffix) override
	{
		assert(!_mshFilePath.empty());
		gmsh::initialize();
		ManageGMSHLog();

		gmsh::open(_mshFilePath);
		if (_mshFileIsTmp)
			remove(_mshFilePath.c_str());

		int viewId = gmsh::view::add(suffix);

		vector<std::string> modelNames;
		gmsh::model::list(modelNames);
		string modelName = modelNames[modelNames.size() - 1];

		vector<size_t> nodeTags;
		vector<vector<double>> values;
		nodeTags.reserve(this->Vertices.size());
		values.reserve(this->Vertices.size());

		for (auto const& [nodeTag, v] : _vertexExternalNumbers)
		{
			nodeTags.push_back(nodeTag);
			values.push_back({ nodeValues[v->Number] });
		}
		gmsh::view::addModelData(viewId, 0, modelName, "NodeData", nodeTags, values);

		string meshFilePath = outputFilePathPrefix + (outputFilePathPrefix.back() == '/' ? "" : "_") + suffix + ".msh";
		string dataFilePath = outputFilePathPrefix + (outputFilePathPrefix.back() == '/' ? "" : "_") + suffix + ".pos";

		gmsh::write(meshFilePath);
		gmsh::view::write(viewId, dataFilePath);

		cout << suffix << " exported for GMSH to " << dataFilePath << endl;
		gmsh::finalize();
	}

	void ExportToGMSH_Nodes(DomFunction f, const string& outputFilePathPrefix, const string& suffix)
	{
		assert(!_mshFilePath.empty());
		gmsh::initialize();
		ManageGMSHLog();

		gmsh::open(_mshFilePath);
		if (_mshFileIsTmp)
			remove(_mshFilePath.c_str());

		int viewId = gmsh::view::add(suffix);

		vector<std::string> modelNames;
		gmsh::model::list(modelNames);
		string modelName = modelNames[modelNames.size() - 1];

		vector<size_t> nodeTags;
		vector<vector<double>> values;
		nodeTags.reserve(this->Vertices.size());
		values.reserve(this->Vertices.size());

		for (auto const& [nodeTag, v] : _vertexExternalNumbers)
		{
			nodeTags.push_back(nodeTag);
			values.push_back({ f(*v) });
		}
		gmsh::view::addModelData(viewId, 0, modelName, "NodeData", nodeTags, values);

		string meshFilePath = outputFilePathPrefix + (outputFilePathPrefix.back() == '/' ? "" : "_") + suffix + ".msh";
		string dataFilePath = outputFilePathPrefix + (outputFilePathPrefix.back() == '/' ? "" : "_") + suffix + ".pos";

		gmsh::write(meshFilePath);
		gmsh::view::write(viewId, dataFilePath);

		cout << suffix << " exported for GMSH to " << dataFilePath << endl;
		gmsh::finalize();
	}

	void ExportToGMSH_Elements(DomFunction f, const string& outputFilePath, const string& dataName)
	{
		assert(!_mshFilePath.empty());
		gmsh::initialize();
		ManageGMSHLog();

		gmsh::open(_mshFilePath);
		if (_mshFileIsTmp)
			remove(_mshFilePath.c_str());

		int viewId = gmsh::view::add(dataName);

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
			double value = f(e->Center());
			values.push_back({ value });
		}
		gmsh::view::addModelData(viewId, 0, modelName, "ElementData", elementTags, values);

		string dataFilePath = outputFilePath + ".pos";

		gmsh::view::write(viewId, dataFilePath);

		cout << dataName << " exported for GMSH to " << dataFilePath << endl;
		gmsh::finalize();
	}

};

//-------------//
// 2D elements //
//-------------//

#ifdef ENABLE_2D

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
		Utils::FatalError("GMSH element type not managed.");
	return e;
}
#endif // ENABLE_2D

//-------------//
// 3D elements //
//-------------//

#ifdef ENABLE_3D

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
	else if (elemType == GMSH_Hexahedron && this->FileNamePart().compare("gmsh_cart") == 0)
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
#endif // ENABLE_3D

//--------------//
//   2D faces   //
//--------------//

#ifdef ENABLE_2D

template<>
void GMSHMesh<2>::AllocateContiguousMemoryForElements(int elemType, BigNumber numberOfElements)
{
	this->Elements.reserve(numberOfElements);

	if (elemType == GMSH_Quadrilateral)
		this->_quadrilateralElements = vector<QuadrilateralElement>(numberOfElements);
	else if (elemType == GMSH_Triangle)
		this->_triangularElements = vector<TriangularElement>(numberOfElements);
	else
		Utils::FatalError("GMSH element type not managed.");
}

template<>
void GMSHMesh<2>::AllocateContiguousMemoryForFaces(GMSHFaceTypes faceType, BigNumber numberOfFaces)
{
	this->Faces.reserve(numberOfFaces);
	this->_edgeFaces.reserve(numberOfFaces);
}

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
		if (!v1)
			Utils::FatalError("No GMSH node of tag " + to_string(edgeNodes[j]));
		if (!v2)
			Utils::FatalError("No GMSH node of tag " + to_string(edgeNodes[j + 1]));

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
#endif // ENABLE_2D

//--------------//
//   3D faces   //
//--------------//

#ifdef ENABLE_3D

template<>
void GMSHMesh<3>::AllocateContiguousMemoryForElements(int elemType, BigNumber numberOfElements)
{
	this->Elements.reserve(numberOfElements);

	if (elemType == GMSH_Tetrahedron)
		this->_tetrahedralElements = vector<TetrahedralElement>(numberOfElements);
	else if (elemType == GMSH_Hexahedron)
		this->_parallelepipedElements = vector<ParallelepipedElement>(numberOfElements);
	else
		Utils::FatalError("GMSH element type not managed.");
}

template<>
void GMSHMesh<3>::AllocateContiguousMemoryForFaces(GMSHFaceTypes faceType, BigNumber numberOfFaces)
{
	this->Faces.reserve(numberOfFaces);

	if (faceType == GMSHFaceTypes::GMSH_TriangleFace)
		this->_triangularFaces.reserve(numberOfFaces);
	else
		Utils::FatalError("GMSH element type not managed.");
}

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
		else if (faceType == GMSHFaceTypes::GMSH_QuadrilateralFace && this->FileNamePart().compare("gmsh_cart") == 0)
		{
			Vertex* v1 = vertices[0]; // (0, 1, 1)
			Vertex* v2 = vertices[1]; // (1, 1, 1)
			Vertex* v3 = vertices[2]; // (1, 1, 0)
			Vertex* v4 = vertices[3]; // (0, 1, 0)

			CartesianShapeOrientation orientation = CartesianShapeOrientation::None;
			Vertex* origin;
			Vertex* vertex1;
			Vertex* vertex2;
			Vertex* vertex3;

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
				vertex3 = v2;
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
				vertex3 = v3;
			}
			else if (minNorm == OV2.norm())
			{
				origin = v2;
				vertex1 = v3;
				vertex2 = v1;
				vertex3 = v4;
			}
			else if (minNorm == OV3.norm())
			{
				origin = v3;
				vertex1 = v2;
				vertex2 = v4;
				vertex3 = v1;
			}
			else if (minNorm == OV4.norm())
			{
				origin = v4;
				vertex1 = v1;
				vertex2 = v3;
				vertex3 = v2;
			}
			else
				assert(false);

			if (neighbours.size() == 1)
				face = new RectangularFace(faceNumber++, origin, vertex1, vertex2, vertex3, neighbours[0], orientation);
			else if (neighbours.size() == 2)
				face = new RectangularFace(faceNumber++, origin, vertex1, vertex2, vertex3, neighbours[0], neighbours[1], orientation);
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
#endif // ENABLE_3D


#ifdef ENABLE_2D
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
	return nullptr;
}
#endif // ENABLE_2D

#ifdef ENABLE_3D
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
	return nullptr;
}
#endif // ENABLE_3D

#ifdef ENABLE_1D
template <>
bool GMSHMesh<1>::GMSHLogEnabled = false;
#endif
#ifdef ENABLE_2D
template <>
bool GMSHMesh<2>::GMSHLogEnabled = false;
#endif
#ifdef ENABLE_3D
template <>
bool GMSHMesh<3>::GMSHLogEnabled = false;
#endif

#ifdef ENABLE_1D
template <>
bool GMSHMesh<1>::UseCache = true;
#endif
#ifdef ENABLE_2D
template <>
bool GMSHMesh<2>::UseCache = true;
#endif
#ifdef ENABLE_3D
template <>
bool GMSHMesh<3>::UseCache = true;
#endif