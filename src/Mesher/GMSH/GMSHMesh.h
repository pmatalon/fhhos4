#pragma once
#include <array>
#include <gmsh.h>
#include "../../Mesh/PolyhedralMesh.h"
#include "../../Utils/Communicator.h"
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
	vector<size_t> _vertexTags;    // the GMSH tag of each vertex, by number
	bool _part = false;            // the part of a process (OpenPart())
	map<size_t, int> _ghostRanks;  // in a part, the process of each ghost cell kept, by GMSH tag
	double _h = -1;
	double _average_h = -1;
	double _regularity = 1;
	BigNumber _N;
private:
	string _mshFilePath;
	bool _mshFileIsTmp = false;
public:
	inline static bool GMSHLogEnabled = false;
	inline static bool UseCache = true;

	// On several processes (MPI), the part of the mesh of the process (GMSH's partition): its own elements and the other
	// elements of their faces (the ghosts). The process of each element and its global number (contiguous by process, in
	// rank order), by number. Empty on one process.
	vector<int> ElementRanks;
	vector<BigNumber> ElementGlobalNumbers;

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

		// On several processes: the partition of the mesh (OpenPart()), in the cache if a previous run made it for as many
		// processes (the first process then doesn't load the mesh)
		const Communicator& comm = Communicator::Current();
		string inCache = Mesh<Dim>::CacheDirectory + _geometryDescription + (n == 0 ? "" : "_N" + to_string(n));
		string partition = inCache + "_" + to_string(comm.Size()) + "parts";
		bool partitionInCache = comm.IsDistributed() && UseCache && comm.Any(comm.IsRoot() && filesystem::exists(partition));

		bool tmpGeoFileCreated = false;
		string mshFile;
		string extension = FileSystem::Extension(gmshFile);
		if (comm.IsDistributed() && (!comm.IsRoot() || partitionInCache))
		{
			// only its part (OpenPart())
		}
		else if (extension.compare("geo") == 0)
		{
			string mshFileInCache = inCache + ".msh";

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

		if (comm.IsDistributed())
			OpenPart(comm, partition, partitionInCache);
		if (buildMesh)
			Build();
		if (comm.IsDistributed())
			ExchangeGhosts(comm);
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

	// On several processes, after the numbering of the faces by the problem: the global numbers of the faces of number
	// < nFaces (those of the system), contiguous by process in rank order, each face belonging to the smaller process of
	// its elements; nOwned: those of the process. The others are asked to their process: first the faces of its elements
	// (both elements of such a face are in the part), then those of ghosts only, to the process of the ghost (which by
	// then knows the numbers of all the faces of its elements). Collective.
	vector<BigNumber> GlobalFaceNumbers(BigNumber nFaces, BigNumber& nOwned)
	{
		const Communicator& comm = Communicator::Current();
		int rank = comm.Rank();
		auto isOwn = [&](Element<Dim>* e) { return e && ElementRanks[e->Number] == rank; };
		auto faceRank = [&](Face<Dim>* f) { return min(ElementRanks[f->Element1->Number], f->Element2 ? ElementRanks[f->Element2->Number] : rank); };
		vector<Face<Dim>*> faces(nFaces); // by number
		for (Face<Dim>* f : this->Faces)
			if (f->Number < nFaces)
				faces[f->Number] = f;
		vector<Face<Dim>*> owned, ofOwnElements, ofGhostsOnly;
		for (Face<Dim>* f : faces)
		{
			if (!isOwn(f->Element1) && !isOwn(f->Element2))
				ofGhostsOnly.push_back(f);
			else if (faceRank(f) == rank)
				owned.push_back(f);
			else
				ofOwnElements.push_back(f);
		}
		nOwned = owned.size();
		BigNumber begin = comm.ExclusiveScan(nOwned);
		vector<BigNumber> global(nFaces, -1);
		for (BigNumber i = 0; i < nOwned; i++)
			global[owned[i]->Number] = begin + i;
		AskFaceNumbers(comm, ofOwnElements, faceRank, global);
		AskFaceNumbers(comm, ofGhostsOnly, [&](Face<Dim>* f) { return ElementRanks[f->Element1->Number]; }, global);
		return global;
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
				Utils::FatalError(FileSystem::DataFileNotFound(geoFile, { geoFile, Mesh<Dim>::MeshDirectory + geoFile,
					Mesh<Dim>::MeshDirectory + "2D/" + geoFile, Mesh<Dim>::MeshDirectory + "3D/" + geoFile }));
		}
		return geoFile;
	}

	// Copy of the GMSH script with N = n, in a temporary file (the meshes directory may be read-only). The relative paths
	// of its Include statements, from the directory of the script, are made absolute.
	static string SetN(string geoFile, BigNumber n)
	{
		assert(n > 0);
		string geoFileWithN = FileSystem::TemporaryFile(FileSystem::FileNameWithoutExtension(geoFile) + "_N" + to_string(n) + ".geo");
		filesystem::path geoDirectory = filesystem::absolute(geoFile).parent_path();
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
			size_t include = line.find_first_not_of(" \t");
			if (line.find("N =") == 0 || line.find("N=") == 0)
			{
				NFound = true;
				gmshScriptWithN << "N = " << n << ";\r" << endl;
			}
			else if (include != string::npos && line.compare(include, 9, "Include \"") == 0) // Include "file";
			{
				size_t begin = include + 9;
				size_t end = line.find('"', begin);
				if (end != string::npos && filesystem::path(line.substr(begin, end - begin)).is_relative())
					line = line.substr(0, begin) + (geoDirectory / line.substr(begin, end - begin)).string() + line.substr(end);
				gmshScriptWithN << line << endl;
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

	// On several processes: the mesh partitioned by GMSH (METIS) into one part per process, with ghost cells (the cells of
	// other parts that share a node with the part), in GMSH's split files in the directory partition (in the cache
	// directory, which must be shared by the nodes: -cache-dir), by the first process unless inCache. Kept for the next
	// runs on as many processes (written in a temporary directory renamed once complete: a run never reads an incomplete
	// partition), or deleted with -no-cache. Each process opens its part, and keeps the ghost cells that share a face
	// with its cells.
	void OpenPart(const Communicator& comm, string partition, bool inCache)
	{
		vector<char> path;
		if (comm.IsRoot() && inCache)
			cout << "Partition of the mesh into " << comm.Size() << " parts read from the cache: " << partition << endl;
		else if (comm.IsRoot())
		{
			cout << "Partition of the mesh into " << comm.Size() << " parts by GMSH..." << endl;
			gmsh::option::setNumber("Mesh.PartitionCreateTopology", 0);
			gmsh::option::setNumber("Mesh.PartitionCreatePhysicals", 0);
			gmsh::option::setNumber("Mesh.PartitionCreateGhostCells", 1);
			gmsh::option::setNumber("Mesh.PartitionSplitMeshFiles", 1);
			gmsh::model::mesh::partition(comm.Size());
			string written = partition + "_tmp" + to_string(getpid());
			filesystem::create_directories(written);
			gmsh::write(written + "/part.msh"); // one file per part: part_1.msh, part_2.msh...
			error_code error;
			if (!UseCache)
				partition = written;
			else
				filesystem::rename(written, partition, error);
			if (error && filesystem::exists(partition)) // written by a run at the same time: kept
				filesystem::remove_all(written);
			else if (error)
				Utils::FatalError("Unable to rename " + written + " to " + partition + ": " + error.message());
		}
		if (comm.IsRoot())
			path.assign(partition.begin(), partition.end());
		comm.Broadcast(path);
		partition.assign(path.begin(), path.end());

		string file = partition + "/part_" + to_string(comm.Rank() + 1) + ".msh";
		bool found = FileSystem::FileExists(file);
		gmsh::clear();
		if (found)
			gmsh::open(file);
		int nParts = found ? gmsh::model::getNumberOfPartitions() : comm.Size();
		int otherNParts = comm.Max(nParts != comm.Size() ? nParts : 0);
		if (otherNParts > 0)
			Utils::FatalError("The partition of the mesh in " + partition + " (from a previous run) has " + to_string(otherNParts) + " parts, but the program runs on " + to_string(comm.Size()) + " processes. Delete this directory.");
		if (comm.Any(!found))
			Utils::FatalError((found ? "Some processes cannot read the partition of the mesh in " + partition : file + " not found") + ". On several nodes, the cache directory must be shared by them (-cache-dir, e.g. on the scratch of a cluster).");
		if (!UseCache)
		{
			comm.Barrier(); // (all have read their part)
			if (comm.IsRoot())
				filesystem::remove_all(partition);
		}
		_part = true;
		KeepGhostsOfFaces();
	}

	// In the part of a process: the ghost cells that share no face with a cell of the part are removed from GMSH's model;
	// the process of the others is kept (_ghostRanks)
	void KeepGhostsOfFaces()
	{
		gmsh::vectorpair entities;
		gmsh::model::getEntities(entities, Dim);
		vector<FaceKey> ownFaces;
		vector<int> ghostEntities;
		for (auto [dim, tag] : entities)
		{
			vector<size_t> ghosts;
			vector<int> parts;
			gmsh::model::mesh::getGhostElements(dim, tag, ghosts, parts);
			for (size_t i = 0; i < ghosts.size(); i++)
				_ghostRanks[ghosts[i]] = parts[i] - 1;
			if (!ghosts.empty())
			{
				ghostEntities.push_back(tag);
				continue;
			}
			vector<int> types;
			gmsh::model::mesh::getElementTypes(types, dim, tag);
			for (int type : types)
			{
				int n;
				vector<size_t> nodes = FaceNodes(type, tag, n);
				for (size_t i = 0; i < nodes.size(); i += n)
					ownFaces.push_back(Key(&nodes[i], n));
			}
		}
		sort(ownFaces.begin(), ownFaces.end());
		for (int tag : ghostEntities)
		{
			vector<size_t> removed;
			vector<int> types;
			gmsh::model::mesh::getElementTypes(types, Dim, tag);
			for (int type : types)
			{
				vector<size_t> elements, elementNodes;
				gmsh::model::mesh::getElementsByType(type, elements, elementNodes, tag);
				int n;
				vector<size_t> nodes = FaceNodes(type, tag, n);
				size_t nFaces = nodes.size() / n / elements.size();
				for (size_t e = 0; e < elements.size(); e++)
				{
					bool ofOwnFace = false;
					for (size_t f = 0; f < nFaces && !ofOwnFace; f++)
						ofOwnFace = binary_search(ownFaces.begin(), ownFaces.end(), Key(&nodes[(e * nFaces + f) * n], n));
					if (!ofOwnFace)
					{
						removed.push_back(elements[e]);
						_ghostRanks.erase(elements[e]);
					}
				}
			}
			gmsh::model::mesh::removeElements(Dim, tag, removed); // (still in GMSH's list of the ghosts: not to be read after)
		}
	}

	// On several processes, after Build(): the global numbers of the elements (own: contiguous by process; ghosts: given
	// by their process), the physical part of the ghosts and those of their faces with one element in the part that are
	// on the boundary of the domain (with their boundary part), asked to the process of each ghost: the others lie
	// between the ghost and an element of another part (not in this one). _h, the regularity and the average h: those of
	// the whole mesh.
	void ExchangeGhosts(const Communicator& comm)
	{
		int rank = comm.Rank();
		ElementRanks.assign(this->Elements.size(), rank);
		for (Element<Dim>* e : this->Elements)
		{
			auto it = _ghostRanks.find(e->Id);
			if (it != _ghostRanks.end())
				ElementRanks[e->Number] = it->second;
		}
		BigNumber nOwned = count(ElementRanks.begin(), ElementRanks.end(), rank);
		BigNumber next = comm.ExclusiveScan(nOwned);
		ElementGlobalNumbers.assign(this->Elements.size(), -1);
		double sum_h = 0;
		vector<vector<BigNumber>> requests(comm.Size());
		vector<vector<Element<Dim>*>> asked(comm.Size());
		for (Element<Dim>* e : this->Elements)
		{
			int r = ElementRanks[e->Number];
			if (r == rank)
			{
				if (!e->PhysicalPart)
					Utils::FatalError("Element " + to_string(e->Id) + " has no physical part. Check GMSH file.");
				ElementGlobalNumbers[e->Number] = next++;
				sum_h += e->Diameter();
			}
			else
			{
				requests[r].push_back(e->Id);
				asked[r].push_back(e);
			}
		}

		// Reply for each element asked: its number, its physical part, its number of faces on the boundary of the domain,
		// then for each of them its key and its boundary part (-1: none)
		vector<vector<BigNumber>> received = comm.AllToAll(requests);
		vector<vector<BigNumber>> replies(comm.Size());
		for (int r = 0; r < comm.Size(); r++)
		{
			vector<BigNumber>& reply = replies[r];
			for (BigNumber tag : received[r])
			{
				Element<Dim>* e = GetElementFromGMSHTag(tag);
				reply.push_back(ElementGlobalNumbers[e->Number]);
				reply.push_back(e->PhysicalPart->Id);
				size_t nBoundaryFaces = reply.size();
				reply.push_back(0);
				for (Face<Dim>* f : e->Faces)
				{
					if (!f->IsDomainBoundary)
						continue;
					FaceKey key = Key(f);
					reply.insert(reply.end(), key.begin(), key.end());
					reply.push_back(f->BoundaryPart ? f->BoundaryPart->Id : -1);
					reply[nBoundaryFaces]++;
				}
			}
		}
		vector<vector<BigNumber>> answers = comm.AllToAll(replies);
		for (int r = 0; r < comm.Size(); r++)
		{
			const BigNumber* a = answers[r].data();
			for (Element<Dim>* e : asked[r])
			{
				ElementGlobalNumbers[e->Number] = *a++;
				e->PhysicalPart = this->GetPhysicalGroup(*a++);
				const BigNumber* boundaryFaces = a + 1;
				a = boundaryFaces + *a * (KeySize + 1);
				for (Face<Dim>* f : e->Faces)
				{
					if (f->Element2) // (between two elements of the part)
						continue;
					FaceKey key = Key(f);
					f->IsDomainBoundary = false;
					f->BoundaryPart = nullptr;
					for (const BigNumber* b = boundaryFaces; b < a; b += KeySize + 1)
					{
						if (equal(key.begin(), key.end(), b))
						{
							f->IsDomainBoundary = true;
							f->BoundaryPart = b[KeySize] < 0 ? nullptr : this->GetBoundaryGroup(b[KeySize]);
						}
					}
				}
			}
		}
		this->BoundaryFaces.clear();
		this->InteriorFaces.clear();
		this->FillBoundaryAndInteriorFaceLists();

		_h = comm.Max(_h);
		_regularity = comm.Min(_regularity);
		_average_h = comm.Sum(sum_h) / comm.Sum(nOwned);
	}

	// A face known by several processes: the GMSH tags of its vertices, sorted (then 0)
	static constexpr int KeySize = 4;
	using FaceKey = array<BigNumber, KeySize>;

	static FaceKey Key(const size_t* nodeTags, int nNodes)
	{
		FaceKey key{};
		copy(nodeTags, nodeTags + nNodes, key.begin());
		sort(key.begin(), key.begin() + nNodes);
		return key;
	}

	FaceKey Key(Face<Dim>* f) const
	{
		vector<size_t> tags;
		for (Vertex* v : f->Vertices())
			tags.push_back(_vertexTags[v->Number]);
		return Key(tags.data(), tags.size());
	}

	// The face whose vertices have the GMSH tags of key (nullptr if none)
	Face<Dim>* FaceWithKey(const BigNumber* key)
	{
		vector<MeshVertex<Dim>*> vertices;
		for (int i = 0; i < KeySize && key[i] != 0; i++)
		{
			MeshVertex<Dim>* v = GetVertexFromGMSHTag(key[i]);
			if (!v)
				return nullptr;
			vertices.push_back(v);
		}
		return vertices.empty() ? nullptr : this->ExistingFaceWithVertices(vertices);
	}

	// The global numbers (global, by face number) of faces, asked to the process owner(f) of each, which knows them.
	// Collective.
	template <class Owner>
	void AskFaceNumbers(const Communicator& comm, const vector<Face<Dim>*>& faces, Owner owner, vector<BigNumber>& global)
	{
		vector<vector<BigNumber>> requests(comm.Size());
		vector<vector<Face<Dim>*>> asked(comm.Size());
		for (Face<Dim>* f : faces)
		{
			int r = owner(f);
			FaceKey key = Key(f);
			requests[r].insert(requests[r].end(), key.begin(), key.end());
			asked[r].push_back(f);
		}
		vector<vector<BigNumber>> received = comm.AllToAll(requests);
		vector<vector<BigNumber>> replies(comm.Size());
		for (int r = 0; r < comm.Size(); r++)
		{
			for (size_t q = 0; q < received[r].size(); q += KeySize)
			{
				Face<Dim>* f = FaceWithKey(&received[r][q]);
				replies[r].push_back(f && f->Number < (BigNumber)global.size() ? global[f->Number] : -1);
			}
		}
		vector<vector<BigNumber>> answers = comm.AllToAll(replies);
		bool unknown = false;
		for (int r = 0; r < comm.Size(); r++)
		{
			for (size_t q = 0; q < asked[r].size(); q++)
			{
				global[asked[r][q]->Number] = answers[r][q];
				unknown = unknown || answers[r][q] < 0;
			}
		}
		if (comm.Any(unknown))
			Utils::FatalError("GMSHMesh: a face of the system is unknown to the process asked for its number.");
	}

	// Defined at the end of the file for each Dim: the GMSH tags of the vertices of the faces of the elements of type
	// elemType of the entity tag (-1: all), by element (in GMSH's order), nFaceNodes per face
	static vector<size_t> FaceNodes(int elemType, int tag, int& nFaceNodes) { assert(false); return {}; }

	// In the part of a process: the vertices of a face (in the order of a cycle) in an order that depends only on their
	// GMSH tags, the same on all the processes, which then orient the face (and its basis) alike: from the smallest tag,
	// towards the smaller of its two neighbours in the cycle. Elsewhere, the order of the first element of the face.
	void OrientAlikeOnAllProcesses(vector<MeshVertex<Dim>*>& vertices) const
	{
		if (!_part)
			return;
		auto smallerTag = [this](MeshVertex<Dim>* a, MeshVertex<Dim>* b) { return _vertexTags[a->Number] < _vertexTags[b->Number]; };
		rotate(vertices.begin(), min_element(vertices.begin(), vertices.end(), smallerTag), vertices.end());
		if (vertices.size() > 2 && smallerTag(vertices.back(), vertices[1]))
			reverse(vertices.begin() + 1, vertices.end());
	}

	void Build()
	{
		//-----------------//
		// Physical groups //
		//-----------------//

		const int DefaultPhyGroup = -1;

		// The part of a process (OpenPart()): its nodes are on entities of all dimensions, the entities of the other parts
		// are empty, the ghost cells have no physical group (given by their process, ExchangeGhosts()), and the nodes of
		// some geometric points are in other parts
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
		gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord, _part ? -1 : Dim, -1, includeBoundary, returnParametricCoord);

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
				_vertexTags.push_back(nodeTags[i]);
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
			if (!v && _part)
				continue;
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
				if (!physicalPart && !_part)
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

			if (entityElementTypes.empty() && _part)
				continue;
			if (entityElementTypes.empty())
				Utils::FatalError("No element type returned by GMSH for entity " + to_string(entityTag) + "!");

			for (int elemType : entityElementTypes)
				meshElementTypes.insert(elemType);

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

				vector<Element<Dim>*> createdElements(elements.size());
				double h = -1;
				double regularity = 2;
				double sum_h = 0;

				#pragma omp parallel for reduction(max:h) reduction(min:regularity) reduction(+:sum_h)
				for (BigNumber j = 0; j < elements.size(); ++j)
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

					createdElements[j] = e;

					h          = max(h,          e->Diameter());
					regularity = min(regularity, e->Regularity());
					sum_h     += e->Diameter();
				}

				for (Element<Dim>* e : createdElements)
				{
					this->AddElement(e, false);
					_elementExternalNumbers.insert({ e->Id, e });
				}

				_h          = max(_h,          h);
				_regularity = min(_regularity, regularity);
				_average_h += sum_h;

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
	{
		if (_mshFileIsTmp)
			remove(_mshFilePath.c_str()); // if not already removed by an export
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

		// One file per fine mesh, which keeps it until its destruction (or the end of the process)
		static int nFineMeshes = 0;
		nFineMeshes++;
		string coarse_mesh_tmp_file = FileSystem::TemporaryFile("coarse.msh");
		string fine_mesh_tmp_file = FileSystem::TemporaryFile("fine" + to_string(nFineMeshes) + ".msh");

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
inline Element<2>* GMSHMesh<2>::CreateElement(int elemType, const vector<size_t>& elementNodes, size_t start, size_t elemIndex)
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
inline Element<3>* GMSHMesh<3>::CreateElement(int elemType, const vector<size_t>& elementNodes, size_t start, size_t elemIndex)
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
		// The corners by their side of the center (front: +x, right: +y, top: +z), whatever GMSH's order of the nodes:
		// corners[x + 2y + 4z]
		Vertex* nodes[8];
		double cx = 0, cy = 0, cz = 0;
		for (int i = 0; i < 8; i++)
		{
			nodes[i] = GetVertexFromGMSHTag(elementNodes[8 * elemIndex + i]);
			cx += nodes[i]->X / 8;
			cy += nodes[i]->Y / 8;
			cz += nodes[i]->Z / 8;
		}
		Vertex* corners[8];
		for (Vertex* v : nodes)
			corners[(v->X > cx) + 2 * (v->Y > cy) + 4 * (v->Z > cz)] = v;

		_parallelepipedElements[start + elemIndex] = ParallelepipedElement(0, corners[0], corners[1], corners[2], corners[4], corners[5], corners[6], corners[3], corners[7]);
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
inline void GMSHMesh<2>::AllocateContiguousMemoryForElements(int elemType, BigNumber numberOfElements)
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
inline void GMSHMesh<2>::AllocateContiguousMemoryForFaces(GMSHFaceTypes faceType, BigNumber numberOfFaces)
{
	this->Faces.reserve(numberOfFaces);
	this->_edgeFaces.reserve(numberOfFaces);
}

template <>
inline vector<size_t> GMSHMesh<2>::FaceNodes(int elemType, int tag, int& nFaceNodes)
{
	vector<size_t> edgeNodes;
	bool onlyPrimaryNodes = true;
	gmsh::model::mesh::getElementEdgeNodes(elemType, edgeNodes, tag, onlyPrimaryNodes);
	nFaceNodes = 2;
	return edgeNodes;
}

template <>
inline void GMSHMesh<2>::CreateFaces(int elemType, BigNumber& faceNumber)
{
	int nFaceNodes;
	vector<size_t> edgeNodes = FaceNodes(elemType, -1, nFaceNodes);

	AllocateContiguousMemoryForFaces(GMSHFaceTypes::GMSH_SegmentFace, edgeNodes.size() / 2);

	for (size_t j = 0; j < edgeNodes.size(); j += 2)
	{
		MeshVertex<2>* v1 = (MeshVertex<2>*)GetVertexFromGMSHTag(edgeNodes[j]);
		MeshVertex<2>* v2 = (MeshVertex<2>*)GetVertexFromGMSHTag(edgeNodes[j + 1]);
		if (!v1)
			Utils::FatalError("No GMSH node of tag " + to_string(edgeNodes[j]));
		if (!v2)
			Utils::FatalError("No GMSH node of tag " + to_string(edgeNodes[j + 1]));
		vector<MeshVertex<2>*> ends = { v1, v2 };
		OrientAlikeOnAllProcesses(ends);
		v1 = ends[0];
		v2 = ends[1];

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
inline void GMSHMesh<3>::AllocateContiguousMemoryForElements(int elemType, BigNumber numberOfElements)
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
inline void GMSHMesh<3>::AllocateContiguousMemoryForFaces(GMSHFaceTypes faceType, BigNumber numberOfFaces)
{
	this->Faces.reserve(numberOfFaces);

	if (faceType == GMSHFaceTypes::GMSH_TriangleFace)
		this->_triangularFaces.reserve(numberOfFaces);
	else if (faceType == GMSHFaceTypes::GMSH_QuadrilateralFace)
		this->_rectangularFaces.reserve(numberOfFaces);
	else
		Utils::FatalError("GMSH element type not managed.");
}

template <>
inline vector<size_t> GMSHMesh<3>::FaceNodes(int elemType, int tag, int& nFaceNodes)
{
	if (elemType != GMSHElementTypes::GMSH_Tetrahedron && elemType != GMSHElementTypes::GMSH_Hexahedron)
		Utils::FatalError("GMSH element type not managed.");
	GMSHFaceTypes faceType = elemType == GMSHElementTypes::GMSH_Tetrahedron ? GMSHFaceTypes::GMSH_TriangleFace : GMSHFaceTypes::GMSH_QuadrilateralFace;
	nFaceNodes = faceType; // (the number of vertices of the face)

	vector<size_t> faceNodes;
	bool onlyPrimaryNodes = true;
	gmsh::model::mesh::getElementFaceNodes(elemType, faceType, faceNodes, tag, onlyPrimaryNodes);
	return faceNodes;
}

template <>
inline void GMSHMesh<3>::CreateFaces(int elemType, BigNumber& faceNumber)
{
	int nFaceVertices;
	vector<size_t> faceNodes = FaceNodes(elemType, -1, nFaceVertices);
	GMSHFaceTypes faceType = (GMSHFaceTypes)nFaceVertices;

	AllocateContiguousMemoryForFaces(faceType, faceNodes.size() / nFaceVertices);

	for (size_t j = 0; j < faceNodes.size(); j += nFaceVertices)
	{
		vector<MeshVertex<3>*> vertices(nFaceVertices);
		for (int k = 0; k < nFaceVertices; k++)
		{
			MeshVertex<3>* v = (MeshVertex<3>*)GetVertexFromGMSHTag(faceNodes[j+k]);
			vertices[k] = v;
		}
		OrientAlikeOnAllProcesses(vertices);

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
			// The plane of the face: its normal axis, along which it is flat (GMSH's coordinates are exact to ~1e-12 only)
			int normalAxis;
			(Vect<3>(vertices[0], vertices[1]).cwiseAbs() + Vect<3>(vertices[0], vertices[2]).cwiseAbs()).minCoeff(&normalAxis);
			CartesianShapeOrientation orientation = normalAxis == 0 ? CartesianShapeOrientation::InYOZ : (normalAxis == 1 ? CartesianShapeOrientation::InXOZ : CartesianShapeOrientation::InXOY);

			// As RectangularFace takes them: the corner of smallest coordinates, its neighbour along the first axis of the
			// plane, its neighbour along the second one, the opposite corner (whatever the direction of the cycle)
			auto smallerCorner = [](MeshVertex<3>* a, MeshVertex<3>* b) { return a->X + a->Y + a->Z < b->X + b->Y + b->Z; };
			int o = min_element(vertices.begin(), vertices.end(), smallerCorner) - vertices.begin();
			Vertex* origin = vertices[o];
			Vertex* vertex1 = vertices[(o + 1) % 4];
			Vertex* vertex2 = vertices[(o + 3) % 4];
			Vertex* vertex3 = vertices[(o + 2) % 4];
			int firstAxis = orientation == CartesianShapeOrientation::InYOZ ? 1 : 0;
			if (abs(Vect<3>(origin, vertex1)[firstAxis]) < abs(Vect<3>(origin, vertex2)[firstAxis]))
				swap(vertex1, vertex2);

			if (neighbours.size() == 1)
				this->_rectangularFaces.emplace_back(faceNumber++, origin, vertex1, vertex2, vertex3, neighbours[0], orientation);
			else if (neighbours.size() == 2)
				this->_rectangularFaces.emplace_back(faceNumber++, origin, vertex1, vertex2, vertex3, neighbours[0], neighbours[1], orientation);
			else
				assert(false);
			face = &this->_rectangularFaces.back();

			// The side of each element the face is on, which gives its outer normal: the side of smaller coordinates if the
			// origin of the face is nearer to the back left bottom corner of the element than to its front right top one
			for (Element<3>* e : neighbours)
			{
				ParallelepipedElement* p = dynamic_cast<ParallelepipedElement*>(e);
				bool smallerSide = abs(Vect<3>(origin, p->BackLeftBottomCorner)[normalAxis]) < abs(Vect<3>(origin, p->FrontRightTopCorner)[normalAxis]);
				Face<3>*& side = normalAxis == 0 ? (smallerSide ? p->BackFace : p->FrontFace) : (normalAxis == 1 ? (smallerSide ? p->LeftFace : p->RightFace) : (smallerSide ? p->BottomFace : p->TopFace));
				side = face;
			}
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
inline Face<2>* GMSHMesh<2>::GetBoundaryFaceFromGMSHNodes(int faceType, const vector<size_t>& faceNodes, size_t& faceNodeIndex)
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
inline Face<3>* GMSHMesh<3>::GetBoundaryFaceFromGMSHNodes(int faceType, const vector<size_t>& faceNodes, size_t& faceNodeIndex)
{
	assert(faceType == GMSHElementTypes::GMSH_Triangle || faceType == GMSHElementTypes::GMSH_Quadrilateral);

	// (three vertices identify a quadrilateral too)
	MeshVertex<3>* v1 = (MeshVertex<3>*)GetVertexFromGMSHTag(faceNodes[faceNodeIndex]);
	MeshVertex<3>* v2 = (MeshVertex<3>*)GetVertexFromGMSHTag(faceNodes[faceNodeIndex + 1]);
	MeshVertex<3>* v3 = (MeshVertex<3>*)GetVertexFromGMSHTag(faceNodes[faceNodeIndex + 2]);
	faceNodeIndex += faceType == GMSHElementTypes::GMSH_Triangle ? 3 : 4;

	for (Face<3>* f : v1->Faces)
	{
		if (f->IsDomainBoundary && f->HasVertex(v2) && f->HasVertex(v3))
			return f;
	}
	assert(false && "Face not found");
	return nullptr;
}
#endif // ENABLE_3D
