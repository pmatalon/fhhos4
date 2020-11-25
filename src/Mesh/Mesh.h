#pragma once
#include "Element.h"
#include "Face.h"
#include "../Utils/ElementParallelLoop.h"
#include "../Utils/MatlabScript.h"
#include "../Utils/FileSystem.h"
using namespace std;

struct CoarseningStrategyDetails
{
	CoarseningStrategy CS;
	int nFineFacesAddedByCoarseElement = -1;
	int nFineElementsByCoarseElement = -1;
	int nFineFacesByKeptCoarseFace = -1;

	CoarseningStrategyDetails(CoarseningStrategy cs)
	{
		CS = cs;
	}
	CoarseningStrategyDetails() : CoarseningStrategyDetails(CoarseningStrategy::None) {}

	bool HasDetails()
	{
		return nFineElementsByCoarseElement != -1 && nFineFacesAddedByCoarseElement != -1 && nFineFacesByKeptCoarseFace != -1;
	}
};

template <int Dim>
class MeshVertex : public Vertex
{
public:
	vector<Element<Dim>*> Elements;
	vector<Face<Dim>*> Faces;

	mutex Mutex;

	MeshVertex(BigNumber number, double x) : Vertex(number, x) {}
	MeshVertex(BigNumber number, double x, double y) : Vertex(number, x, y) {}
	MeshVertex(BigNumber number, double x, double y, double z) : Vertex(number, x, y, z) {}
	MeshVertex(BigNumber number, const DomPoint& p) : Vertex(number, p) {}
	MeshVertex(const Vertex v) : Vertex(v) {}

	bool IsVertexOf(Element<Dim>* e)
	{
		for (Element<3>* e2 : this->Elements)
		{
			if (e2 == e)
				return true;
		}
		return false;
	}

	~MeshVertex() override
	{
		Elements.clear();
		Faces.clear();
	}
};

template <int Dim>
class Mesh
{
private:
	atomic<uint32_t> _currentElementId{ 0 };
	double _coarseningFactor = 0;
protected:
	// Vertices to keep at all levels to avoid "eroding" the geometry
	set<Vertex*> _geometricVertices;
public:
	vector<Vertex*> Vertices;
	vector<Element<Dim>*> Elements;
	vector<Face<Dim>*> Faces;
	vector<Face<Dim>*> BoundaryFaces;
	vector<Face<Dim>*> InteriorFaces;
	vector<Face<Dim>*> DirichletFaces;
	vector<Face<Dim>*> NeumannFaces;

	vector<PhysicalGroup<Dim>*> PhysicalParts;
	vector<BoundaryGroup*> BoundaryParts;

	Mesh<Dim>* CoarseMesh = nullptr;
	Mesh<Dim>* FineMesh = nullptr;
	CoarseningStrategyDetails ComesFrom;

	mutex MutexElements;
	mutex MutexFaces;

	static string MeshDirectory;
	static string CacheDirectory;

	Mesh() {}

	virtual string Description() = 0;
	virtual string FileNamePart() = 0;
	virtual string GeometryDescription() = 0;
	virtual double H() = 0;
	virtual double Regularity() = 0;

	static void SetDirectories()
	{
		MeshDirectory = FileSystem::RootPath() + "/data/meshes/";
		CacheDirectory = FileSystem::RootPath() + "/data/cache/";
	}

	double CoarseningFactor()
	{
		return _coarseningFactor;
	}

	BigNumber NewElementId()
	{
		return _currentElementId++;
	}

	void AddElement(Element<Dim>* e, bool lock = true)
	{
		if (lock)
			MutexElements.lock();

		if (e->Id == 0)
			e->Id = this->NewElementId();
		if (e->Number == -1 || e->Number == 0)
			e->Number = this->Elements.size();
		this->Elements.push_back(e);

		if (lock)
			MutexElements.unlock();
	}

	void RemoveElement(Element<Dim>* e)
	{
		e->IsDeleted = true;
	}

	void ReplaceAndDeleteElement(Element<Dim>* oldElement, Element<Dim>* newElement)
	{
		newElement->Number = oldElement->Number;
		AddElement(newElement);
		oldElement->IsDeleted = true;
	}

	void AddFace(Face<Dim>* f, bool addToBoundaryAndIteriorLists = true)
	{
		MutexFaces.lock();

		if (f->Number == -1 || f->Number == 0)
			f->Number = this->Faces.size();

		int s = this->Faces.size();
		this->Faces.push_back(f);
		assert(this->Faces.size() == s + 1);

		if (addToBoundaryAndIteriorLists)
		{
			if (f->IsDomainBoundary)
				this->BoundaryFaces.push_back(f);
			else
				this->InteriorFaces.push_back(f);
		}

		MutexFaces.unlock();
	}

	void RemoveFace(Face<Dim>* f, bool removeFromBoundaryAndIteriorLists = true)
	{
		f->IsDeleted = true;

		if (removeFromBoundaryAndIteriorLists)
		{
			// TODO: optimize
			BigNumber i = 0;
			if (f->IsDomainBoundary)
			{
				for (i = 0; i < this->BoundaryFaces.size(); i++)
				{
					if (this->BoundaryFaces[i] == f)
						break;
				}
				this->BoundaryFaces.erase(this->BoundaryFaces.begin() + i);
			}
			else
			{
				for (i = 0; i < this->InteriorFaces.size(); i++)
				{
					if (this->InteriorFaces[i] == f)
						break;
				}
				this->InteriorFaces.erase(this->InteriorFaces.begin() + i);
			}
		}
	}

	void FinalizeCreation()
	{
		vector<Element<Dim>*> tmpElements(this->Elements.begin(), this->Elements.end());
		this->Elements.clear();
		for (Element<Dim>* e : tmpElements)
		{
			if (e->IsDeleted)
				delete e;
			else
				this->Elements.push_back(e);
		}

		vector<Face<Dim>*> tmpFaces(this->Faces.begin(), this->Faces.end());
		this->Faces.clear();
		for (Face<Dim>* f : tmpFaces)
		{
			if (f->IsDeleted)
				delete f;
			else
				this->Faces.push_back(f);
		}

		BigNumber number = 0;
		for (Element<Dim>* e : this->Elements)
			e->Number = number++;
		number = 0;
		for (Face<Dim>* f : this->Faces)
			f->Number = number++;

		this->InitFaceLocalNumbering();
	}

	virtual size_t MemoryUsage()
	{
		return 0;
	}

	double SkeletonMeasure()
	{
		double measure = 0;
		for (Face<Dim>* face : this->Faces)
			measure += face->Measure();
		return measure;
	}

	Face<Dim>* ExistingFaceWithVertices(const vector<MeshVertex<Dim>*>& vertices)
	{
		for (Face<Dim>* f : vertices[0]->Faces)
		{
			bool thisFaceHasAllVertices = true;
			for (int k = 1; k < vertices.size(); k++)
			{
				if (!f->HasVertex(vertices[k]))
				{
					thisFaceHasAllVertices = false;
					break;
				}
			}
			if (thisFaceHasAllVertices)
				return f;
		}
		return nullptr;
	}

	void SetDiffusionField(DiffusionField<Dim>* diffusionField)
	{
		for (PhysicalGroup<Dim>* phyPart : PhysicalParts)
		{
			phyPart->ConstantDiffTensor = diffusionField->ConstantDiffTensor(phyPart);
		}
	}

	void SetBoundaryConditions(BoundaryConditions* bc)
	{
		for (BoundaryGroup* boundary : BoundaryParts)
		{
			boundary->Condition = bc->GetBoundaryConditionType(boundary);
			boundary->ConditionFunction = boundary->Condition == BoundaryConditionType::Neumann ? bc->NeumannFunction : bc->DirichletFunction;
		}
		
		FillDirichletAndNeumannFaceLists();
	}

	void ExportToMatlab(string outputDirectory)
	{
		string faceFilePath = outputDirectory + "/mesh_faces_" + to_string(Dim) + "D_" + FileNamePart() + ".dat";
		ExportFacesToMatlab(faceFilePath);
		string elemFilePath = outputDirectory + "/mesh_elements_" + to_string(Dim) + "D_" + FileNamePart() + ".dat";
		ExportElementCentersToMatlab(elemFilePath);
	}
	virtual void ExportFacesToMatlab(string filePath)
	{
		FILE* file = fopen(filePath.c_str(), "w");
		fprintf(file, "Number x1    y1    x2    y2 IsDomainBoundary IsRemovedOnCoarserGrid\n");
		for (Face<Dim>* f : this->Faces)
			f->ExportFaceToMatlab(file);
		fclose(file);
		cout << "Faces exported to         " << filePath << endl;
	}
	void ExportElementCentersToMatlab(string filePath)
	{
		MatlabScript script(filePath);
		for (auto e : this->Elements)
			script.PlotText(e->Center(), to_string(e->Number), "r");
	}

	virtual void ExportToGMSH(FunctionalBasis<Dim>* basis, const Vector &coeffs, const string& outputFilePathPrefix, const string& suffix)
	{
		Utils::Warning("Impossible to export the solution to GMSH because this mesh does not come from GMSH.");
	}

	virtual void CoarsenMesh(CoarseningStrategy strategy, int coarseningFactor)
	{
		if (Utils::IsRefinementStrategy(strategy))
			return;
		Utils::FatalError("Unmanaged coarsening strategy");
	}
	virtual void RefineMesh(CoarseningStrategy strategy)
	{
		Utils::FatalError("Unmanaged refinement strategy");
	}

protected:
	void InitFaceLocalNumbering()
	{
		ElementParallelLoop<Dim> parallelLoop(this->Elements);
		parallelLoop.Execute([](Element<Dim>* e)
			{
				e->InitFaceLocalNumbering();
			});
	}

	void InitFinerElementsLocalNumbering()
	{
		ElementParallelLoop<Dim> parallelLoop(this->Elements);
		parallelLoop.Execute([](Element<Dim>* e)
			{
				e->InitFinerElementsLocalNumbering();
			});
	}

	virtual void InitializeCoarsening(Mesh<Dim>* coarseMesh)
	{
		this->CoarseMesh = coarseMesh;
		coarseMesh->FineMesh = this;

		coarseMesh->PhysicalParts = this->PhysicalParts;
		coarseMesh->BoundaryParts = this->BoundaryParts;

		coarseMesh->_geometricVertices = this->_geometricVertices;
	}

	virtual void InitializeRefinement(Mesh<Dim>* fineMesh)
	{
		this->FineMesh = fineMesh;
		fineMesh->CoarseMesh = this;

		fineMesh->PhysicalParts = this->PhysicalParts;
		fineMesh->BoundaryParts = this->BoundaryParts;

		fineMesh->_geometricVertices = this->_geometricVertices;
	}

	virtual void FinalizeCoarsening()
	{
		CoarseMesh->FinalizeCreation();
		CoarseMesh->FillBoundaryAndInteriorFaceLists();
		CoarseMesh->FillDirichletAndNeumannFaceLists();
		CoarseMesh->InitFinerElementsLocalNumbering();
		CoarseMesh->_coarseningFactor = CoarseMesh->H() / this->H();
	}

	virtual void FinalizeRefinement()
	{
		this->InitFinerElementsLocalNumbering();
		this->_coarseningFactor = this->H() / FineMesh->H();
	}

	void FillBoundaryAndInteriorFaceLists()
	{
		bool fillBoundaryFaces = this->BoundaryFaces.empty();
		bool fillInteriorFaces = this->InteriorFaces.empty();
		if (!fillBoundaryFaces && !fillInteriorFaces)
			return;

		if (fillInteriorFaces)
			this->InteriorFaces.reserve(this->Faces.size());

		for (Face<Dim>* f : this->Faces)
		{
			if (fillBoundaryFaces && f->IsDomainBoundary)
				this->BoundaryFaces.push_back(f);
			else if (fillInteriorFaces && !f->IsDomainBoundary)
				this->InteriorFaces.push_back(f);
		}
	}

	void FillDirichletAndNeumannFaceLists()
	{
		if (!this->DirichletFaces.empty() || !this->NeumannFaces.empty())
			return;

		assert(!this->BoundaryFaces.empty());

		for (Face<Dim>* f : this->BoundaryFaces)
		{
			if (f->HasDirichletBC())
				this->DirichletFaces.push_back(f);
			else if (f->HasNeumannBC())
				this->NeumannFaces.push_back(f);
			else
				assert(false);
		}

		if (CoarseMesh)
			CoarseMesh->FillDirichletAndNeumannFaceLists();
	}


	void LinkFacesToCoarseFaces()
	{
		// Continue linking
		for (Face<Dim>* fineFace : this->Faces)
		{
			// Boundary face
			if (fineFace->IsDomainBoundary)
			{
				fineFace->IsRemovedOnCoarserGrid = false;
				for (Face<Dim>* coarseFace : fineFace->Element1->CoarserElement->Faces)
				{
					if (coarseFace->IsDomainBoundary && coarseFace->Contains(fineFace->Center()))
					{
						fineFace->CoarseFace = coarseFace;
						coarseFace->FinerFaces.push_back(fineFace);
						fineFace->BoundaryPart = coarseFace->BoundaryPart;
						break;
					}
				}
				// If no coarse face has been found, it may be because the geometry is curved and the refinement yields a non-nested mesh.
				// In that case, we take the closest one w.r.t. the centers
				if (!fineFace->CoarseFace)
				{
					Face<Dim>* closestCoarseFace = fineFace->ClosestFaceAmongst(fineFace->Element1->CoarserElement->Faces, true);
					assert(closestCoarseFace && "A coarse face should have been found.");

					fineFace->CoarseFace = closestCoarseFace;
					closestCoarseFace->FinerFaces.push_back(fineFace);
					fineFace->BoundaryPart = closestCoarseFace->BoundaryPart;
				}
			}
			// Interior face
			else
			{
				Element<Dim>* coarseElement1 = fineFace->Element1->CoarserElement;
				Element<Dim>* coarseElement2 = fineFace->Element2->CoarserElement;
				fineFace->IsRemovedOnCoarserGrid = coarseElement1 == coarseElement2;
				if (!fineFace->IsRemovedOnCoarserGrid)
				{
					vector<Face<Dim>*> coarseFaces = coarseElement1->InterfaceWith(coarseElement2);
					assert(!coarseFaces.empty());
					Face<Dim>* closestCoarseFace = fineFace->ClosestFaceAmongst(coarseFaces, false);
					assert(closestCoarseFace && "A coarse face should have been found.");

					fineFace->CoarseFace = closestCoarseFace;
					closestCoarseFace->FinerFaces.push_back(fineFace);
					fineFace->BoundaryPart = closestCoarseFace->BoundaryPart;
				}
				else
					coarseElement1->FinerFacesRemoved.push_back(fineFace);
			}
		}
	}

public:
	virtual void Serialize(ostream& os) const
	{
		for (auto element : this->Elements)
			os << *element << endl;
		for (auto face : this->Faces)
			os << *face << endl;
	}

	friend ostream& operator<<(ostream& os, const Mesh<Dim>& s)
	{
		s.Serialize(os);
		return os;
	}

	// Performs sanity checks to verify that the mesh correctly built
	virtual void SanityCheck()
	{
		for (auto e : this->Elements)
			e->UnitTests();

		for (auto f : this->Faces)
			assert(!f->IsDeleted);

		for (auto e : this->Elements)
		{
			assert(!e->IsDeleted);
			for (auto f : e->Faces)
			{
				assert(!f->IsDeleted);

				auto n = e->OuterNormalVector(f);

				auto neighbour = f->GetNeighbour(e);
				if (neighbour != nullptr)
				{
					auto n2 = neighbour->OuterNormalVector(f);
					if (abs(n.dot(n2) + 1) >= Utils::NumericalZero)
					{
						// Analysis of the problem
						cout << "Problem with normal vectors at " << *f << endl;
						cout << "One of them is not pointing outwards." << endl;
						double d = abs(n.dot(n2) + 1);
						auto c = e->Center();
						cout << "%------- Matlab script to plot element " << e->Number << ":" << endl;
						e->ExportToMatlab("r", true);
						cout << endl;
						cout << "%------- Matlab script to plot element " << neighbour->Number << ":" << endl;
						neighbour->ExportToMatlab("b", true);
						cout << endl;
						MatlabScript s;
						s.PlotSegment(f->Vertices()[0], f->Vertices()[1], "k-", 2);

						n = e->OuterNormalVector(f);
						n2 = neighbour->OuterNormalVector(f);
						e->UnitTests();
						Utils::FatalError("Problem with normal vectors: one of them is not pointing outwards.");
					}
				}
				else
					assert(e->IsOnBoundary() && "This element has no neighbour on the other side of this face, althoug it's not supposed to be on the boundary.");
			}
		}

		for (auto f : this->BoundaryFaces)
		{
			assert(f->IsDomainBoundary && "This face is in the BoundaryFaces list but has not the flag IsDomainBoundary.");
			if (this->ComesFrom.CS != CoarseningStrategy::FaceCoarsening)
			{
				assert(f->Element1 && "This face is on the boundary but has no Element1.");
				assert(!f->Element2 && "This face is on the boundary but connects two elements.");
			}
		}

		for (auto f : this->InteriorFaces)
		{
			assert(!f->IsDomainBoundary);
			if (this->ComesFrom.CS != CoarseningStrategy::FaceCoarsening)
			{
				assert(f->Element1 && "This face is in the interior but has no Element1.");
				assert(f->Element2 && "This face is in the interior but has no Element2.");
			}
		}

		assert(this->BoundaryFaces.size() + this->InteriorFaces.size() == this->Faces.size() && "Boundary faces + interior faces does NOT add up to the total number of faces.");
				
		if (CoarseMesh)
		{
			CoarseMesh->SanityCheck();

			// Links between coarse and fine meshes

			assert(this->PhysicalParts.size() == CoarseMesh->PhysicalParts.size() && "Discrepency in the PhysicalParts of the fine and coarse meshes.");
			assert(this->BoundaryParts.size() == CoarseMesh->BoundaryParts.size() && "Discrepency in the BoundaryParts of the fine and coarse meshes.");

			for (int i = 0; i < this->PhysicalParts.size(); i++)
				assert(this->PhysicalParts[i] == CoarseMesh->PhysicalParts[i] && "The PhysicalParts in fine and coarse meshes must have the same addresses.");
			for (int i = 0; i < this->BoundaryParts.size(); i++)
				assert(this->BoundaryParts[i] == CoarseMesh->BoundaryParts[i] && "The BoundaryParts in fine and coarse meshes must have the same addresses.");

			assert(this->NeumannFaces.empty() == CoarseMesh->NeumannFaces.empty() && "The fine or coarse mesh has NeumannFaces while the other has none.");

			if (CoarseMesh->ComesFrom.CS != CoarseningStrategy::FaceCoarsening)
			{
				for (Element<Dim>* fe : this->Elements)
				{
					assert(fe->CoarserElement != nullptr && "This fine element has no coarse element.");

					if (fe->PhysicalPart)
						assert(fe->PhysicalPart == fe->CoarserElement->PhysicalPart && "Associated fine and coarse element have different PhysicalParts.");
					else
						assert(!fe->CoarserElement->PhysicalPart && "This fine element has no PhysicalPart while its associated coarse element has one.");

					bool feIsReferenced = false;
					for (Element<Dim>* e : fe->CoarserElement->FinerElements)
					{
						assert(e->CoarserElement == fe->CoarserElement && "Discrepency in the link between fine and coarse elements.");
						if (e == fe)
						{
							feIsReferenced = true;
							break;
						}
					}
					assert(feIsReferenced && "This fine element is not referenced in the list FinerElements of its associated coarse element.");
				}

				BigNumber totalFinerElements = 0;
				for (Element<Dim>* ce : CoarseMesh->Elements)
				{
					totalFinerElements += ce->FinerElements.size();

					if (Utils::BuildsNestedMeshHierarchy(CoarseMesh->ComesFrom.CS))
						assert(ce->FinerElements.size() > 0 && "This coarse element doesn't have any finer element.");

					for (Face<Dim>* ff : ce->FinerFacesRemoved)
						assert(ff->IsRemovedOnCoarserGrid && "This face is supposed to be removed inside this coarse element, but the flag IsRemovedOnCoarserGrid is not set.");
				}

				assert(totalFinerElements == this->Elements.size() && "Discrepency between the number of elements in the fine mesh and the number of fine elements associated to coarse elements.");
			}

			if (CoarseMesh->ComesFrom.CS != CoarseningStrategy::IndependentRemeshing)
			{
				for (Face<Dim>* cf : CoarseMesh->Faces)
				{
					assert(cf->FinerFaces.size() > 0 && "This coarse face doesn't have any finer face.");
					for (Face<Dim>* ff : cf->FinerFaces)
					{
						assert(!ff->IsRemovedOnCoarserGrid && "This fine face is referenced in the list FinerFaces of a coarse face, but is flagged IsRemovedOnCoarserGrid.");
						assert(ff->CoarseFace == cf && "Discrepency in the link between fine and coarse faces.");
					}
				}

				for (Face<Dim>* ff : this->Faces)
				{
					if (ff->IsRemovedOnCoarserGrid && Utils::BuildsNestedMeshHierarchy(CoarseMesh->ComesFrom.CS))
					{
						Element<Dim>* ce1 = ff->Element1->CoarserElement;
						Element<Dim>* ce2 = ff->Element2->CoarserElement;
						assert(ce1 == ce2 && "This fine face is flagged IsRemovedOnCoarserGrid but Element1->CoarserElement != Element2->CoarserElement.");

						bool ffIsReferenced = false;
						for (Face<Dim>* removedFace : ce1->FinerFacesRemoved)
						{
							if (ff == removedFace)
							{
								ffIsReferenced = true;
								break;
							}
						}
						assert(ffIsReferenced && "This fine face is flagged IsRemovedOnCoarserGrid but is not referenced in the list FinerFacesRemoved of the coarse element it's supposed to be inside.");
					}

					if (ff->IsDomainBoundary)
						assert(ff->CoarseFace && "This boundary face has no coarse face.");
				}

				for (Face<Dim>* ff : this->Faces)
				{
					if (ff->IsRemovedOnCoarserGrid)
						continue;
					if (ff->BoundaryPart)
					{
						assert(ff->CoarseFace->BoundaryPart && "This fine face has a BoundaryPart but its coarsened one has none.");
						assert(ff->CoarseFace->BoundaryPart == ff->BoundaryPart && "This face and its coarsened one have different BoundaryParts.");
					}
					else
						assert(!ff->CoarseFace->BoundaryPart && "This fine face has no BoundaryPart but its coarsened one has one.");
				}
			}


			struct ChunkResult { double total = 0; };
			ParallelLoop<Element<Dim>*, ChunkResult> parallelLoopFine(this->Elements);
			parallelLoopFine.Execute([](Element<Dim>* e, ParallelChunk<ChunkResult>* chunk)
				{
					chunk->Results.total += e->Measure();
				});
			double fineMeshTotalMeasure = 0;
			parallelLoopFine.AggregateChunkResults([&fineMeshTotalMeasure](ChunkResult& chunk)
				{
					fineMeshTotalMeasure += chunk.total;
				});
			
			ParallelLoop<Element<Dim>*, ChunkResult> parallelLoopCoarse(CoarseMesh->Elements);
			parallelLoopCoarse.Execute([](Element<Dim>* e, ParallelChunk<ChunkResult>* chunk)
				{
					chunk->Results.total += e->Measure();
				});
			double coarseMeshTotalMeasure = 0;
			parallelLoopCoarse.AggregateChunkResults([&coarseMeshTotalMeasure](ChunkResult& chunk)
				{
					coarseMeshTotalMeasure += chunk.total;
				});

			if (abs(coarseMeshTotalMeasure - fineMeshTotalMeasure) > Utils::NumericalZero)
				Utils::Error("Fine and coarse meshes should have the same total measure.");

			
			Mesh<Dim>* meshToGetInfo = nullptr;
			if (Utils::IsRefinementStrategy(this->ComesFrom.CS))
				meshToGetInfo = this;
			else if (CoarseMesh->ComesFrom.CS == CoarseningStrategy::StandardCoarsening || CoarseMesh->ComesFrom.CS == CoarseningStrategy::AgglomerationCoarsening)
				meshToGetInfo = CoarseMesh;

			if (meshToGetInfo != nullptr && meshToGetInfo->ComesFrom.HasDetails())
			{
				int nFineFacesAddedByCoarseElement = meshToGetInfo->ComesFrom.nFineFacesAddedByCoarseElement;
				int nFineElementsByCoarseElement = meshToGetInfo->ComesFrom.nFineElementsByCoarseElement;
				int nFineFacesByKeptCoarseFace = meshToGetInfo->ComesFrom.nFineFacesByKeptCoarseFace;

				for (Element<Dim>* ce : CoarseMesh->Elements)
				{
					assert(ce->FinerElements.size() == nFineElementsByCoarseElement && "This coarse element does not have the expected number of fine elements.");
					assert(ce->FinerFacesRemoved.size() == nFineFacesAddedByCoarseElement && "This coarse element does not have the expected number of FinerFacesRemoved.");
				}
				BigNumber finerFacesAdded = nFineFacesAddedByCoarseElement * CoarseMesh->Elements.size();
				assert(this->Elements.size() == nFineElementsByCoarseElement * CoarseMesh->Elements.size());
				assert(this->Faces.size() == nFineFacesByKeptCoarseFace * CoarseMesh->Faces.size() + finerFacesAdded);
			}
		}
	}

	Mesh<Dim>* RefineNTimes(int nRefinements, CoarseningStrategy strategy)
	{
		Mesh<Dim>* mesh = this;
		for (int i = 0; i < nRefinements; i++)
		{
			mesh->RefineMesh(strategy);
			mesh = mesh->FineMesh;
		}
		return mesh;
	}

	Mesh<Dim>* RefineUntilNElements(BigNumber nElements, CoarseningStrategy strategy)
	{
		cout << "Building fine mesh by successive refinements of the coarse mesh" << endl;
		Mesh<Dim>* mesh = this;
		while (mesh->Elements.size() < nElements)
		{
			mesh->RefineMesh(strategy);
			mesh = mesh->FineMesh;
		}
		return mesh;
	}

	void DeleteCoarseMeshes()
	{
		delete this->CoarseMesh;
		this->CoarseMesh = nullptr;

		for (Element<Dim>* e : this->Elements)
			e->CoarserElement = nullptr;
		for (Face<Dim>* f : this->Faces)
		{
			f->CoarseFace = nullptr;
			f->IsRemovedOnCoarserGrid = false;
		}

		this->ComesFrom = CoarseningStrategyDetails(CoarseningStrategy::None);
	}

	virtual void SetOverlappingFineElements()
	{
		assert(false);
	}

	virtual void SetOverlappingFineElementsSubTriangles()
	{
		assert(false);
	}

	void DeleteOverlappingFineElementsInformation()
	{
		ElementParallelLoop<Dim> parallelLoop(this->Elements);
		parallelLoop.Execute([](Element<Dim>* e)
			{
				// Deletion of the intersections
				for (auto it = e->OverlappingFineElements.begin(); it != e->OverlappingFineElements.end(); it++)
				{
					vector<PhysicalShape<Dim>*> intersectionCoarseFine = it->second;
					for (PhysicalShape<Dim>* intersection : intersectionCoarseFine)
						delete intersection;
				}
				e->OverlappingFineElements.clear();
			});
	}

	void ExportToMatlab2(string filePath = "")
	{
		//vector<string> colors = { "r", "b", "m", "g", "y", "c" };
		RotatingList<string> colors({ "r", "b", "m", "k", "g", "y", "c" });
		MatlabScript s(filePath);
		s.OpenFigure();
		s.Add("axis(axes, 'equal');");
		s.Add("axes.Visible = false;");

		map<PhysicalGroup<Dim>*, string> phyColors;
		for (PhysicalGroup<Dim>* pp : this->PhysicalParts)
			phyColors.insert({ pp, colors.GetAndMoveNext() });

		for (Element<Dim>* e : this->Elements)
		{
			if (!e->PhysicalPart || e->PhysicalPart->Name.compare("domain") == 0)
			{
				s.PlotPolygonEdges(e->Shape()->Vertices(), "k");
			}
			else
			{
				string color = phyColors.at(e->PhysicalPart);//colors[e->PhysicalPart->Id % colors.size()];
				s.PlotPolygon(e->Shape()->Vertices(), color);
			}
			//s.PlotText(e->Center(), to_string(e->Number));
		}

		for (Face<Dim>* f : this->Faces)
		{
			bool isAtPhysicalBoundary = !f->IsDomainBoundary && !f->Element1->IsInSamePhysicalPartAs(f->Element2);
			if (f->IsDomainBoundary || isAtPhysicalBoundary)
				s.PlotSegment(f->Vertices()[0], f->Vertices()[1], "k", 1);
		}
	}

	void PlotClustersForApproxL2(bool useOverlappingElements)
	{
		//vector<string> colors = { "r", "b", "g", "k", "m", "y", "c" };
		RotatingList<string> colors({ "r", "b", "g", "k", "m", "y", "c" });
		MatlabScript s("/mnt/c/Users/pierr/Desktop/approx_L2_1.m");
		s.OpenFigure();
		s.Add("axis(axes, 'equal');");
		s.Add("axes.Visible = false;");

		map<Element<Dim>*, string> elementColors;

		for (Element<Dim>* coarse : this->Elements)
		{
			s.Comment("---------------------------------------");
			int nTestedColors = 0;
			string color = colors.GetAndMoveNext();
			//string color = colors[coarse->Number % colors.size()];
			bool colorDifferentFromNeighbours = true;
			do
			{
				colorDifferentFromNeighbours = true;
				for (Element<Dim>* neighbour : coarse->Neighbours(false))
				{
					auto it = elementColors.find(neighbour);
					if (it != elementColors.end())
					{
						string neighbourColor = it->second;
						if (color.compare(neighbourColor) == 0)
						{
							colorDifferentFromNeighbours = false;
							nTestedColors++;
							color = colors.GetAndMoveNext();
							break;
						}
					}
				}
			} while (!colorDifferentFromNeighbours && nTestedColors < colors.Size());

			elementColors.insert({ coarse, color });

			s.Comment("----- coarse " + to_string(coarse->Number));
			if (useOverlappingElements)
			{
				for (auto it = coarse->OverlappingFineElements.begin(); it != coarse->OverlappingFineElements.end(); it++)
				{
					Element<Dim>* fe = it->first;
					vector<PhysicalShape<Dim>*> subShapesInCoarse = it->second;
					s.Comment("----- fine " + to_string(fe->Number));
					for (auto ss : subShapesInCoarse)
						s.PlotPolygon(ss->Vertices(), color, "-");
				}
			}
			else
			{
				for (auto fe : coarse->FinerElements)
				{
					s.Comment("----- fine " + to_string(fe->Number));
					s.PlotPolygon(fe->Shape()->Vertices(), color, "-");
				}
			}
		}
		s.Comment("---------------------------------------------");
		s.Comment("--------------------COARSE-------------------");
		s.Comment("---------------------------------------------");
		for (Element<Dim>* coarse : this->Elements)
		{
			s.Comment("---------------------------------------");
			s.PlotPolygonEdges(coarse->Shape()->Vertices(), "k", 3);
			s.PlotText(coarse->Center(), to_string(coarse->Number));
		}
	}

	virtual ~Mesh() 
	{
		if (CoarseMesh)
			delete CoarseMesh;

		for (size_t i = 0; i < this->Vertices.size(); ++i)
			delete this->Vertices[i];
		this->Vertices.clear();

		if (!FineMesh)
		{
			for (size_t i = 0; i < this->PhysicalParts.size(); ++i)
				delete this->PhysicalParts[i];
			for (size_t i = 0; i < this->BoundaryParts.size(); ++i)
				delete this->BoundaryParts[i];
		}
	}
};

template <int Dim>
string Mesh<Dim>::MeshDirectory = "";
template <int Dim>
string Mesh<Dim>::CacheDirectory = "";