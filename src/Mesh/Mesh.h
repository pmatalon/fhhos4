#pragma once
#include <vector>
#include "Element.h"
#include "Face.h"
#include "../Utils/ParallelLoop.h"
#include "../Utils/MatlabScript.h"
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

	MeshVertex(BigNumber number, double x) : Vertex(number, x) {}
	MeshVertex(BigNumber number, double x, double y) : Vertex(number, x, y) {}
	MeshVertex(BigNumber number, double x, double y, double z) : Vertex(number, x, y, z) {}
	MeshVertex(BigNumber number, DomPoint p) : Vertex(number, p) {}
	MeshVertex(const Vertex v) : Vertex(v) {}

	~MeshVertex() override
	{
		Elements.clear();
		Faces.clear();
	}
};

template <int Dim>
class Mesh
{
public:
	vector<Vertex*> Vertices;
	vector<Element<Dim>*> Elements;
	vector<Face<Dim>*> Faces;
	vector<Face<Dim>*> BoundaryFaces;
	vector<Face<Dim>*> InteriorFaces;
	vector<Face<Dim>*> DirichletFaces;
	vector<Face<Dim>*> NeumannFaces;

	DiffusionPartition<Dim>* _diffusionPartition = nullptr;
	BoundaryConditions* _boundaryConditions = nullptr;

	Mesh<Dim>* CoarseMesh = nullptr;
	Mesh<Dim>* FineMesh = nullptr;
	CoarseningStrategyDetails ComesFrom;

	static string MeshDirectory;

	Mesh() {}

	virtual string Description() = 0;
	virtual string FileNamePart() = 0;
	virtual double H() = 0;
	virtual double Regularity() = 0;

	virtual void CoarsenMesh(CoarseningStrategy strategy)
	{
		Utils::FatalError("Unmanaged coarsening strategy");
	}
	virtual void RefineMesh(CoarseningStrategy strategy)
	{
		Utils::FatalError("Unmanaged refinement strategy");
	}

	void AddFace(Face<Dim>* f, bool addToBoundaryAndIteriorLists = true)
	{
		if (f->Number == -1)
			f->Number = this->Faces.size();

		this->Faces.push_back(f);

		if (addToBoundaryAndIteriorLists)
		{
			if (f->IsDomainBoundary)
				this->BoundaryFaces.push_back(f);
			else
				this->InteriorFaces.push_back(f);
		}
	}

	void RemoveElement(Element<Dim>* e)
	{
		BigNumber iToRemove = e->Number;
		assert(this->Elements[iToRemove] == e);
		this->Elements.erase(this->Elements.begin() + e->Number);
		for (int i = iToRemove; i < this->Elements.size(); i++)
			this->Elements[i]->Number = i;
	}

	void RemoveFace(Face<Dim>* f, bool removeFromBoundaryAndIteriorLists = true)
	{
		BigNumber iToRemove = f->Number;
		assert(this->Faces[iToRemove] == f);
		this->Faces.erase(this->Faces.begin() + f->Number);
		for (int i = iToRemove; i < this->Faces.size(); i++)
			this->Faces[i]->Number = i;

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

	void FillBoundaryAndInteriorFaceLists()
	{
		assert((this->BoundaryFaces.empty() && this->InteriorFaces.empty()) || (!this->BoundaryFaces.empty() && !this->InteriorFaces.empty()));

		if (!this->BoundaryFaces.empty() || !this->InteriorFaces.empty())
			return;

		this->InteriorFaces.reserve(this->Faces.size());

		for (Face<Dim>* f : this->Faces)
		{
			if (f->IsDomainBoundary)
				this->BoundaryFaces.push_back(f);
			else
				this->InteriorFaces.push_back(f);
		}
	}

	double SkeletonMeasure()
	{
		double measure = 0;
		for (Face<Dim>* face : this->Faces)
			measure += face->Measure();
		return measure;
	}

	Face<Dim>* ExistingFaceWithVertices(vector<MeshVertex<Dim>*> vertices)
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

	void SetDiffusionCoefficient(DiffusionPartition<Dim>* diffusionPartition)
	{
		this->_diffusionPartition = diffusionPartition;

		ParallelLoop<Element<Dim>*, EmptyResultChunk> parallelLoop(this->Elements);
		parallelLoop.Execute([&diffusionPartition](Element<Dim>* e, ParallelChunk<EmptyResultChunk>* chunk)
			{
				e->SetDiffusionCoefficient(diffusionPartition); // For DG
				e->SetDiffusionTensor(diffusionPartition);
			});

		if (this->CoarseMesh)
			this->CoarseMesh->SetDiffusionCoefficient(diffusionPartition);
	}

	void SetBoundaryConditions(BoundaryConditions* bc)
	{
		if (!this->DirichletFaces.empty() || !this->NeumannFaces.empty())
			return;

		this->_boundaryConditions = bc;

		ParallelLoop<Face<Dim>*, EmptyResultChunk> parallelLoop(this->BoundaryFaces);
		parallelLoop.Execute([bc](Face<Dim>* f, ParallelChunk<EmptyResultChunk>* chunk)
			{
				f->SetBoundaryConditions(bc);
			});
		
		for (auto f : this->BoundaryFaces)
		{
			if (f->HasDirichletBC())
				this->DirichletFaces.push_back(f);
			else if (f->HasNeumannBC())
				this->NeumannFaces.push_back(f);
		}

		if (this->CoarseMesh)
			this->CoarseMesh->SetBoundaryConditions(bc);
	}

	void ExportFacesToMatlab(string outputDirectory, bool dummy)
	{
		string filePath = outputDirectory + "/faces" + to_string(Dim) + "D_" + FileNamePart() + ".dat";
		ExportFacesToMatlab(filePath);
	}
	virtual void ExportFacesToMatlab(string filePath)
	{
		FILE* file = fopen(filePath.c_str(), "w");
		fprintf(file, "Number x1    y1    x2    y2 IsDomainBoundary IsRemovedOnCoarserGrid\n");
		for (Face<Dim>* f : this->Faces)
			f->ExportFaceToMatlab(file);
		fclose(file);
		cout << "Faces exported to \t" << filePath << endl;
	}
	void ExportElementCentersToMatlab(string filePath)
	{
		MatlabScript script(filePath);
		for (auto e : this->Elements)
			script.PlotText(e->Center(), to_string(e->Number), "r");
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
						break;
					}
				}
				// If no coarse face has been found, it may be because the geometry is curved and the refinement yields a non-nested mesh.
				// In that case, we take the closest one w.r.t. the centers
				if (!fineFace->CoarseFace)
				{
					double smallestDistance = -1;
					Face<Dim>* closestCoarseFace = nullptr;
					for (Face<Dim>* coarseFace : fineFace->Element1->CoarserElement->Faces)
					{
						if (coarseFace->IsDomainBoundary)
						{
							double distance = Vect<Dim>(fineFace->Center(), coarseFace->Center()).norm();
							if (!closestCoarseFace || distance < smallestDistance)
							{
								smallestDistance = distance;
								closestCoarseFace = coarseFace;
							}
						}
					}
					if (closestCoarseFace)
					{
						fineFace->CoarseFace = closestCoarseFace;
						closestCoarseFace->FinerFaces.push_back(fineFace);
					}
					else
						assert(false && "A coarse face should have been found.");
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
					// TODO: change CommonFaceWith so it returns a vector of faces (local refinement)
					Face<Dim>* coarseFace = coarseElement1->CommonFaceWith(coarseElement2);
					assert(coarseFace != nullptr);
					fineFace->CoarseFace = coarseFace;
					coarseFace->FinerFaces.push_back(fineFace);
				}
				else
					coarseElement1->FinerFacesRemoved.push_back(fineFace);
			}
		}
	}

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

	virtual void SanityCheck()
	{
		for (auto e : this->Elements)
			e->UnitTests();

		for (auto e : this->Elements)
		{
			for (auto f : e->Faces)
			{
				auto n = e->OuterNormalVector(f);

				auto neighbour = f->GetNeighbour(e);
				if (neighbour != nullptr)
				{
					auto n2 = neighbour->OuterNormalVector(f);
					if (abs(n.dot(n2) + 1) >= 1e-15)
					{
						// Analysis of the problem
						cout << "Problem with normal vectors at " << *f << endl;
						cout << "One of them is not pointing outwards." << endl;
						double d = abs(n.dot(n2) + 1);
						auto c = e->Center();
						cout << "Matlab script to plot element " << e->Number << ":" << endl;
						e->ExportToMatlab();
						cout << endl;
						cout << "Matlab script to plot element " << neighbour->Number << ":" << endl;
						neighbour->ExportToMatlab();
						cout << endl;

						n = e->OuterNormalVector(f);
						n2 = neighbour->OuterNormalVector(f);
						e->UnitTests();
						assert(false && "Problem with normal vectors: one of them is not pointing outwards.");
					}
					assert(abs(n.dot(n2) + 1) < 1e-15);
				}
				else
					assert(e->IsOnBoundary());
			}
		}

		for (auto f : this->BoundaryFaces)
		{
			assert(f->IsDomainBoundary);
			if (this->ComesFrom.CS != CoarseningStrategy::FaceCoarsening)
			{
				assert(f->Element1);
				assert(!f->Element2);
			}
		}

		for (auto f : this->InteriorFaces)
		{
			assert(!f->IsDomainBoundary);
			if (this->ComesFrom.CS != CoarseningStrategy::FaceCoarsening)
			{
				assert(f->Element1);
				assert(f->Element2);
			}
		}

		assert(this->BoundaryFaces.size() + this->InteriorFaces.size() == this->Faces.size());
		
		if (CoarseMesh)
		{
			CoarseMesh->SanityCheck();

			if (CoarseMesh->ComesFrom.CS != CoarseningStrategy::FaceCoarsening)
			{
				for (Element<Dim>* fe : this->Elements)
				{
					assert(fe->CoarserElement != nullptr);

					bool feIsReferenced = false;
					for (Element<Dim>* e : fe->CoarserElement->FinerElements)
					{
						assert(e->CoarserElement == fe->CoarserElement);
						if (e == fe)
						{
							feIsReferenced = true;
							break;
						}
					}
					assert(feIsReferenced);
				}

				for (Element<Dim>* ce : CoarseMesh->Elements)
				{
					assert(ce->FinerElements.size() > 0);
					for (Face<Dim>* ff : ce->FinerFacesRemoved)
						assert(ff->IsRemovedOnCoarserGrid);
				}
			}

			if (CoarseMesh->ComesFrom.CS != CoarseningStrategy::IndependentRemeshing)
			{
				for (Face<Dim>* cf : CoarseMesh->Faces)
				{
					assert(cf->FinerFaces.size() > 0);
					for (Face<Dim>* ff : cf->FinerFaces)
					{
						assert(!ff->IsRemovedOnCoarserGrid);
						assert(ff->CoarseFace == cf);
					}
				}

				for (Face<Dim>* ff : this->Faces)
				{
					if (ff->IsRemovedOnCoarserGrid)
					{
						Element<Dim>* ce1 = ff->Element1->CoarserElement;
						Element<Dim>* ce2 = ff->Element2->CoarserElement;
						assert(ce1 == ce2);

						bool ffIsReferenced = false;
						for (Face<Dim>* removedFace : ce1->FinerFacesRemoved)
						{
							if (ff == removedFace)
							{
								ffIsReferenced = true;
								break;
							}
						}
						assert(ffIsReferenced);
					}
				}
			}

			Mesh<Dim>* meshToGetInfo = nullptr;
			if (this->ComesFrom.CS == CoarseningStrategy::SplittingRefinement || this->ComesFrom.CS == CoarseningStrategy::BeyRefinement)
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
					assert(ce->FinerElements.size() == nFineElementsByCoarseElement);
					assert(ce->FinerFacesRemoved.size() == nFineFacesAddedByCoarseElement);
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

	virtual ~Mesh() 
	{
		for (size_t i = 0; i < this->Elements.size(); ++i)
			delete this->Elements[i];
		this->Elements.clear();

		for (size_t i = 0; i < this->Faces.size(); ++i)
			delete this->Faces[i];
		this->Faces.clear();

		if (CoarseMesh)
			delete CoarseMesh;

		for (size_t i = 0; i < this->Vertices.size(); ++i)
			delete this->Vertices[i];
		this->Vertices.clear();
	}
};

template <int Dim>
string Mesh<Dim>::MeshDirectory = "";