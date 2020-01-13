#pragma once
#include <vector>
#include "Element.h"
#include "Face.h"
#include "../Utils/ParallelLoop.h"
using namespace std;

enum class CoarseningStrategy : unsigned
{
	Standard = 0,
	Agglomeration,
	StructuredRefinement
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

	static string MeshDirectory;

	Mesh() {}

	virtual string Description() = 0;
	virtual string FileNamePart() = 0;
	virtual double H() = 0;
	virtual void CoarsenMesh(CoarseningStrategy strategy) = 0;

	void AddFace(Face<Dim>* f)
	{
		this->Faces.push_back(f);
		if (f->IsDomainBoundary)
			this->BoundaryFaces.push_back(f);
		else
			this->InteriorFaces.push_back(f);
	}

	double SkeletonMeasure()
	{
		double measure = 0;
		for (Face<Dim>* face : this->Faces)
			measure += face->Measure();
		return measure;
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
		fprintf(file, "Number x1    y1    x2    y2 IsDomainBoundary\n");
		for (Face<Dim>* f : this->Faces)
			f->ExportFaceToMatlab(file);
		fclose(file);
		cout << "Faces exported to \t" << filePath << endl;
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
		{
			e->UnitTests();

			for (auto f : e->Faces)
			{
				auto n = e->OuterNormalVector(f);

				auto neighbour = f->GetNeighbour(e);
				if (neighbour != nullptr)
				{
					auto n2 = neighbour->OuterNormalVector(f);
					assert(abs(n.dot(n2) + 1) < 1e-15);
				}
				else
					assert(e->IsOnBoundary());
			}
		}

		for (auto f : this->BoundaryFaces)
			assert(f->IsDomainBoundary);

		for (auto f : this->InteriorFaces)
			assert(!f->IsDomainBoundary);

		assert(this->BoundaryFaces.size() + this->InteriorFaces.size() == this->Faces.size());

		if (CoarseMesh)
		{
			CoarseMesh->SanityCheck();

			for (Element<Dim>* fe : this->Elements)
			{
				assert(fe->CoarserElement != nullptr);

				bool feIsReferenced = false;
				for (Element<Dim>* e : fe->CoarserElement->FinerElements)
				{
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
			
			for (Face<Dim>* cf : CoarseMesh->Faces)
			{
				assert(cf->FinerFaces.size() > 0);
				for (Face<Dim>* ff : cf->FinerFaces)
					assert(!ff->IsRemovedOnCoarserGrid);
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