#pragma once
#include <vector>
#include "Element.h"
#include "Face.h"
#include "../Utils/ParallelLoop.h"
using namespace std;

enum class CoarseningStrategy : unsigned
{
	Standard = 0,
	Agglomeration = 1
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

	DiffusionPartition<Dim>* _diffusionPartition = nullptr;

	Mesh<Dim>* CoarseMesh = NULL;

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
	}

	void ExportFacesToMatlab(string outputDirectory, bool dummy)
	{
		string filePath = outputDirectory + "/faces" + to_string(Dim) + "D_" + FileNamePart() + ".dat";
		ExportFacesToMatlab(filePath);
	}
	void ExportFacesToMatlab(string filePath)
	{
		FILE* file = fopen(filePath.c_str(), "w");
		fprintf(file, "Number OriginX OriginY OriginZ WidthX WidthY WidthZ Orientation Boundary\n");
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