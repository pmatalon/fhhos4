#pragma once
#include <vector>
#include "Element.h"
#include "Face.h"
using namespace std;

enum class CoarseningStrategy : unsigned
{
	AgglomerationAndMergeColinearFaces = 0,
	AgglomerationAndKeepFineFaces = 1
};

template <int Dim>
class Mesh
{
public:
	vector<Element<Dim>*> Elements;
	vector<Face<Dim>*> Faces;
	vector<Face<Dim>*> BoundaryFaces;
	vector<Face<Dim>*> InteriorFaces;

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

	void ExportFacesToMatlab(string outputDirectory)
	{
		string filePath = outputDirectory + "/faces" + to_string(Dim) + "D_" + FileNamePart() + ".dat";
		FILE* file = fopen(filePath.c_str(), "w");
		fprintf(file, "Number OriginX OriginY OriginZ Orientation Boundary\n");
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
	}
};