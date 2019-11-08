#pragma once
#include "Polygon.h"
#include "../Mesh.h"
using namespace std;

class PolygonalMesh : public Mesh<2>
{
public:
	const int AgglomerateSize = 4;

	PolygonalMesh() : Mesh()
	{}

	string Description()
	{
		return "Polygonal";
	}

	string FileNamePart()
	{
		return "polygonal";
	}

	double H()
	{
		assert(false);
	}

	void CoarsenMesh(CoarseningStrategy strategy)
	{
		if (strategy == CoarseningStrategy::Standard)
			CoarsenByAgglomerationAndMergeColinearFaces();
		else if (strategy == CoarseningStrategy::Agglomeration)
			CoarsenByAgglomerationAndKeepFineFaces();
		else
			assert(false && "Coarsening strategy not implemented!");
		this->CoarseMesh->SetDiffusionCoefficient(this->_diffusionPartition);
		this->CoarseMesh->SetBoundaryConditions(this->_boundaryConditions);
	}

	void CoarsenByAgglomerationAndMergeColinearFaces()
	{
		if (this->Elements.size() <= AgglomerateSize)
		{
			cout << "Error: impossible to build coarse mesh. Only " << this->Elements.size() << " element(s) left." << endl;
		}
		else
		{
			PolygonalMesh* coarseMesh = new PolygonalMesh();


			this->CoarseMesh = coarseMesh;
		}
	}



	void CoarsenByAgglomerationAndKeepFineFaces()
	{
		if (this->Elements.size() <= AgglomerateSize)
		{
			cout << "Error: impossible to build coarse mesh. Only " << this->Elements.size() << " element(s) left." << endl;
		}
		else
		{
			PolygonalMesh* coarseMesh = new PolygonalMesh();

			this->CoarseMesh = coarseMesh;
		}
	}
};