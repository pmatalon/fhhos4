#pragma once
#include "../GMSHMesh.h"
using namespace std;

class Square_GMSHUnstructuredTriangularMesh : public GMSHMesh<2>
{
public:
	Square_GMSHUnstructuredTriangularMesh() : GMSHMesh("square_unstruct_tri.msh")
	{
		this->_description = "GMSH unstructured triangular";
		this->_fileNamePart = "gmsh-uns-tri";
	}

	void RefineMeshBySplitting() override
	{
		GMSHMesh::RefineMeshBySplitting();
		this->FineMesh->ComesFrom.nFineElementsByCoarseElement = 4;
		this->FineMesh->ComesFrom.nFineFacesAddedByCoarseElement = 3;
		this->FineMesh->ComesFrom.nFineFacesByKeptCoarseFace = 2;
	}
};