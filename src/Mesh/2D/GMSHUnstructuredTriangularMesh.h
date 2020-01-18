#pragma once
#include "../GMSHMesh.h"
using namespace std;

class GMSHUnstructuredTriangularMesh : public GMSHMesh<2>
{
public:
	GMSHUnstructuredTriangularMesh() : GMSHMesh("square_unstruct_tri.msh")
	{
		this->_description = "GMSH unstructured triangular";
		this->_fileNamePart = "gmsh-uns-tri";
	}

	void RefineMesh() override
	{
		GMSHMesh::RefineMesh();
		this->FineMesh->ComesFrom.nFineElementsByCoarseElement = 4;
		this->FineMesh->ComesFrom.nFineFacesAddedByCoarseElement = 3;
		this->FineMesh->ComesFrom.nFineFacesByKeptCoarseFace = 2;
	}
};