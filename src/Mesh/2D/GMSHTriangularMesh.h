#pragma once
#include "../GMSHMesh.h"
using namespace std;

class GMSHTriangularMesh : public GMSHMesh<2>
{
public:
	GMSHTriangularMesh() : GMSHMesh("square4tri.geo")
	{
		this->_description = "GMSH triangular";
		this->_fileNamePart = "gmsh-tri";
	}

	void RefineMesh() override
	{
		GMSHMesh::RefineMesh();
		this->FineMesh->ComesFrom.nFineElementsByCoarseElement = 4;
		this->FineMesh->ComesFrom.nFineFacesAddedByCoarseElement = 3;
		this->FineMesh->ComesFrom.nFineFacesByKeptCoarseFace = 2;
	}
};