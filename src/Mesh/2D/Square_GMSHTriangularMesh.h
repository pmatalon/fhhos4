#pragma once
#include "../GMSHMesh.h"
using namespace std;

class Square_GMSHTriangularMesh : public GMSHMesh<2>
{
public:
	Square_GMSHTriangularMesh() : GMSHMesh("2D/square4tri.geo")
	{
		this->_description = "GMSH triangular";
		this->_fileNamePart = "gmsh-tri";
	}

	void RefineMeshBySplitting() override
	{
		GMSHMesh::RefineMeshBySplitting();
		this->FineMesh->ComesFrom.nFineElementsByCoarseElement = 4;
		this->FineMesh->ComesFrom.nFineFacesAddedByCoarseElement = 3;
		this->FineMesh->ComesFrom.nFineFacesByKeptCoarseFace = 2;
	}
};