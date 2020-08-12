#pragma once
#include "../GMSHMesh.h"
using namespace std;

class Square_GMSHTriangularMesh : public GMSHMesh<2>
{
public:
	Square_GMSHTriangularMesh() : GMSHMesh("2D/square_stri.geo")
	{
		this->_description = "GMSH structured triangular";
		this->_fileNamePart = "square-gmsh-stri";
		this->_geometryDescription = "Square";
	}

	void RefineMeshBySplitting() override
	{
		GMSHMesh::RefineMeshBySplitting();
		this->FineMesh->ComesFrom.nFineElementsByCoarseElement = 4;
		this->FineMesh->ComesFrom.nFineFacesAddedByCoarseElement = 3;
		this->FineMesh->ComesFrom.nFineFacesByKeptCoarseFace = 2;
	}
};