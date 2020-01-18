#pragma once
#include "../GMSHMesh.h"
using namespace std;

class GMSHTetrahedralMesh : public GMSHMesh<3>
{
public:
	GMSHTetrahedralMesh() : GMSHMesh("cube6.geo")
	{
		this->_description = "GMSH tetrahedral";
		this->_fileNamePart = "gmsh-tetra";
	}

	void RefineMesh() override
	{
		GMSHMesh::RefineMesh();
		this->FineMesh->ComesFrom.nFineElementsByCoarseElement = 8;
		this->FineMesh->ComesFrom.nFineFacesAddedByCoarseElement = 8;
		this->FineMesh->ComesFrom.nFineFacesByKeptCoarseFace = 4;
	}
};