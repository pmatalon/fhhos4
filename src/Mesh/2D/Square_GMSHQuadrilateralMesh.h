#pragma once
#include "../GMSHMesh.h"
using namespace std;

class Square_GMSHQuadrilateralMesh : public GMSHMesh<2>
{
public:
	Square_GMSHQuadrilateralMesh() : GMSHMesh("square_quad.msh")
	{
		this->_description = "GMSH quadrilateral";
		this->_fileNamePart = "gmsh-quad";
	}

	void RefineMeshBySplitting() override
	{
		GMSHMesh::RefineMeshBySplitting();
		this->FineMesh->ComesFrom.nFineElementsByCoarseElement = 4;
		this->FineMesh->ComesFrom.nFineFacesAddedByCoarseElement = 4;
		this->FineMesh->ComesFrom.nFineFacesByKeptCoarseFace = 2;
	}

};