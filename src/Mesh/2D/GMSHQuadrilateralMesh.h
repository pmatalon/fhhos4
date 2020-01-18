#pragma once
#include "../GMSHMesh.h"
using namespace std;

class GMSHQuadrilateralMesh : public GMSHMesh<2>
{
public:
	GMSHQuadrilateralMesh() : GMSHMesh("square_quad.msh")
	{
		this->_description = "GMSH quadrilateral";
		this->_fileNamePart = "gmsh-quad";
	}

	void RefineMesh() override
	{
		GMSHMesh::RefineMesh();
		this->FineMesh->ComesFrom.nFineElementsByCoarseElement = 4;
		this->FineMesh->ComesFrom.nFineFacesAddedByCoarseElement = 4;
		this->FineMesh->ComesFrom.nFineFacesByKeptCoarseFace = 2;
	}

};