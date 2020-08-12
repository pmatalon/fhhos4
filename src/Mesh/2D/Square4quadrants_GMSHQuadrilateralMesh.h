#pragma once
#include "../GMSHMesh.h"
using namespace std;

class Square4quadrants_GMSHQuadrilateralMesh : public GMSHMesh<2>
{
public:
	Square4quadrants_GMSHQuadrilateralMesh() : GMSHMesh("2D/square4quadrants_quad.msh")
	{
		this->_description = "GMSH quadrilateral";
		this->_fileNamePart = "square4quadrants-gmsh-quad";
		this->_geometryDescription = "Square 4 quadrants";
	}

	void RefineMeshBySplitting() override
	{
		GMSHMesh::RefineMeshBySplitting();
		this->FineMesh->ComesFrom.nFineElementsByCoarseElement = 4;
		this->FineMesh->ComesFrom.nFineFacesAddedByCoarseElement = 4;
		this->FineMesh->ComesFrom.nFineFacesByKeptCoarseFace = 2;
	}

};