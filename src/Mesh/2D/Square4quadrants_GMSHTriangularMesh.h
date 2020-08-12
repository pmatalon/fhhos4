#pragma once
#include "../GMSHMesh.h"
using namespace std;

class Square4quadrants_GMSHTriangularMesh : public GMSHMesh<2>
{
public:
	Square4quadrants_GMSHTriangularMesh(BigNumber n) : GMSHMesh("2D/square4quadrants_stri.geo", 1.0/n)
	{
		this->_description = "GMSH structured triangular";
		this->_fileNamePart = "square4quadrants-gmsh-stri";
		this->_geometryDescription = "Square 4 quadrants";
	}

	void RefineMeshBySplitting() override
	{
		GMSHMesh::RefineMeshBySplitting();
		this->FineMesh->ComesFrom.nFineElementsByCoarseElement = 4;
		this->FineMesh->ComesFrom.nFineFacesAddedByCoarseElement = 3;
		this->FineMesh->ComesFrom.nFineFacesByKeptCoarseFace = 2;
	}
};