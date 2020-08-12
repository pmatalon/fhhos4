#pragma once
#include "../GMSHMesh.h"
using namespace std;

class Square4quadrants_GMSHUnstructTriangularMesh : public GMSHMesh<2>
{
public:
	Square4quadrants_GMSHUnstructTriangularMesh(BigNumber n) : GMSHMesh("2D/square.geo", 1.0 / n)
	{
		this->_description = "Square 4 quandrants - GMSH unstructured triangular";
		this->_fileNamePart = "square4quadrants-gmsh-tri";
		this->_geometryDescription = "Square 4 quadrants";
	}

	Square4quadrants_GMSHUnstructTriangularMesh() : GMSHMesh("2D/square4quadrants_tri.msh")
	{
		this->_description = "Square 4 quandrants - GMSH unstructured triangular";
		this->_fileNamePart = "square4quadrants-gmsh-tri";
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