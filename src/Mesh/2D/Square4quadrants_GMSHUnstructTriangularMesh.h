#pragma once
#include "../GMSHMesh.h"
using namespace std;

class Square4quadrants_GMSHUnstructTriangularMesh : public GMSHMesh<2>
{
public:
	Square4quadrants_GMSHUnstructTriangularMesh() : GMSHMesh("2D/square4quadrants_uns_tri.msh")
	{
		this->_description = "Square 4 quandrants - GMSH unstructured triangular";
		this->_fileNamePart = "sq4quadrants-gmsh-uns-tri";
	}

	void RefineMeshBySplitting() override
	{
		GMSHMesh::RefineMeshBySplitting();
		this->FineMesh->ComesFrom.nFineElementsByCoarseElement = 4;
		this->FineMesh->ComesFrom.nFineFacesAddedByCoarseElement = 3;
		this->FineMesh->ComesFrom.nFineFacesByKeptCoarseFace = 2;
	}
};