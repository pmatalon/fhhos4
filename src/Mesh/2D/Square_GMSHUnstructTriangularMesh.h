#pragma once
#include "../GMSHMesh.h"
using namespace std;

class Square_GMSHUnstructTriangularMesh : public GMSHMesh<2>
{
public:
	Square_GMSHUnstructTriangularMesh(BigNumber n) : GMSHMesh("2D/square.geo", "GMSH unstructured triangular", "gmsh-uns-tri", 1.0/n)
	{}

	void RefineMeshBySplitting() override
	{
		GMSHMesh::RefineMeshBySplitting();
		this->FineMesh->ComesFrom.nFineElementsByCoarseElement = 4;
		this->FineMesh->ComesFrom.nFineFacesAddedByCoarseElement = 3;
		this->FineMesh->ComesFrom.nFineFacesByKeptCoarseFace = 2;
	}
};