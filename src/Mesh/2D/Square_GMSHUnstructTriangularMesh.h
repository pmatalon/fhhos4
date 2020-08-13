#pragma once
#include "../GMSHMesh.h"
using namespace std;

class Square_GMSHUnstructTriangularMesh : public GMSHMesh<2>
{
public:
	Square_GMSHUnstructTriangularMesh(BigNumber n) :
		GMSHMesh("2D/square_tri.geo", "GMSH unstructured triangular", "square-gmsh-tri", "Square", n)
	{}

	void RefineMeshBySplitting() override
	{
		GMSHMesh::RefineMeshBySplitting();
		this->FineMesh->ComesFrom.nFineElementsByCoarseElement = 4;
		this->FineMesh->ComesFrom.nFineFacesAddedByCoarseElement = 3;
		this->FineMesh->ComesFrom.nFineFacesByKeptCoarseFace = 2;
	}
};