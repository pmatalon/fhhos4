#pragma once
#include "GMSHMesh.h"
using namespace std;

class Square_GMSHTriangularMesh : public GMSHMesh<2>
{
public:
	Square_GMSHTriangularMesh(BigNumber n) : 
		GMSHMesh("2D/square_stri.geo", "GMSH structured triangular", "square_gmsh_stri", "Square", n)
	{}

	void RefineMeshBySplitting() override
	{
		GMSHMesh::RefineMeshBySplitting();
		this->FineMesh->ComesFrom.nFineElementsByCoarseElement = 4;
		this->FineMesh->ComesFrom.nFineFacesAddedByCoarseElement = 3;
		this->FineMesh->ComesFrom.nFineFacesByKeptCoarseFace = 2;
	}
};