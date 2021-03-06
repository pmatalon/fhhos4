#pragma once
#include "GMSHMesh.h"
using namespace std;

class Square_GMSHQuadrilateralMesh : public GMSHMesh<2>
{
public:
	Square_GMSHQuadrilateralMesh(BigNumber n) : 
		GMSHMesh(nullptr, "2D/square_quad.geo", "GMSH quadrilateral", "square_gmsh_quad", "Square", n)
	{}

	void RefineMeshBySplitting() override
	{
		GMSHMesh::RefineMeshBySplitting();
		this->FineMesh->ComesFrom.nFineElementsByCoarseElement = 4;
		this->FineMesh->ComesFrom.nFineFacesAddedByCoarseElement = 4;
		this->FineMesh->ComesFrom.nFineFacesByKeptCoarseFace = 2;
	}

};