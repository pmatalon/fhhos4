#pragma once
#include "GMSHMesh.h"
using namespace std;

class Square_GMSHCartesianMesh : public GMSHMesh<2>
{
public:
	Square_GMSHCartesianMesh(BigNumber n) :
		GMSHMesh(nullptr, "2D/square_cart.geo", "GMSH Cartesian", "square_gmsh_cart", "Square", n)
	{}

	void RefineMeshBySplitting() override
	{
		GMSHMesh::RefineMeshBySplitting();
		this->FineMesh->ComesFrom.nFineElementsByCoarseElement = 4;
		this->FineMesh->ComesFrom.nFineFacesAddedByCoarseElement = 4;
		this->FineMesh->ComesFrom.nFineFacesByKeptCoarseFace = 2;
	}

};