#pragma once
#include "../GMSHMesh.h"
using namespace std;

class Square_GMSHCartesianMesh : public GMSHMesh<2>
{
public:
	Square_GMSHCartesianMesh() : GMSHMesh("square4cart.geo")
	{
		this->_description = "GMSH Cartesian";
		this->_fileNamePart = "gmsh-cart";
	}

	void RefineMeshBySplitting() override
	{
		GMSHMesh::RefineMeshBySplitting();
		this->FineMesh->ComesFrom.nFineElementsByCoarseElement = 4;
		this->FineMesh->ComesFrom.nFineFacesAddedByCoarseElement = 4;
		this->FineMesh->ComesFrom.nFineFacesByKeptCoarseFace = 2;
	}

};