#pragma once
#include "../GMSHMesh.h"
using namespace std;

class Square_GMSHCartesianMesh : public GMSHMesh<2>
{
public:
	Square_GMSHCartesianMesh() : GMSHMesh("2D/square_cart.geo")
	{
		this->_description = "GMSH Cartesian";
		this->_fileNamePart = "square-gmsh-cart";
		this->_geometryDescription = "Square";
	}

	void RefineMeshBySplitting() override
	{
		GMSHMesh::RefineMeshBySplitting();
		this->FineMesh->ComesFrom.nFineElementsByCoarseElement = 4;
		this->FineMesh->ComesFrom.nFineFacesAddedByCoarseElement = 4;
		this->FineMesh->ComesFrom.nFineFacesByKeptCoarseFace = 2;
	}

};