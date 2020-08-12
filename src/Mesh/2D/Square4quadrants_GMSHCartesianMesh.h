#pragma once
#include "../GMSHMesh.h"
using namespace std;

class Square4quadrants_GMSHCartesianMesh : public GMSHMesh<2>
{
public:
	Square4quadrants_GMSHCartesianMesh() : GMSHMesh("2D/square4quadrants_cart.geo")
	{
		this->_description = "GMSH Cartesian";
		this->_fileNamePart = "square4quadrants-gmsh-cart";
		this->_geometryDescription = "Square 4 quadrants";
	}

	void RefineMeshBySplitting() override
	{
		GMSHMesh::RefineMeshBySplitting();
		this->FineMesh->ComesFrom.nFineElementsByCoarseElement = 4;
		this->FineMesh->ComesFrom.nFineFacesAddedByCoarseElement = 4;
		this->FineMesh->ComesFrom.nFineFacesByKeptCoarseFace = 2;
	}

};