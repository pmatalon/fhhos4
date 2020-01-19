#pragma once
#include "../GMSHMesh.h"
using namespace std;

class GMSHCartesianMesh : public GMSHMesh<2>
{
public:
	GMSHCartesianMesh() : GMSHMesh("square_cart_n2.msh")
	{
		this->_description = "GMSH Cartesian";
		this->_fileNamePart = "gmsh-cart";
	}

	void RefineMesh() override
	{
		GMSHMesh::RefineMesh();
		this->FineMesh->ComesFrom.nFineElementsByCoarseElement = 4;
		this->FineMesh->ComesFrom.nFineFacesAddedByCoarseElement = 4;
		this->FineMesh->ComesFrom.nFineFacesByKeptCoarseFace = 2;
	}

};