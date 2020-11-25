#pragma once
#include "GMSHMesh.h"
using namespace std;

class Square4quadrants_GMSHTriangularMesh : public GMSHMesh<2>
{
public:
	Square4quadrants_GMSHTriangularMesh(BigNumber n) : 
		GMSHMesh("2D/square4quadrants_stri.geo", "GMSH structured triangular", "square4quadrants_gmsh_stri", "Square 4 quadrants", n / 2) // n/2 because it builds n subdivisions in each quadrant
	{}

	void RefineMeshBySplitting() override
	{
		GMSHMesh::RefineMeshBySplitting();
		this->FineMesh->ComesFrom.nFineElementsByCoarseElement = 4;
		this->FineMesh->ComesFrom.nFineFacesAddedByCoarseElement = 3;
		this->FineMesh->ComesFrom.nFineFacesByKeptCoarseFace = 2;
	}
};