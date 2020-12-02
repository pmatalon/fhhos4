#pragma once
#include "GMSHMesh.h"
using namespace std;

class Square4quadrants_GMSHQuadrilateralMesh : public GMSHMesh<2>
{
public:
	Square4quadrants_GMSHQuadrilateralMesh(BigNumber n) : 
		GMSHMesh(nullptr, "2D/square4quadrants_quad.geo", "GMSH quadrilateral", "square4quadrants_gmsh_quad", "Square 4 quadrants", n / 2) // n/2 because it builds n subdivisions in each quadrant
	{}

	void RefineMeshBySplitting() override
	{
		GMSHMesh::RefineMeshBySplitting();
		this->FineMesh->ComesFrom.nFineElementsByCoarseElement = 4;
		this->FineMesh->ComesFrom.nFineFacesAddedByCoarseElement = 4;
		this->FineMesh->ComesFrom.nFineFacesByKeptCoarseFace = 2;
	}
};