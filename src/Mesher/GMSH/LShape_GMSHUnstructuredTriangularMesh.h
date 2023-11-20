#pragma once
#include "GMSHMesh.h"
using namespace std;

class LShape_GMSHUnstructuredTriangularMesh : public GMSHMesh<2>
{
public:
	LShape_GMSHUnstructuredTriangularMesh(BigNumber n) : 
		GMSHMesh(nullptr, "2D/L_shape_tri.geo", "GMSH unstructured triangular", "L_shape_gmsh_tri", "L-shape", n/4)  // divided by 2 because it builds n subdivisions in each quadrant + divided by 2 because the domain is (-1,1) instead of (0,1)
	{}

	void RefineMeshBySplitting() override
	{
		GMSHMesh::RefineMeshBySplitting();
		this->FineMesh->ComesFrom.nFineElementsByCoarseElement = 4;
		this->FineMesh->ComesFrom.nFineFacesAddedByCoarseElement = 3;
		this->FineMesh->ComesFrom.nFineFacesByKeptCoarseFace = 2;
	}
};