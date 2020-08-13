#pragma once
#include "../GMSHMesh.h"
using namespace std;

class Square4quadrants_GMSHUnstructTriangularMesh : public GMSHMesh<2>
{
public:
	Square4quadrants_GMSHUnstructTriangularMesh(BigNumber n) : 
		GMSHMesh("2D/square4quadrants_tri.geo", "GMSH unstructured triangular", "square4quadrants-gmsh-tri", "Square 4 quadrants", n / 2) // n/2 because it builds n subdivisions in each quadrant
	{}

	void RefineMeshBySplitting() override
	{
		GMSHMesh::RefineMeshBySplitting();
		this->FineMesh->ComesFrom.nFineElementsByCoarseElement = 4;
		this->FineMesh->ComesFrom.nFineFacesAddedByCoarseElement = 3;
		this->FineMesh->ComesFrom.nFineFacesByKeptCoarseFace = 2;
	}
};