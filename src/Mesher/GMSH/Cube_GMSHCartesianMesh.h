#pragma once
#include "GMSHMesh.h"
using namespace std;

class Cube_GMSHCartesianMesh : public GMSHMesh<3>
{
public:
	Cube_GMSHCartesianMesh(BigNumber n) : 
		GMSHMesh(nullptr, "3D/cube_cart.geo", "GMSH Cartesian", "gmsh_cart", "Cube", n)
	{}
private:
	Cube_GMSHCartesianMesh(string description, string fileNamePart, string geometryDescription) :
		GMSHMesh(nullptr, description, fileNamePart, geometryDescription)
	{}

public:
	void RefineMeshBySplitting() override
	{
		GMSHMesh::RefineMeshBySplitting();
		this->FineMesh->ComesFrom.nFineElementsByCoarseElement = 8;
		this->FineMesh->ComesFrom.nFineFacesAddedByCoarseElement = 12;
		this->FineMesh->ComesFrom.nFineFacesByKeptCoarseFace = 4;
	}

protected:
	virtual GMSHMesh<3>* CreateNewGMSHMesh() override
	{
		return new Cube_GMSHCartesianMesh(this->_description, this->_fileNamePart, this->_geometryDescription);
	}
};