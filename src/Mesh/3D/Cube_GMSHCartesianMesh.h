#pragma once
#include "../GMSHMesh.h"
using namespace std;

class Cube_GMSHCartesianMesh : public GMSHMesh<3>
{
public:
	Cube_GMSHCartesianMesh() : GMSHMesh("3D/cube4cart.geo", "GMSH Cartesian", "gmsh-cart")
	{}
private:
	Cube_GMSHCartesianMesh(string description, string fileNamePart) : GMSHMesh(description, fileNamePart)
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
		return new Cube_GMSHCartesianMesh(this->_description, this->_fileNamePart);
	}
};