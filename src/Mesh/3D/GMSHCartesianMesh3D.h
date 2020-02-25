#pragma once
#include "../GMSHMesh.h"
using namespace std;

class GMSHCartesianMesh3D : public GMSHMesh<3>
{
public:
	GMSHCartesianMesh3D() : GMSHMesh("cube4cart.geo", "GMSH Cartesian", "gmsh-cart")
	{}
private:
	GMSHCartesianMesh3D(string description, string fileNamePart) : GMSHMesh(description, fileNamePart)
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
	virtual GMSHMesh<3>* CreateEmptyGMSHMesh() override
	{
		return new GMSHCartesianMesh3D(this->_description, this->_fileNamePart);
	}
};