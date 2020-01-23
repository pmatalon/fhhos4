#pragma once
#include "../GMSHMesh.h"
using namespace std;

class GMSHTetrahedralMesh : public GMSHMesh<3>
{
public:
	GMSHTetrahedralMesh() : GMSHMesh("cube4tetra.geo", "GMSH tetrahedral", "gmsh-tetra")
	{}
private:
	GMSHTetrahedralMesh(string description, string fileNamePart) : GMSHMesh(description, fileNamePart)
	{}

public:
	void RefineMesh() override
	{
		GMSHMesh::RefineMesh();
		this->FineMesh->ComesFrom.nFineElementsByCoarseElement = 8;
		this->FineMesh->ComesFrom.nFineFacesAddedByCoarseElement = 8;
		this->FineMesh->ComesFrom.nFineFacesByKeptCoarseFace = 4;
	}

protected:
	virtual GMSHMesh<3>* CreateEmptyMesh() override
	{
		return new GMSHTetrahedralMesh(this->_description, this->_fileNamePart);
	}
};