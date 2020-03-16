#pragma once
#include "GMSHTetrahedralMesh.h"
using namespace std;

class Cube_GMSHTetrahedralMesh : public GMSHTetrahedralMesh
{
public:
	Cube_GMSHTetrahedralMesh() :
		GMSHTetrahedralMesh("cube4tetra.geo", "GMSH tetrahedral", "gmsh-tetra")
	{}

private:
	Cube_GMSHTetrahedralMesh(string description, string fileNamePart) : GMSHTetrahedralMesh(description, fileNamePart)
	{}

protected:
	virtual GMSHMesh<3>* CreateEmptyGMSHMesh() override
	{
		return new Cube_GMSHTetrahedralMesh(this->_description, this->_fileNamePart);
	}
};