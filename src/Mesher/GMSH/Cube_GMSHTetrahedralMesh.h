#pragma once
#include "GMSHTetrahedralMesh.h"
using namespace std;

class Cube_GMSHTetrahedralMesh : public GMSHTetrahedralMesh
{
public:
	Cube_GMSHTetrahedralMesh(BigNumber n) :
		GMSHTetrahedralMesh("3D/cube_tetra.geo", "GMSH unstructured tetrahedral", "gmsh_tetra", "Cube", n)
	{}

private:
	Cube_GMSHTetrahedralMesh(string description, string fileNamePart, string geometryDescription) : 
		GMSHTetrahedralMesh(description, fileNamePart, geometryDescription)
	{}

protected:
	virtual GMSHMesh<3>* CreateNewGMSHMesh() override
	{
		return new Cube_GMSHTetrahedralMesh(this->_description, this->_fileNamePart, this->_geometryDescription);
	}
};