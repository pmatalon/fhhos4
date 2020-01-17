#pragma once
#include "../GMSHMesh.h"
using namespace std;

class GMSHTetrahedralMesh : public GMSHMesh<3>
{
public:
	GMSHTetrahedralMesh() : GMSHMesh("cube6.geo")
	{
		this->_description = "GMSH tetrahedral";
		this->_fileNamePart = "gmsh-tetra";
	}

};