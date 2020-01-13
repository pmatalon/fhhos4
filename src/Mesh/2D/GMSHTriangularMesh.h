#pragma once
#include "../GMSHMesh.h"
using namespace std;

class GMSHTriangularMesh : public GMSHMesh<2>
{
public:
	GMSHTriangularMesh() : GMSHMesh("square_tri_n2.msh")
	{
		this->_description = "GMSH triangular";
		this->_fileNamePart = "gmsh-tri";
	}

};