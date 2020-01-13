#pragma once
#include "../GMSHMesh.h"
using namespace std;

class GMSHUnstructuredTriangularMesh : public GMSHMesh<2>
{
public:
	GMSHUnstructuredTriangularMesh() : GMSHMesh("square_unstruct_tri.msh")
	{
		this->_description = "GMSH unstructured triangular";
		this->_fileNamePart = "gmsh-uns-tri";
	}

};