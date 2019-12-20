#pragma once
#include "Triangle.h"
#include "../PolyhedralMesh.h"
#include <gmsh.h>
using namespace std;

enum GMSHElementTypes
{
	Segment = 1,
	Triangle = 2,
	Quadrilateral = 3,
	Tetrahedron = 4
};

class GMSHTriangularMesh : public PolyhedralMesh<2>
{
public:
	GMSHTriangularMesh() : PolyhedralMesh()
	{
		gmsh::initialize();
		gmsh::open("/mnt/c/Users/pierr/Documents/Source/Repos/dghho/data/meshes/square.msh");

		vector<int> elementTypes;
		vector<vector<size_t>> elementTags;
		vector<vector<size_t>> nodeTags;
		gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags, 2);

		for (int i = 0; i< elementTypes.size(); i++)
		{
			int elemType = elementTypes[i];
			vector<size_t> elements = elementTags[i];
			
		}
	}

	~GMSHTriangularMesh()
	{
		gmsh::finalize();
	}

	string Description() override
	{
		return "GMSH triangular";
	}

	string FileNamePart() override
	{
		return "GMSH_tri";
	}

	double H() override
	{
		//return sqrt(2) / (double)this->Nx;
	}

	void CoarsenMesh(CoarseningStrategy strategy) override
	{
		assert(false);
	}

};