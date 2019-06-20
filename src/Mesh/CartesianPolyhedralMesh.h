#pragma once
#include "Mesh.h"
using namespace std;

template <int Dim>
class CartesianPolyhedralMesh : public Mesh<Dim>
{
public:
	CartesianPolyhedralMesh() : Mesh<Dim>() {}

	virtual string Description() { return "CartesianPolyhedralMesh"; }
	virtual string FileNamePart() { return "nX"; };
	virtual double H() { assert(false); };

	virtual void CoarsenMesh(CoarseningStrategy strategy)
	{

	}

};