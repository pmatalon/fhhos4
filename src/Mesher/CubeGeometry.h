#pragma once
#include "../TestCases/Diffusion/DiffusionField.h"
using namespace std;

class CubeGeometry
{
public:
	static vector<PhysicalGroup<3>*> PhysicalParts()
	{
		vector<PhysicalGroup<3>*> physicalParts;
		physicalParts.push_back(new PhysicalGroup<3>(1, "domain"));
		return physicalParts;
	}

	static vector<BoundaryGroup*> BoundaryParts()
	{
		vector<BoundaryGroup*> boundaryParts;
		return boundaryParts;
	}

	static DiffusionField<3> DiffField(double anisotropyRatio, double anisotropyAngle)
	{
		return DiffusionField<3>(anisotropyRatio, anisotropyAngle);
	}
};