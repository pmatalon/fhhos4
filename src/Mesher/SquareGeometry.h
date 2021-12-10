#pragma once
#include "../TestCases/Diffusion/DiffusionField.h"
using namespace std;

class SquareGeometry
{
public:
	static vector<PhysicalGroup<2>*> PhysicalParts()
	{
		vector<PhysicalGroup<2>*> physicalParts;
		physicalParts.push_back(new PhysicalGroup<2>(1, "domain"));
		return physicalParts;
	}

	static vector<BoundaryGroup*> BoundaryParts()
	{
		vector<BoundaryGroup*> boundaryParts;
		boundaryParts.push_back(new BoundaryGroup(1, "bottomBoundary"));
		boundaryParts.push_back(new BoundaryGroup(2, "rightBoundary"));
		boundaryParts.push_back(new BoundaryGroup(3, "topBoundary"));
		boundaryParts.push_back(new BoundaryGroup(4, "leftBoundary"));
		return boundaryParts;
	}

	static DiffusionField<2> DiffField(double anisotropyRatio, double anisotropyAngle)
	{
		return DiffusionField<2>(anisotropyRatio, anisotropyAngle);
	}
};