#pragma once
#include "../TestCases/Diffusion/DiffusionField.h"
using namespace std;

class SegmentGeometry
{
public:
	static vector<PhysicalGroup<1>*> PhysicalParts()
	{
		vector<PhysicalGroup<1>*> physicalParts;
		physicalParts.push_back(new PhysicalGroup<1>(1, "left"));
		physicalParts.push_back(new PhysicalGroup<1>(2, "right"));
		return physicalParts;
	}

	static DiffusionField<1> DiffField(double heterogeneityRatio)
	{
		Tensor<1> tensorLeft(heterogeneityRatio);
		Tensor<1> tensorRight(1);
		return DiffusionField<1>("left", tensorLeft, "right", tensorRight);
	}
};