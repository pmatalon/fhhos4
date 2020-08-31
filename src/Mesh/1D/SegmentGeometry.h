#pragma once
#include "../../Problem/DiffusionField.h"
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
		Tensor<1>* tensorLeft = new Tensor<1>(heterogeneityRatio);
		Tensor<1>* tensorRight = new Tensor<1>(1);
		return DiffusionField<1>("left", tensorLeft, "right", tensorRight);
	}
};