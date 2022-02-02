#pragma once
#include "../TestCases/Diffusion/DiffusionField.h"
using namespace std;

class Square4quadrantsGeometry
{
public:
	static vector<PhysicalGroup<2>*> PhysicalParts()
	{
		vector<PhysicalGroup<2>*> physicalParts;
		physicalParts.push_back(new PhysicalGroup<2>(1, "quadrantBottomLeft"));
		physicalParts.push_back(new PhysicalGroup<2>(2, "quadrantBottomRight"));
		physicalParts.push_back(new PhysicalGroup<2>(3, "quadrantTopRight"));
		physicalParts.push_back(new PhysicalGroup<2>(4, "quadrantTopLeft"));
		return physicalParts;
	}

	static DiffusionField<2> DiffField(double heterogeneityRatio, double anisotropyRatio, double anisotropyAngle)
	{
		double kappa1 = heterogeneityRatio;
		double kappa2 = 1;
		Tensor<2> diffTensor1(kappa1, anisotropyRatio, anisotropyAngle);
		Tensor<2> diffTensor2(kappa2, anisotropyRatio, anisotropyAngle);

		map<string, Tensor<2>> tensors;
		tensors.insert({ "quadrantTopLeft", diffTensor1 });
		tensors.insert({ "quadrantTopRight", diffTensor2 });
		tensors.insert({ "quadrantBottomLeft", diffTensor2 });
		tensors.insert({ "quadrantBottomRight", diffTensor1 });

		return DiffusionField<2>(tensors);
	}
};