#pragma once
#include "SquareTestCase.h"
using namespace std;

class RoomWithWallTestCase : public SquareTestCase
{
public:
	RoomWithWallTestCase(ProblemArguments pb) :
		SquareTestCase(pb)
	{
		// Diffusion field
		double kappa1 = pb.HeterogeneityRatio;
		double kappa2 = 1;
		Tensor<2>* diffTensor1 = new Tensor<2>(kappa1, pb.AnisotropyRatio, pb.AnisotropyAngle);
		Tensor<2>* diffTensor2 = new Tensor<2>(kappa2, pb.AnisotropyRatio, pb.AnisotropyAngle);

		map<string, Tensor<2>*> tensors;
		tensors.insert({ "wall", diffTensor1 });
		tensors.insert({ "air", diffTensor2 });

		this->DiffField = DiffusionField<2>(tensors);
	}
};
