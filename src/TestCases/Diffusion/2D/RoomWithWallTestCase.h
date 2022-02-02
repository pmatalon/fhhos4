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
		Tensor<2> diffTensor1(pb.HeterogeneityRatio, pb.AnisotropyRatio, pb.AnisotropyAngle);
		Tensor<2> diffTensor2(                    1, pb.AnisotropyRatio, pb.AnisotropyAngle);

		map<string, Tensor<2>> tensors;
		tensors.insert({ "wall", diffTensor1 });
		tensors.insert({ "air", diffTensor2 });

		this->DiffField = DiffusionField<2>(tensors);
	}
};
