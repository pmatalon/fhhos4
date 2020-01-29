#pragma once
#include "Smoother.h"
using namespace std;

class SmootherFactory
{

public:
	static Smoother* Create(string smootherCode, int nSmootherIterations, int blockSize)
	{
		if (smootherCode.compare(GaussSeidel::Code()) == 0)
			return new GaussSeidelSmoother(nSmootherIterations);
		if (smootherCode.compare(ReverseGaussSeidel::Code()) == 0)
			return new ReverseGaussSeidelSmoother(nSmootherIterations);
		if (smootherCode.compare(BlockGaussSeidel::Code()) == 0)
			return new BlockGaussSeidelSmoother(blockSize, nSmootherIterations);
		if (smootherCode.compare(ReverseBlockGaussSeidel::Code()) == 0)
			return new ReverseBlockGaussSeidelSmoother(blockSize, nSmootherIterations);
		if (smootherCode.compare(Jacobi::Code()) == 0)
			return new JacobiSmoother(nSmootherIterations);
		if (smootherCode.compare(BlockJacobi::Code()) == 0)
			return new BlockJacobiSmoother(blockSize, nSmootherIterations);
		if (smootherCode.compare(BlockDampedJacobi23::Code()) == 0)
			return new BlockDampedJacobi23Smoother(blockSize, nSmootherIterations);

		Utils::FatalError("Unmanaged smoother code: '" + smootherCode + "'");
	}
};