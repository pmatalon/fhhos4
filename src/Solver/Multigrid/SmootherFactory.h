#pragma once
#include "Smoother.h"
using namespace std;

class SmootherFactory
{

public:
	static Smoother* Create(string smootherCode, int nSmootherIterations, int blockSize, double omega)
	{
		if (smootherCode.compare("gs") == 0)
			return new GaussSeidelSmoother(Direction::Forward, nSmootherIterations);
		if (smootherCode.compare("rgs") == 0)
			return new GaussSeidelSmoother(Direction::Backward, nSmootherIterations);

		if (smootherCode.compare("sor") == 0)
		{
			if (omega == 1)
				return new GaussSeidelSmoother(Direction::Forward, nSmootherIterations);
			else
				return new BlockSORSmoother(1, omega, Direction::Forward, nSmootherIterations);
		}
		if (smootherCode.compare("rsor") == 0)
		{
			if (omega == 1)
				return new GaussSeidelSmoother(Direction::Backward, nSmootherIterations);
			else
				return new BlockSORSmoother(1, omega, Direction::Backward, nSmootherIterations);
		}

		if (smootherCode.compare("bsor") == 0 || smootherCode.compare("bgs") == 0)
		{
			if (blockSize == 1 && omega == 1)
				return new GaussSeidelSmoother(Direction::Forward, nSmootherIterations);
			else if (omega == 1)
				return new BlockGaussSeidelSmoother(blockSize, Direction::Forward, false, nSmootherIterations);
			else
				return new BlockSORSmoother(blockSize, omega, Direction::Forward, nSmootherIterations);
		}
		if (smootherCode.compare("rbsor") == 0 || smootherCode.compare("rbgs") == 0)
		{
			if (blockSize == 1 && omega == 1)
				return new GaussSeidelSmoother(Direction::Backward, nSmootherIterations);
			else if (omega == 1)
				return new BlockGaussSeidelSmoother(blockSize, Direction::Backward, false, nSmootherIterations);
			else
				return new BlockSORSmoother(blockSize, omega, Direction::Backward, nSmootherIterations);
		}
		if (smootherCode.compare("hbgs") == 0)
		{
			return new BlockGaussSeidelSmoother(blockSize, Direction::Forward, true, nSmootherIterations);
		}
		if (smootherCode.compare("hrbgs") == 0)
		{
			return new BlockGaussSeidelSmoother(blockSize, Direction::Backward, true, nSmootherIterations);
		}


		if (smootherCode.compare("sbsor") == 0 || smootherCode.compare("sbgs") == 0)
		{
			if (blockSize == 1 && omega == 1)
				return new GaussSeidelSmoother(Direction::Symmetric, nSmootherIterations);
			else if (omega == 1)
				return new BlockGaussSeidelSmoother(blockSize, Direction::Symmetric, false, nSmootherIterations);
			else
				return new BlockSORSmoother(blockSize, omega, Direction::Symmetric, nSmootherIterations);
		}

		if (smootherCode.compare("j") == 0)
			return new BlockJacobiSmoother(1, omega, nSmootherIterations);
		if (smootherCode.compare("j23") == 0)
			return new BlockJacobiSmoother(1, 2.0/3.0, nSmootherIterations);
		if (smootherCode.compare("bj") == 0)
			return new BlockJacobiSmoother(blockSize, omega, nSmootherIterations);
		if (smootherCode.compare("bj23") == 0)
			return new BlockJacobiSmoother(blockSize, 2.0/3.0, nSmootherIterations);

		Utils::FatalError("Unmanaged smoother code: '" + smootherCode + "'");
		return nullptr;
	}
};