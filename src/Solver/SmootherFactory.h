#pragma once
#include "Smoother.h"
using namespace std;

class SmootherFactory
{

public:
	static Smoother* Create(string smootherCode, int nSmootherIterations, int blockSize, double omega)
	{
		if (smootherCode.compare("sor") == 0 || smootherCode.compare("gs") == 0)
			return new BlockSORSmoother(1, omega, Direction::Forward, nSmootherIterations);
		if (smootherCode.compare("rsor") == 0 || smootherCode.compare("rgs") == 0)
			return new BlockSORSmoother(1, omega, Direction::Backward, nSmootherIterations);

		if (smootherCode.compare("bsor") == 0 || smootherCode.compare("bgs") == 0)
			return new BlockSORSmoother(blockSize, omega, Direction::Forward, nSmootherIterations);
		if (smootherCode.compare("rbsor") == 0 || smootherCode.compare("rbgs") == 0)
			return new BlockSORSmoother(blockSize, omega, Direction::Backward, nSmootherIterations);

		if (smootherCode.compare("sbsor") == 0 || smootherCode.compare("sbgs") == 0)
			return new BlockSORSmoother(blockSize, omega, Direction::Symmetric, nSmootherIterations);

		if (smootherCode.compare("j") == 0)
			return new BlockJacobiSmoother(1, omega, nSmootherIterations);
		if (smootherCode.compare("bj") == 0)
			return new BlockJacobiSmoother(blockSize, omega, nSmootherIterations);
		if (smootherCode.compare("bj23") == 0)
			return new BlockJacobiSmoother(blockSize, 2.0/3.0, nSmootherIterations);

		Utils::FatalError("Unmanaged smoother code: '" + smootherCode + "'");
		return nullptr;
	}
};