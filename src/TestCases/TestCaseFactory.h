#pragma once
#include "1D/SineSolution1DTestCase.h"
#include "1D/PolySolution1DTestCase.h"
#include "1D/Heterogeneity1DTestCase.h"
#include "2D/SineSolution2DTestCase.h"
#include "2D/PolySolution2DTestCase.h"
#include "2D/ExpSolution2DTestCase.h"
#include "2D/ZeroSolution2DTestCase.h"
#include "2D/OneSolution2DTestCase.h"
#include "2D/XSolution2DTestCase.h"
#include "2D/KelloggTestCase.h"
#include "3D/SineSolution3DTestCase.h"
#include "3D/PolySolution3DTestCase.h"
#include "3D/ExpSolution3DTestCase.h"
using namespace std;

template <int Dim>
class TestCaseFactory
{
public:
	static TestCase<Dim>* Create(string tcCode, DiffusionPartition<Dim>* diffusionPartition, string bcCode) { assert(false); }
};

template <>
TestCase<1>* TestCaseFactory<1>::Create(string tcCode, DiffusionPartition<1>* diffusionPartition, string bcCode)
{
	if (tcCode.compare("sine") == 0)
		return new SineSolution1DTestCase(diffusionPartition, bcCode);
	if (tcCode.compare("poly") == 0)
		return new PolySolution1DTestCase(diffusionPartition, bcCode);
	if (tcCode.compare("heterog") == 0)
		return new Heterogeneity1DTestCase(diffusionPartition, bcCode);

	Utils::FatalError("Test case '" + tcCode + "' is unknown or not implemented in 1D. Check -rhs argument.");
}

template <>
TestCase<2>* TestCaseFactory<2>::Create(string tcCode, DiffusionPartition<2>* diffusionPartition, string bcCode)
{
	if (tcCode.compare("sine") == 0)
		return new SineSolution2DTestCase(diffusionPartition, bcCode);
	if (tcCode.compare("poly") == 0)
		return new PolySolution2DTestCase(diffusionPartition, bcCode);
	if (tcCode.compare("exp") == 0)
		return new ExpSolution2DTestCase(diffusionPartition, bcCode);
	if (tcCode.compare("zero") == 0)
		return new ZeroSolution2DTestCase(diffusionPartition, bcCode);
	if (tcCode.compare("one") == 0)
		return new OneSolution2DTestCase(diffusionPartition, bcCode);
	if (tcCode.compare("x") == 0)
		return new XSolution2DTestCase(diffusionPartition, bcCode);
	if (tcCode.compare("kellogg") == 0)
		return new KelloggTestCase(diffusionPartition, bcCode);

	Utils::FatalError("Test case '" + tcCode + "' is unknown or not implemented in 2D. Check -rhs argument.");
}

template <>
TestCase<3>* TestCaseFactory<3>::Create(string tcCode, DiffusionPartition<3>* diffusionPartition, string bcCode)
{
	if (tcCode.compare("sine") == 0)
		return new SineSolution3DTestCase(diffusionPartition, bcCode);
	if (tcCode.compare("poly") == 0)
		return new PolySolution3DTestCase(diffusionPartition, bcCode);
	if (tcCode.compare("exp") == 0)
		return new ExpSolution3DTestCase(diffusionPartition, bcCode);

	Utils::FatalError("Test case '" + tcCode + "' is unknown or not implemented in 3D. Check -rhs argument.");
}
