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
	static TestCase<Dim>* Create(string tcCode, DiffusionField<Dim>* diffusionField, string bcCode) { assert(false); }
};

template <>
TestCase<1>* TestCaseFactory<1>::Create(string tcCode, DiffusionField<1>* diffusionField, string bcCode)
{
	if (tcCode.compare("sine") == 0)
		return new SineSolution1DTestCase(diffusionField, bcCode);
	if (tcCode.compare("poly") == 0)
		return new PolySolution1DTestCase(diffusionField, bcCode);
	if (tcCode.compare("heterog") == 0)
		return new Heterogeneity1DTestCase(diffusionField, bcCode);

	Utils::FatalError("Test case '" + tcCode + "' is unknown or not implemented in 1D. Check -rhs argument.");
}

template <>
TestCase<2>* TestCaseFactory<2>::Create(string tcCode, DiffusionField<2>* diffusionField, string bcCode)
{
	if (tcCode.compare("sine") == 0)
		return new SineSolution2DTestCase(diffusionField, bcCode);
	if (tcCode.compare("poly") == 0)
		return new PolySolution2DTestCase(diffusionField, bcCode);
	if (tcCode.compare("exp") == 0)
		return new ExpSolution2DTestCase(diffusionField, bcCode);
	if (tcCode.compare("zero") == 0)
		return new ZeroSolution2DTestCase(diffusionField, bcCode);
	if (tcCode.compare("one") == 0)
		return new OneSolution2DTestCase(diffusionField, bcCode);
	if (tcCode.compare("x") == 0)
		return new XSolution2DTestCase(diffusionField, bcCode);
	if (tcCode.compare("kellogg") == 0)
		return new KelloggTestCase(diffusionField, bcCode);

	Utils::FatalError("Test case '" + tcCode + "' is unknown or not implemented in 2D. Check -rhs argument.");
}

template <>
TestCase<3>* TestCaseFactory<3>::Create(string tcCode, DiffusionField<3>* diffusionField, string bcCode)
{
	if (tcCode.compare("sine") == 0)
		return new SineSolution3DTestCase(diffusionField, bcCode);
	if (tcCode.compare("poly") == 0)
		return new PolySolution3DTestCase(diffusionField, bcCode);
	if (tcCode.compare("exp") == 0)
		return new ExpSolution3DTestCase(diffusionField, bcCode);

	Utils::FatalError("Test case '" + tcCode + "' is unknown or not implemented in 3D. Check -rhs argument.");
}
