#pragma once
#include "DefaultTestCase.h"
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
#include "2D/BarWith4HolesTestCase.h"
#include "2D/SquareCircleTestCase.h"
#include "2D/SquareHolesTestCase.h"
#include "2D/EDFTestCase.h"
#include "3D/SineSolution3DTestCase.h"
#include "3D/PolySolution3DTestCase.h"
#include "3D/ExpSolution3DTestCase.h"
using namespace std;

template <int Dim>
class TestCaseFactory
{
public:
	static TestCase<Dim>* Create(ProblemArguments pb) { assert(false); }
};

template <>
TestCase<1>* TestCaseFactory<1>::Create(ProblemArguments pb)
{
	if (pb.TestCaseCode.compare("default") == 0)
		return new DefaultTestCase<1>(pb);
	if (pb.TestCaseCode.compare("sine") == 0)
		return new SineSolution1DTestCase(pb);
	if (pb.TestCaseCode.compare("poly") == 0)
		return new PolySolution1DTestCase(pb);
	if (pb.TestCaseCode.compare("heterog") == 0)
		return new Heterogeneity1DTestCase(pb);

	Utils::FatalError("Test case '" + pb.TestCaseCode + "' is unknown or not implemented in 3D. Check -tc argument.");
}

template <>
TestCase<2>* TestCaseFactory<2>::Create(ProblemArguments pb)
{
	if (pb.TestCaseCode.compare("default") == 0)
		return new DefaultTestCase<2>(pb);
	if (pb.TestCaseCode.compare("sine") == 0)
		return new SineSolution2DTestCase(pb);
	if (pb.TestCaseCode.compare("poly") == 0)
		return new PolySolution2DTestCase(pb);
	if (pb.TestCaseCode.compare("exp") == 0)
		return new ExpSolution2DTestCase(pb);
	if (pb.TestCaseCode.compare("zero") == 0)
		return new ZeroSolution2DTestCase(pb);
	if (pb.TestCaseCode.compare("one") == 0)
		return new OneSolution2DTestCase(pb);
	if (pb.TestCaseCode.compare("x") == 0)
		return new XSolution2DTestCase(pb);
	if (pb.TestCaseCode.compare("kellogg") == 0)
		return new KelloggTestCase(pb);
	if (pb.TestCaseCode.compare("barwith4holes") == 0)
		return new BarWith4HolesTestCase(pb);
	if (pb.TestCaseCode.compare("squarecircle") == 0)
		return new SquareCircleTestCase(pb);
	if (pb.TestCaseCode.compare("squareholes") == 0)
		return new SquareHolesTestCase(pb);
	if (pb.TestCaseCode.compare("edf") == 0)
		return new EDFTestCase(pb);

	Utils::FatalError("Test case '" + pb.TestCaseCode + "' is unknown or not implemented in 3D. Check -tc argument.");
}

template <>
TestCase<3>* TestCaseFactory<3>::Create(ProblemArguments pb)
{
	if (pb.TestCaseCode.compare("default") == 0)
		return new DefaultTestCase<3>(pb);
	if (pb.TestCaseCode.compare("sine") == 0)
		return new SineSolution3DTestCase(pb);
	if (pb.TestCaseCode.compare("poly") == 0)
		return new PolySolution3DTestCase(pb);
	if (pb.TestCaseCode.compare("exp") == 0)
		return new ExpSolution3DTestCase(pb);

	Utils::FatalError("Test case '" + pb.TestCaseCode + "' is unknown or not implemented in 3D. Check -tc argument.");
}
