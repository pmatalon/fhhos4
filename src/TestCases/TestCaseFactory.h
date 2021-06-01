#pragma once
#include "DefaultTestCase.h"
#include "1D/SegmentTestCase.h"
#include "1D/Heterogeneity1DTestCase.h"
#include "2D/SquareTestCase.h"
#include "2D/KelloggTestCase.h"
#include "2D/BarWith4HolesTestCase.h"
#include "2D/SquareCircleTestCase.h"
#include "2D/SquareCornerSquareTestCase.h"
#include "2D/SquareCenterSquareTestCase.h"
#include "2D/SquareHolesTestCase.h"
#include "2D/RoomWithWallTestCase.h"
#include "2D/EDFTestCase.h"
#include "2D/HybridMeshTestCase.h"
#include "2D/MagnetismTestCase.h"
#include "3D/CubeTestCase.h"
#include "3D/PlateWith4HolesTestCase.h"
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
	if (pb.TestCaseCode.compare("segment") == 0)
		return new SegmentTestCase(pb);
	if (pb.TestCaseCode.compare("heterog") == 0)
		return new Heterogeneity1DTestCase(pb);

	Utils::FatalError("Test case '" + pb.TestCaseCode + "' is unknown or not implemented in 1D. Check -tc argument or use '-tc default'.");
}

template <>
TestCase<2>* TestCaseFactory<2>::Create(ProblemArguments pb)
{
	if (pb.TestCaseCode.compare("default") == 0)
		return new DefaultTestCase<2>(pb);
	if (pb.TestCaseCode.compare("square") == 0)
		return new SquareTestCase(pb);
	if (pb.TestCaseCode.compare("square4quadrants") == 0)
		return new SquareTestCase(pb);
	if (pb.TestCaseCode.compare("kellogg") == 0)
		return new KelloggTestCase(pb);
	if (pb.TestCaseCode.compare("barwith4holes") == 0)
		return new BarWith4HolesTestCase(pb);
	if (pb.TestCaseCode.compare("squarecircle") == 0)
		return new SquareCircleTestCase(pb);
	if (pb.TestCaseCode.compare("squarecornersquare") == 0)
		return new SquareCornerSquareTestCase(pb);
	if (pb.TestCaseCode.compare("squarecentersquare") == 0)
		return new SquareCenterSquareTestCase(pb);
	if (pb.TestCaseCode.compare("squareholes") == 0)
		return new SquareHolesTestCase(pb);
	if (pb.TestCaseCode.compare("edf") == 0)
		return new EDFTestCase(pb);
	if (pb.TestCaseCode.compare("hybridmesh") == 0)
		return new HybridMeshTestCase(pb);
	if (pb.TestCaseCode.compare("magnetism") == 0)
		return new MagnetismTestCase(pb);
	if (pb.TestCaseCode.compare("roomwithwall") == 0)
		return new RoomWithWallTestCase(pb);

	Utils::FatalError("Test case '" + pb.TestCaseCode + "' is unknown or not implemented in 2D. Check -tc argument or use '-tc default'.");
}

template <>
TestCase<3>* TestCaseFactory<3>::Create(ProblemArguments pb)
{
	if (pb.TestCaseCode.compare("default") == 0)
		return new DefaultTestCase<3>(pb);
	if (pb.TestCaseCode.compare("cube") == 0)
		return new CubeTestCase(pb);
	if (pb.TestCaseCode.compare("platewith4holes") == 0)
		return new PlateWith4HolesTestCase(pb);

	Utils::FatalError("Test case '" + pb.TestCaseCode + "' is unknown or not implemented in 3D. Check -tc argument or use '-tc default'.");
}
