#pragma once
#include "2D/SquareBiHarDDTestCase.h"
using namespace std;

template <int Dim>
class BiHarDDTestCaseFactory
{
public:
	static BiHarmonicDDTestCase<Dim>* Create(ProblemArguments pb) { assert(false); }
};

template <>
BiHarmonicDDTestCase<1>* BiHarDDTestCaseFactory<1>::Create(ProblemArguments pb)
{
	Utils::FatalError("Test case '" + pb.TestCaseCode + "' is unknown or not implemented in 1D. Check -tc argument.");
	return nullptr;
}

template <>
BiHarmonicDDTestCase<2>* BiHarDDTestCaseFactory<2>::Create(ProblemArguments pb)
{
	if (pb.TestCaseCode.compare("square") == 0)
		return new SquareBiHarDDTestCase(pb);

	Utils::FatalError("Test case '" + pb.TestCaseCode + "' is unknown or not implemented in 2D. Check -tc argument.");
	return nullptr;
}

template <>
BiHarmonicDDTestCase<3>* BiHarDDTestCaseFactory<3>::Create(ProblemArguments pb)
{
	Utils::FatalError("Test case '" + pb.TestCaseCode + "' is unknown or not implemented in 3D. Check -tc argument.");
	return nullptr;
}
