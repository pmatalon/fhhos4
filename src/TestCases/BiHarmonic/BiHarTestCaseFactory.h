#pragma once
#include "2D/SquareBiHarTestCase.h"
using namespace std;

template <int Dim>
class BiHarTestCaseFactory
{
public:
	static BiHarmonicTestCase<Dim>* Create(ProblemArguments pb) { assert(false); }
};


#ifdef ENABLE_1D
template <>
BiHarmonicTestCase<1>* BiHarTestCaseFactory<1>::Create(ProblemArguments pb)
{
	Utils::FatalError("Test case '" + pb.TestCaseCode + "' is unknown or not implemented in 1D. Check -tc argument.");
	return nullptr;
}
#endif // ENABLE_1D

template <>
BiHarmonicTestCase<2>* BiHarTestCaseFactory<2>::Create(ProblemArguments pb)
{
	if (pb.TestCaseCode.compare("square") == 0)
		return new SquareBiHarTestCase(pb);

	Utils::FatalError("Test case '" + pb.TestCaseCode + "' is unknown or not implemented in 2D. Check -tc argument.");
	return nullptr;
}

#ifdef ENABLE_3D
template <>
BiHarmonicTestCase<3>* BiHarTestCaseFactory<3>::Create(ProblemArguments pb)
{
	Utils::FatalError("Test case '" + pb.TestCaseCode + "' is unknown or not implemented in 3D. Check -tc argument.");
	return nullptr;
}
#endif // ENABLE_3D
