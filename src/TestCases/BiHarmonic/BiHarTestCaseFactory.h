#pragma once
#include "BiHarDefaultTestCase.h"
#ifdef ENABLE_2D
	#include "2D/SquareBiHarTestCase.h"
#endif
#ifdef ENABLE_3D
	#include "3D/CubeBiHarTestCase.h"
#endif
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


#ifdef ENABLE_2D
template <>
BiHarmonicTestCase<2>* BiHarTestCaseFactory<2>::Create(ProblemArguments pb)
{
	if (pb.TestCaseCode.compare("default") == 0)
		return new BiHarDefaultTestCase<2>(pb);
	if (pb.TestCaseCode.compare("square") == 0)
		return new SquareBiHarTestCase(pb);

	Utils::FatalError("Test case '" + pb.TestCaseCode + "' is unknown or not implemented in 2D. Check -tc argument or use '-tc default'.");
	return nullptr;
}
#endif // ENABLE_2D


#ifdef ENABLE_3D
template <>
BiHarmonicTestCase<3>* BiHarTestCaseFactory<3>::Create(ProblemArguments pb)
{
	if (pb.TestCaseCode.compare("default") == 0)
		return new BiHarDefaultTestCase<3>(pb);
	if (pb.TestCaseCode.compare("cube") == 0)
		return new CubeBiHarTestCase(pb);

	Utils::FatalError("Test case '" + pb.TestCaseCode + "' is unknown or not implemented in 3D. Check -tc argument or use '-tc default'.");
	return nullptr;
}
#endif // ENABLE_3D
