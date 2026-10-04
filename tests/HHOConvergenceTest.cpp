// A-priori convergence of the HHO discretization, independent of any solver study: the linear
// system is solved directly (sparse Cholesky) so that the L2 error is the discretization error.
#include <gtest/gtest.h>
#include <string>
#include <tuple>
#include <vector>
#include "support/RunHelper.h"

#ifdef ENABLE_2D
using namespace fhhos4_tests;

namespace
{
	ProgramArguments SquareArgs(const std::string& meshCode, int k, int n)
	{
		ProgramArguments args;
		args.Problem.GeoCode = "square";
		args.Discretization.Mesher = "inhouse";
		args.Discretization.MeshCode = meshCode;
		args.Discretization.N = n;
		args.Discretization.PolyDegree = k + 1;
		args.Solver.SolverCode = "ch";
		return args;
	}

	const std::vector<int> MeshSizes = { 8, 16, 32 };
}

// The L2-error convergence order should match the theoretical rate the code itself encodes in
// Diffusion_HHO::AssertSchemeConvergence (h^2 for k=0, h^(k+2) for k>=1).
class ConvergenceOrderTest : public ::testing::TestWithParam<std::tuple<std::string, int>>
{
};

//   ./bin/fhhos4 -geo square -mesh {cart|stri} -mesher inhouse -s ch -k {0|1|2|3} -n {8|16|32}
TEST_P(ConvergenceOrderTest, MatchesTheoreticalOrder)
{
	auto [meshCode, k] = GetParam();

	std::vector<double> h;
	std::vector<double> errors;
	for (int n : MeshSizes)
	{
		ProgramResults results = RunDiffusionHHO(SquareArgs(meshCode, k, n));
		ASSERT_GT(results.L2Error, 0);
		h.push_back(1.0 / n);
		errors.push_back(results.L2Error);
	}

	double order = EstimateConvergenceOrder(h, errors);
	double expectedOrder = (k == 0) ? 2.0 : static_cast<double>(k + 2);
	EXPECT_NEAR(order, expectedOrder, 0.3);
}

INSTANTIATE_TEST_SUITE_P(Square, ConvergenceOrderTest,
	::testing::Combine(::testing::Values(std::string("cart"), std::string("stri")), ::testing::Values(0, 1, 2, 3)));

// In-house quadrilateral mesh (trapezoids: -stretch 0.5 by default). Its elements had no physical
// part, and the assembly crashed. Only k=0: the quadrilateral elements use the polynomials of the
// reference square mapped by the bilinear map, which are not polynomials on a trapezoid, and the
// orders k>=1 don't converge at their theoretical rates.
//   ./bin/fhhos4 -geo square -mesh quad -mesher inhouse -s ch -k 0 -n {8|16|32}
INSTANTIATE_TEST_SUITE_P(SquareQuad, ConvergenceOrderTest,
	::testing::Combine(::testing::Values(std::string("quad")), ::testing::Values(0)));
#endif // ENABLE_2D
