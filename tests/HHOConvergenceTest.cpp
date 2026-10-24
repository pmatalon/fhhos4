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

//------------------------------------------- The space Q -------------------------------------------//

// HHO in Q on Cartesian meshes (-poly-space q: cells and faces Q^k, reconstruction Q^k + P^{k+1}): the orders of P.
// Before 2026-10-23 (reconstruction in Q^{k+1}, bases of Q not ordered by degree), it didn't converge for k >= 1.
//   ./bin/fhhos4 -geo square -mesh cart -mesher inhouse -poly-space q -s ch -k {0|1|2|3} -n {8|16|32}
class ConvergenceOrderQTest : public ::testing::TestWithParam<int>
{
};

TEST_P(ConvergenceOrderQTest, MatchesTheoreticalOrder)
{
	int k = GetParam();

	std::vector<double> h;
	std::vector<double> errors;
	for (int n : MeshSizes)
	{
		ProgramArguments args = SquareArgs("cart", k, n);
		args.Discretization.UsePolynomialSpaceQ = true;
		ProgramResults results = RunDiffusionHHO(args);
		ASSERT_GT(results.L2Error, 0);
		h.push_back(1.0 / n);
		errors.push_back(results.L2Error);
	}

	double order = EstimateConvergenceOrder(h, errors);
	double expectedOrder = (k == 0) ? 2.0 : static_cast<double>(k + 2);
	EXPECT_NEAR(order, expectedOrder, 0.3);
}

INSTANTIATE_TEST_SUITE_P(Square, ConvergenceOrderQTest, ::testing::Values(0, 1, 2, 3));
#endif // ENABLE_2D

using namespace fhhos4_tests;

// The solution of -source poly, x(1-x)y(1-y) (times z(1-z) in 3D), is in Q^2: HHO in Q reproduces it at k = 2, to
// rounding (the reconstruction of the interpolate of a function of the reconstruction space is the function, and the
// stabilization vanishes on it); also with the bases orthogonalized on each cell (-e-basis monomials -e-ogb 1: the cell
// basis made of the first functions of the reconstruction basis). In P, it needs k + 1 >= 2d.
//   ./bin/fhhos4 -geo {square|cube} -mesh cart -mesher inhouse -poly-space q -source poly -s ch -k 2 -n 2 [-e-basis monomials -e-ogb 1]
class ExactnessQTest : public ::testing::TestWithParam<std::tuple<std::string, bool>>
{
};

TEST_P(ExactnessQTest, ReproducesASolutionOfQ2)
{
	auto [geo, orthogonalizedMonomials] = GetParam();
	ProgramArguments args;
	args.Problem.GeoCode = geo;
	args.Problem.SourceCode = "poly";
	args.Discretization.Mesher = "inhouse";
	args.Discretization.MeshCode = "cart";
	args.Discretization.N = 2; // exact on any mesh: the smallest with interior faces (the orthogonalization is long in 3D)
	args.Discretization.PolyDegree = 3;
	args.Discretization.UsePolynomialSpaceQ = true;
	if (orthogonalizedMonomials)
	{
		args.Discretization.ElemBasisCode = "monomials";
		args.Discretization.OrthogonalizeElemBasesCode = 1;
	}
	args.Solver.SolverCode = "ch";
	ProgramResults results = RunDiffusionHHO(args);
	EXPECT_LT(results.L2Error, 1e-12);
}

#ifdef ENABLE_2D
INSTANTIATE_TEST_SUITE_P(Square, ExactnessQTest, ::testing::Combine(::testing::Values(std::string("square")), ::testing::Bool()));
#endif // ENABLE_2D

#ifdef ENABLE_3D
INSTANTIATE_TEST_SUITE_P(Cube, ExactnessQTest, ::testing::Combine(::testing::Values(std::string("cube")), ::testing::Bool()));
#endif // ENABLE_3D
