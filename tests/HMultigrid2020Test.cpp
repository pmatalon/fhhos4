// Validation tests for "An h-multigrid method for Hybrid High-Order discretizations"
// (Di Pietro, Hulsemann, Matalon, Mycek, Rude, Ruiz, SIAM J. Sci. Comput. 2021),
// reproducibility/2020_MG_for_HHO.md.
#include <gtest/gtest.h>
#include <algorithm>
#include <string>
#include <tuple>
#include <vector>
#include "support/RunHelper.h"

using namespace fhhos4_tests;

namespace
{
	ProgramArguments SquareVCycleArgs(const std::string& meshCode, int k, int n)
	{
		ProgramArguments args;
		args.Problem.GeoCode = "square";
		args.Discretization.Mesher = "inhouse";
		args.Discretization.MeshCode = meshCode;
		args.Discretization.N = n;
		args.Discretization.PolyDegree = k + 1;
		args.Solver.SolverCode = "mg";
		args.Solver.MG.CycleLetter = 'V';
		args.Solver.MG.PreSmoothingIterations = 1;
		args.Solver.MG.PostSmoothingIterations = 1;
		return args;
	}

	const std::vector<int> QuickMeshSizes = { 8, 16, 32, 64 };

	// n=8/16 sit in a degenerate pre-asymptotic regime (the multigrid hierarchy barely has any
	// levels, so iteration counts are artificially tiny, e.g. 1) - not representative of the
	// mesh-independence claim. The paper's own Figure 4.1 range for this case starts at n=32.
	const std::vector<int> MeshIndependenceSizes = { 32, 64, 128 };
}

// Figure 4.1: L2-error convergence order should match the theoretical rate the code itself
// encodes in Diffusion_HHO::AssertSchemeConvergence (h^2 for k=0, h^(k+2) for k>=1).
class ConvergenceOrderTest : public ::testing::TestWithParam<std::tuple<std::string, int>>
{
};

TEST_P(ConvergenceOrderTest, MatchesTheoreticalOrder)
{
	auto [meshCode, k] = GetParam();

	std::vector<double> h;
	std::vector<double> errors;
	for (int n : QuickMeshSizes)
	{
		ProgramResults results = RunDiffusionHHO(SquareVCycleArgs(meshCode, k, n));
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

// Figure 4.1: the V(1,1)-cycle iteration count should be essentially mesh-independent.
class MeshIndependenceTest : public ::testing::TestWithParam<std::tuple<std::string, int>>
{
};

TEST_P(MeshIndependenceTest, IterationCountStaysBounded)
{
	auto [meshCode, k] = GetParam();

	std::vector<int> iterationCounts;
	for (int n : MeshIndependenceSizes)
	{
		ProgramResults results = RunDiffusionHHO(SquareVCycleArgs(meshCode, k, n));
		ASSERT_GT(results.IterationCount, 0);
		iterationCounts.push_back(results.IterationCount);
	}

	int smallest = iterationCounts.front();
	int largest = iterationCounts.back();
	// Real V(1,1) behavior has a mild pre-asymptotic increase before leveling off (e.g. observed
	// 19 -> 24 -> 27 for the unstructured mesh at k=0), not a hard plateau - bound the growth
	// ratio generously rather than requiring near-exact constancy.
	EXPECT_LE(largest, smallest * 2) << "iteration count should stay roughly bounded across mesh refinement";
	EXPECT_LE(largest, 40);
}

INSTANTIATE_TEST_SUITE_P(Square, MeshIndependenceTest,
	::testing::Combine(::testing::Values(std::string("cart"), std::string("stri")), ::testing::Values(0, 1, 2, 3)));

// Figure 4.1 explicitly documents this configuration as diverging (as of release 1.0). On the
// current code it still fails, though the failure mode has drifted: mesh construction now hits
// "Unmanaged refinement strategy" before the solver even runs, rather than iterative divergence.
// Either way this is a known-broken configuration; regression-documents it so a future fix must
// touch this test deliberately rather than silently change behavior unnoticed.
TEST(HMultigrid2020, FailsAsDocumented_StructuredTetra_K0)
{
	ProgramArguments args;
	args.Problem.GeoCode = "cube";
	args.Discretization.Mesher = "inhouse";
	args.Discretization.MeshCode = "stetra";
	args.Discretization.N = 16;
	args.Discretization.PolyDegree = 1; // k = 0
	args.Solver.SolverCode = "mg";
	args.Solver.MG.CycleLetter = 'V';
	args.Solver.MG.PreSmoothingIterations = 2;
	args.Solver.MG.PostSmoothingIterations = 2;

	// Utils::FatalError prints to stdout (not stderr), so GTest's death-test output-matching
	// (which checks stderr) can't see the message; just check the exit code.
	EXPECT_EXIT(RunDiffusionHHO(args), ::testing::ExitedWithCode(EXIT_FAILURE), "");
}

// Figure 4.7: robustness on the heterogeneous Kellogg benchmark.
class KelloggTest : public ::testing::TestWithParam<int>
{
};

TEST_P(KelloggTest, ConvergesWithinBound)
{
	int k = GetParam();

	ProgramArguments args;
	args.Problem.GeoCode = "square4quadrants";
	args.Problem.TestCaseCode = "kellogg";
	args.Discretization.Mesher = "inhouse";
	args.Discretization.MeshCode = "cart";
	args.Discretization.N = 256;
	args.Discretization.PolyDegree = k + 1;
	args.Solver.SolverCode = "mg";
	args.Solver.MG.CycleLetter = 'V';
	args.Solver.MG.PreSmoothingIterations = 1;
	args.Solver.MG.PostSmoothingIterations = 1;

	ProgramResults results = RunDiffusionHHO(args);
	EXPECT_GT(results.IterationCount, 0);
	EXPECT_LE(results.IterationCount, 30);
}

INSTANTIATE_TEST_SUITE_P(SquareFourQuadrants, KelloggTest, ::testing::Values(0, 1, 2, 3));

// Figure 4.9: iteration count should stay roughly bounded as the heterogeneity ratio grows
// from 1e0 to 1e8, with the heterogeneous weighting enabled (-g 1).
class HeterogeneityRatioSweepTest : public ::testing::TestWithParam<int>
{
};

TEST_P(HeterogeneityRatioSweepTest, IterationCountRatioBounded)
{
	int k = GetParam();
	std::vector<double> ratios = { 1e0, 1e2, 1e4, 1e6, 1e8 };

	std::vector<int> iterationCounts;
	for (double ratio : ratios)
	{
		ProgramArguments args;
		args.Problem.GeoCode = "square4quadrants";
		args.Problem.HeterogeneityRatio = ratio;
		args.Discretization.Mesher = "inhouse";
		args.Discretization.MeshCode = "cart";
		args.Discretization.N = 64;
		args.Discretization.PolyDegree = k + 1;
		args.Solver.SolverCode = "mg";
		args.Solver.MG.UseGalerkinOperator = true;
		args.Solver.MG.CycleLetter = 'V';
		args.Solver.MG.PreSmoothingIterations = 0;
		args.Solver.MG.PostSmoothingIterations = 3;

		ProgramResults results = RunDiffusionHHO(args);
		ASSERT_GT(results.IterationCount, 0);
		iterationCounts.push_back(results.IterationCount);
	}

	int minIterations = *std::min_element(iterationCounts.begin(), iterationCounts.end());
	int maxIterations = *std::max_element(iterationCounts.begin(), iterationCounts.end());
	EXPECT_LE(maxIterations, minIterations * 3) << "iteration count should stay roughly bounded across the heterogeneity sweep";
}

INSTANTIATE_TEST_SUITE_P(SquareFourQuadrants, HeterogeneityRatioSweepTest, ::testing::Values(0, 1, 2, 3));
