// Validation tests for "High-order multigrid strategies for HHO discretizations of elliptic
// equations" (Di Pietro, Matalon, Mycek, Rude, Numer. Linear Algebra Appl. 2022),
// reproducibility/2022_high_order_strategies.md.
#include <gtest/gtest.h>
#include <string>
#include <vector>
#include "support/RunHelper.h"

using namespace fhhos4_tests;

namespace
{
	// Section 3.4.1: a locally-refined mesh on the square, driven through the generic "square"
	// sine test case (the .geo file is not one of the special-cased geometries, so it's routed
	// via -tc square rather than by its own geometry code).
	ProgramArguments LocalRefArgs(int elemBasisOrthogonalizeCode, int n, const std::string& solverCode)
	{
		ProgramArguments args;
		args.Problem.GeoCode = "square4quadrants_tri_localref";
		args.Problem.TestCaseCode = "square";
		args.Discretization.N = n;
		args.Discretization.PolyDegree = 2; // k = 1
		args.Discretization.OrthogonalizeElemBasesCode = elemBasisOrthogonalizeCode;
		args.Solver.MG.H_CS = H_CoarsStgy::GMSHSplittingRefinement; // -cs r
		args.Solver.SolverCode = solverCode;
		args.Actions.UseCache = false; // -no-cache
		return args;
	}

	// Mirrors main.cpp's OPT_HPConfig switch (src/main.cpp) for -hp-config 1..4.
	void ApplyHPConfig(ProgramArguments& args, int code)
	{
		switch (code)
		{
		case 1:
			args.Solver.MG.HP_CS = HP_CoarsStgy::H_only;
			break;
		case 2:
			args.Solver.MG.HP_CS = HP_CoarsStgy::P_then_H;
			args.Solver.MG.P_CS = P_CoarsStgy::Minus2;
			args.Solver.MG.GMG_P_Prolong = GMG_P_Prolongation::Injection;
			args.Solver.MG.GMG_P_Restrict = GMG_P_Restriction::RemoveHigherOrders;
			break;
		case 3:
			args.Solver.MG.HP_CS = HP_CoarsStgy::P_then_H;
			args.Solver.MG.P_CS = P_CoarsStgy::Minus2;
			args.Solver.MG.GMG_P_Prolong = GMG_P_Prolongation::H_Prolongation;
			args.Solver.MG.GMG_P_Restrict = GMG_P_Restriction::P_Transpose;
			break;
		case 4:
			args.Solver.MG.HP_CS = HP_CoarsStgy::HP_then_H;
			args.Solver.MG.P_CS = P_CoarsStgy::Minus1;
			break;
		default:
			throw std::runtime_error("ApplyHPConfig: unsupported hp-config code " + std::to_string(code));
		}
	}

	ProgramArguments SquareCartArgs(int k, int n, const std::string& solverCode, double tolerance)
	{
		ProgramArguments args;
		args.Problem.GeoCode = "square";
		args.Discretization.MeshCode = "cart";
		args.Discretization.N = n;
		args.Discretization.PolyDegree = k + 1;
		args.Solver.MG.H_CS = H_CoarsStgy::GMSHSplittingRefinement; // -cs r
		args.Solver.SolverCode = solverCode;
		args.Solver.Tolerance = tolerance;
		return args;
	}
}

// Section 3.4.1 (basis normalization): with local refinement, orthonormalized element bases
// (-e-ogb 3) make the multigrid diverge, while orthogonalization without normalization
// (-e-ogb 1) converges. Observed at n=32: divergence vs. 13 iterations.
//   fhhos4 -geo square4quadrants_tri_localref -no-cache -tc square -cs r -k 1 -n 32 -e-ogb {3|1}
TEST(HpStrategies2022, BasisNormalization_OrthonormalDiverges)
{
	ProgramArguments args = LocalRefArgs(/*elemBasisOrthogonalizeCode*/ 3, /*n*/ 32, "mg");
	EXPECT_EXIT(RunDiffusionHHO(args), ::testing::ExitedWithCode(EXIT_FAILURE), "");
}

TEST(HpStrategies2022, BasisNormalization_OrthogonalConverges)
{
	ProgramArguments args = LocalRefArgs(/*elemBasisOrthogonalizeCode*/ 1, /*n*/ 32, "mg");
	ProgramResults results = RunDiffusionHHO(args);
	EXPECT_GT(results.IterationCount, 0);
	EXPECT_LE(results.IterationCount, 30);
}

// Figures 7-18 (a representative slice): every hp-multigrid strategy (1-4), used both as a
// stand-alone solver and as an FCG preconditioner, should converge.
class HPConfigTest : public ::testing::TestWithParam<std::tuple<std::string, int>>
{
};

TEST_P(HPConfigTest, AllConverge)
{
	auto [solverCode, hpConfig] = GetParam();

	ProgramArguments args = SquareCartArgs(/*k*/ 5, /*n*/ 32, solverCode, /*tolerance*/ 1e-10);
	ApplyHPConfig(args, hpConfig);

	ProgramResults results = RunDiffusionHHO(args);
	EXPECT_GT(results.IterationCount, 0);
	EXPECT_LE(results.IterationCount, 100);
}

INSTANTIATE_TEST_SUITE_P(SquareCart, HPConfigTest,
	::testing::Combine(::testing::Values(std::string("mg"), std::string("fcgmg")), ::testing::Values(1, 2, 3, 4)));

// Figure 3: at high order (hp-config 2, tight tolerance), the L2 error should still follow the
// theoretical h^(k+2) convergence order.
class ConvergenceOrderHighOrderTest : public ::testing::TestWithParam<int>
{
};

TEST_P(ConvergenceOrderHighOrderTest, MatchesTheoreticalOrder)
{
	int k = GetParam();
	std::vector<int> ns = { 16, 32 };

	std::vector<double> h;
	std::vector<double> errors;
	for (int n : ns)
	{
		ProgramArguments args = SquareCartArgs(k, n, "fcgmg", /*tolerance*/ 1e-12);
		ApplyHPConfig(args, /*hp-config*/ 2);

		ProgramResults results = RunDiffusionHHO(args);
		ASSERT_GT(results.L2Error, 0);
		h.push_back(1.0 / n);
		errors.push_back(results.L2Error);
	}

	double order = EstimateConvergenceOrder(h, errors);
	EXPECT_NEAR(order, static_cast<double>(k + 2), 0.4);
}

INSTANTIATE_TEST_SUITE_P(SquareCart, ConvergenceOrderHighOrderTest, ::testing::Values(2, 3, 4, 5));
