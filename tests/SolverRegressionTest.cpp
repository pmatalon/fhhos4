// Regression tests for solver configurations that used to crash, and for the argument
// defaulting. Unlike the paper tests, they don't check reference iteration counts: they check
// that the runs complete and converge, on small in-house meshes.
#include <gtest/gtest.h>
#include <set>
#include <string>
#include <tuple>
#include <vector>
#include "support/RunHelper.h"

using namespace fhhos4_tests;

#ifdef ENABLE_2D
namespace
{
	ProgramArguments SquareCartArgs(int k, int n)
	{
		ProgramArguments args;
		args.Problem.GeoCode = "square";
		args.Discretization.Mesher = "inhouse";
		args.Discretization.MeshCode = "cart";
		args.Discretization.N = n;
		args.Discretization.PolyDegree = k + 1;
		return args;
	}
}

// With the K-cycle (and FCG preconditioned by a multigrid), the post-smoother must return Ax.
// The forward Gauss-Seidel (used for 'gs', and for 'bgs' when the block size is 1, i.e. k=0
// on faces in 2D) returned an empty Ax, which made FCG segfault. The default post-smoother
// ('rbgs', backward) was not affected. A small coarse-size gives several levels, so that the
// FCG of the K-cycle is actually used (with 2 levels, the coarse level is solved directly).
class ForwardPostSmootherTest : public ::testing::TestWithParam<std::string>
{
};

//   ./bin/fhhos4 -geo square -mesh cart -mesher inhouse -k 0 -n 64 -coarse-size 10 -smoothers gs,gs -s {mg -cycle K,1,1|fcgmg|uamg|aggregamg}
TEST_P(ForwardPostSmootherTest, Converges)
{
	std::string solverCode = GetParam();

	ProgramArguments args = SquareCartArgs(0, 64);
	args.Solver.SolverCode = solverCode;
	args.Solver.MG.MatrixMaxSizeForCoarsestLevel = 10;
	args.Solver.MG.PreSmootherCode = "gs";
	args.Solver.MG.PostSmootherCode = "gs";
	bool defaultCycle = true; // the AMGs default to the K-cycle
	if (solverCode == "mg")
	{
		args.Solver.MG.CycleLetter = 'K';
		args.Solver.MG.PreSmoothingIterations = 1;
		args.Solver.MG.PostSmoothingIterations = 1;
		defaultCycle = false;
	}

	ProgramResults results = RunDiffusionHHO(args, defaultCycle);
	EXPECT_GT(results.IterationCount, 0);
	EXPECT_LE(results.IterationCount, 50);
}

INSTANTIATE_TEST_SUITE_P(SquareCart, ForwardPostSmootherTest,
	::testing::Values("mg", "fcgmg", "uamg", "aggregamg"));

// The drawing of the W-cycle in the console used a fixed width of 50 steps, and corrupted the
// heap beyond (7 levels).
//   ./bin/fhhos4 -geo square -mesh cart -mesher inhouse -s mg -cycle W,1,1 -k 0 -n 128 -coarse-size 10
TEST(SolverRegression, WCycleWithManyLevels)
{
	ProgramArguments args = SquareCartArgs(0, 128);
	args.Solver.SolverCode = "mg";
	args.Solver.MG.CycleLetter = 'W';
	args.Solver.MG.WLoops = 2;
	args.Solver.MG.PreSmoothingIterations = 1;
	args.Solver.MG.PostSmoothingIterations = 1;
	args.Solver.MG.MatrixMaxSizeForCoarsestLevel = 10;

	ProgramResults results = RunDiffusionHHO(args, false);
	EXPECT_GT(results.IterationCount, 0);
	EXPECT_LE(results.IterationCount, 20);
}

// A direct solver as preconditioner was cast into an IterativeSolver, which made FCG segfault.
// Being exact, it makes (F)CG converge in 1 iteration.
//   ./bin/fhhos4 -geo square -mesh cart -mesher inhouse -k 1 -n 16 -s {cg|fcg} -preconditioner {ch|lu}
class DirectPreconditionerTest : public ::testing::TestWithParam<std::tuple<std::string, std::string>>
{
};

TEST_P(DirectPreconditionerTest, ConvergesInOneIteration)
{
	auto [solverCode, preconditionerCode] = GetParam();

	ProgramArguments args = SquareCartArgs(1, 16);
	args.Solver.SolverCode = solverCode;
	args.Solver.PreconditionerCode = preconditionerCode;

	ProgramResults results = RunDiffusionHHO(args);
	EXPECT_EQ(results.IterationCount, 1);
}

INSTANTIATE_TEST_SUITE_P(SquareCart, DirectPreconditionerTest, ::testing::Combine(
	::testing::Values("cg", "fcg"),
	::testing::Values("ch", "lu")));

// With the default solver, the prolongations requiring the Galerkin operator (-prolong 4, 5)
// must enable it (with a warning), instead of failing (4) or running without it (5).
//   ./bin/fhhos4 -geo square -mesh cart -mesher inhouse -prolong {4|5}
class GalerkinProlongationDefaultsTest : public ::testing::TestWithParam<int>
{
};

TEST_P(GalerkinProlongationDefaultsTest, DefaultSolverEnablesGalerkinOperator)
{
	ProgramArguments args = SquareCartArgs(0, 16);
	args.Solver.MG.ProlongationCode = GetParam();

	ApplyProgramArgumentDefaults(args);

	EXPECT_EQ(args.Solver.SolverCode, "mg");
	EXPECT_TRUE(args.Solver.MG.UseGalerkinOperator);
	EXPECT_EQ(static_cast<int>(args.Solver.MG.GMG_H_Prolong), GetParam());
}

// With an explicit -s mg, they must be rejected without -g 1.
//   ./bin/fhhos4 -geo square -mesh cart -mesher inhouse -s mg -prolong {4|5}
TEST_P(GalerkinProlongationDefaultsTest, ExplicitMultigridRequiresGalerkinOperator)
{
	ProgramArguments args = SquareCartArgs(0, 16);
	args.Solver.SolverCode = "mg";
	args.Solver.MG.ProlongationCode = GetParam();

	EXPECT_EXIT(ApplyProgramArgumentDefaults(args), ::testing::ExitedWithCode(EXIT_FAILURE), "");
}

INSTANTIATE_TEST_SUITE_P(Prolongation, GalerkinProlongationDefaultsTest, ::testing::Values(4, 5));
#endif // ENABLE_2D

// ThreadLocal: with schedule(static), the thread-local results taken in thread order follow the
// iteration order (the GMSH mesh import and the matrix assemblies rely on the static schedule).
// ThreadLocalCoeffs: each coefficient added in the loop is in the matrix.
class ThreadLocalTest : public ::testing::TestWithParam<std::tuple<BigNumber, int>>
{
protected:
	void TearDown() override { Parallelism::SetNThreads(0); }
};

TEST_P(ThreadLocalTest, ResultsFollowIterationOrder)
{
	auto [loopSize, nThreads] = GetParam();
	Parallelism::SetNThreads(nThreads);

	ThreadLocal<std::vector<BigNumber>> lists;
	#pragma omp parallel for schedule(static)
	for (BigNumber i = 0; i < loopSize; i++)
		lists.Local().push_back(i);

	std::vector<BigNumber> all;
	for (std::vector<BigNumber>& list : lists)
		all.insert(all.end(), list.begin(), list.end());
	ASSERT_EQ(all.size(), loopSize);
	for (BigNumber i = 0; i < loopSize; i++)
		EXPECT_EQ(all[i], i);
}

TEST_P(ThreadLocalTest, CoeffsFillMatrix)
{
	auto [n, nThreads] = GetParam();
	Parallelism::SetNThreads(nThreads);

	ThreadLocalCoeffs coeffs(n, 1);
	#pragma omp parallel for
	for (BigNumber i = 0; i < n; i++)
		coeffs.Local().Add(i, i, i + 1);

	SparseMatrix M(n, n);
	coeffs.Fill(M);
	EXPECT_EQ(M.nonZeros(), n);
	for (BigNumber i = 0; i < n; i++)
		EXPECT_DOUBLE_EQ(M.coeff(i, i), i + 1);
}

INSTANTIATE_TEST_SUITE_P(Sizes, ThreadLocalTest, ::testing::Combine(
	::testing::Values<BigNumber>(0, 3, 1000),
	::testing::Values(1, 4, 16)));

// A parallel loop nested in another one runs on the calling thread only: a ThreadLocal created
// inside the outer loop must have a value for it.
TEST(ThreadLocalNestedTest, CoeffsFillMatrix)
{
	const int nOuter = 32;
	const BigNumber n = 100;
	std::vector<double> traces(nOuter, 0);
	#pragma omp parallel for
	for (int k = 0; k < nOuter; k++)
	{
		ThreadLocalCoeffs coeffs(n, 1);
		#pragma omp parallel for
		for (BigNumber i = 0; i < n; i++)
			coeffs.Local().Add(i, i, k + 1);

		SparseMatrix M(n, n);
		coeffs.Fill(M);
		traces[k] = M.sum();
	}
	for (int k = 0; k < nOuter; k++)
		EXPECT_DOUBLE_EQ(traces[k], n * (k + 1)) << "k=" << k;
}
