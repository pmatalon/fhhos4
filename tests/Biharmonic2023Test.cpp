// Validation tests for "Iterative solution to the biharmonic equation in mixed form discretized
// by the Hybrid High-Order method" (Antonietti, Matalon, Verani, Comput. Math. Appl. 2023),
// reproducibility/2023_biharmonic_problem.md.
//
// The meshes are built by GMSH. The tests check the current iteration counts; the paper's values
// (Tables 1 and 3 of revision 2) are given for comparison, and are all reproduced except one.
// To keep the suite fast, only the smallest mesh sizes and k = 0, 1 are tested.
#include <gtest/gtest.h>
#include <string>
#include <tuple>
#include <vector>
#include "support/RunHelper.h"

using namespace fhhos4_tests;

// Table 1: square, Cartesian mesh. Number of FCG iterations preconditioned by the patch
// preconditioner with neighbourhood depth 8 (-bihar-prec s), or not preconditioned (-bihar-prec no).
class BiharSquareCartTest : public ::testing::TestWithParam<std::tuple<std::string, int, std::vector<ExpectedIterations>>>
{
};

//   ./bin/fhhos4 -pb bihar -geo square -source exp -mesh cart -cs r -s ch -nbh-depth 8 -tol 1e-8 -bihar-prec {s|no} -k {0|1} -n {32|64} -no-cache
TEST_P(BiharSquareCartTest, IterationCounts)
{
	auto [preconditionerCode, k, expected] = GetParam();

	for (const ExpectedIterations& e : expected)
	{
		ProgramArguments args;
		args.Problem.GeoCode = "square";
		args.Problem.SourceCode = "exp";
		args.Discretization.MeshCode = "cart";
		args.Discretization.N = e.N;
		args.Discretization.PolyDegree = k + 1;
		args.Solver.MG.H_CS = H_CoarsStgy::GMSHSplittingRefinement; // -cs r
		args.Solver.SolverCode = "ch";
		args.Solver.NeighbourhoodDepth = 8;
		args.Solver.Tolerance = 1e-8;
		args.Solver.BiHarmonicPreconditionerCode = preconditionerCode;
		args.Actions.UseCache = false; // -no-cache

		ProgramResults results = RunBiHarmonicHHO(args);
		EXPECT_EQ(results.IterationCount, e.Iterations) << "N=" << e.N << " (paper: " << e.PaperIterations << ")";
	}
}

INSTANTIATE_TEST_SUITE_P(Square, BiharSquareCartTest, ::testing::Values(
	std::make_tuple(std::string("s"),  0, std::vector<ExpectedIterations>{ { 32, 13, 13 }, { 64, 19, 19 } }),
	std::make_tuple(std::string("s"),  1, std::vector<ExpectedIterations>{ { 32, 13, 13 }, { 64, 19, 19 } }),
	std::make_tuple(std::string("no"), 0, std::vector<ExpectedIterations>{ { 32, 19, 19 }, { 64, 25, 25 } }),
	// N=32: 30 iterations, as in the paper, when the mesh is loaded from the GMSH cache, whose vertex
	// coordinates differ by ~1e-14 (see tests/README.md): the unpreconditioned FCG is sensitive to it.
	std::make_tuple(std::string("no"), 1, std::vector<ExpectedIterations>{ { 32, 28, 30 }, { 64, 36, 36 } })));

// Table 3: cube, unstructured tetrahedral mesh, k=0, first mesh size (h_0, 3373 elements in the
// paper, 3411 with the current GMSH). Neighbourhood depth 2, and the Laplacian problems solved by
// FCG preconditioned by U-AMG.
class BiharCubeTetTest : public ::testing::TestWithParam<std::tuple<std::string, ExpectedIterations>>
{
};

//   ./bin/fhhos4 -pb bihar -geo cube -source exp -mesh tetra -not-compute-errors -s fcguamg -hp-cs p_h -nbh-depth 2 -bihar-prec-solver bicgstab -tol 1e-8 -k 0 -bihar-prec {s|no} -n 8 -no-cache
TEST_P(BiharCubeTetTest, IterationCounts)
{
	auto [preconditionerCode, e] = GetParam();

	ProgramArguments args;
	args.Problem.GeoCode = "cube";
	args.Problem.SourceCode = "exp";
	args.Discretization.MeshCode = "tetra";
	args.Discretization.N = e.N;
	args.Discretization.PolyDegree = 1; // k = 0
	args.Solver.SolverCode = "fcguamg";
	args.Solver.MG.HP_CS = HP_CoarsStgy::P_then_H; // -hp-cs p_h
	args.Solver.NeighbourhoodDepth = 2;
	args.Solver.BiHarmonicPrecSolverCode = "bicgstab";
	args.Solver.Tolerance = 1e-8;
	args.Solver.BiHarmonicPreconditionerCode = preconditionerCode;
	args.Actions.ComputeErrors = false; // -not-compute-errors
	args.Actions.UseCache = false; // -no-cache

	ProgramResults results = RunBiHarmonicHHO(args);
	EXPECT_EQ(results.IterationCount, e.Iterations) << "N=" << e.N << " (paper: " << e.PaperIterations << ")";
}

INSTANTIATE_TEST_SUITE_P(Cube, BiharCubeTetTest, ::testing::Values(
	std::make_tuple(std::string("s"),  ExpectedIterations{ 8, 14, 14 }),
	std::make_tuple(std::string("no"), ExpectedIterations{ 8, 29, 29 })));
