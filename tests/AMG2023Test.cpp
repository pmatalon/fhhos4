// Validation tests for "Algebraic multigrid preconditioner for statically condensed systems
// arising from lowest-order hybrid discretizations" (Di Pietro, Hulsemann, Matalon, Mycek, Rude,
// SIAM J. Sci. Comput. 2023), reproducibility/2023_AMG_for_hybrid_methods.md.
//
// The meshes are built by GMSH, whose current version produces slightly different meshes from
// the paper's, so the iteration counts may differ a little from the paper's. The tests check
// the current counts; the paper's values, from the CSV files of the accepted version (Results/),
// are given for comparison. To keep the suite fast, only the smallest mesh sizes are tested.
// AGMG, the third solver of the comparison, is an external library that is not built by default.
#include <gtest/gtest.h>
#include <string>
#include <tuple>
#include <vector>
#include "support/RunHelper.h"

#ifdef ENABLE_3D
using namespace fhhos4_tests;

// Figure 4.4 (Cube-tet, k=0): FCG preconditioned by U-AMG (-s fcguamg) or C-AMG (-s fcgaggregamg).
// Reference: cube_fcgcamg.csv (U-AMG), cube_fcgaggregamg.csv (C-AMG).
class AMGCubeTetTest : public ::testing::TestWithParam<std::tuple<std::string, std::vector<ExpectedIterations>>>
{
};

//   ./bin/fhhos4 -geo cube -mesh tetra -k 0 -n {16|32} -s {fcguamg|fcgaggregamg} -no-cache
TEST_P(AMGCubeTetTest, IterationCounts)
{
	auto [solverCode, expected] = GetParam();

	for (const ExpectedIterations& e : expected)
	{
		ProgramArguments args;
		args.Problem.GeoCode = "cube";
		args.Discretization.MeshCode = "tetra";
		args.Discretization.N = e.N;
		args.Discretization.PolyDegree = 1; // k = 0
		args.Solver.SolverCode = solverCode;
		args.Actions.UseCache = false; // -no-cache

		ProgramResults results = RunDiffusionHHO(args);
		EXPECT_EQ(results.IterationCount, e.Iterations) << "N=" << e.N << " (paper: " << e.PaperIterations << ")";
	}
}

INSTANTIATE_TEST_SUITE_P(Cube, AMGCubeTetTest, ::testing::Values(
	std::make_tuple(std::string("fcguamg"),      std::vector<ExpectedIterations>{ { 16, 26, 26 }, { 32, 28, 28 } }),
	std::make_tuple(std::string("fcgaggregamg"), std::vector<ExpectedIterations>{ { 16, 23, 22 }, { 32, 25, 25 } })));
#endif // ENABLE_3D
