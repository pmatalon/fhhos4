// Validation tests for "An h-multigrid method for Hybrid High-Order discretizations"
// (Di Pietro, Hulsemann, Matalon, Mycek, Rude, Ruiz, SIAM J. Sci. Comput. 2021),
// reproducibility/2020_MG_for_HHO.md.
//
// The expected iteration counts are the paper's, taken from the CSV files of the accepted
// version (Multigrid_for_HHO_SISC.zip, Results/). In the CSV file names, p = k+1.
// To keep the suite fast, only the smallest mesh sizes and k = 0, 1 are tested.
#include <gtest/gtest.h>
#include <string>
#include <tuple>
#include <vector>
#include "support/RunHelper.h"

using namespace fhhos4_tests;

#ifdef ENABLE_2D
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

	const std::vector<int> MeshSizes = { 32, 64 };
}

// Figure 4.1: V(1,1)-cycle iteration counts on the unit square.
// Reference: 2D_scalability_homogeneous_V11_g0_p{1|2}_{cart|tri}.csv.
class IterationCountTest : public ::testing::TestWithParam<std::tuple<std::string, int, std::vector<int>>>
{
};

//   ./bin/fhhos4 -geo square -mesh {cart|stri} -mesher inhouse -s mg -cycle V,1,1 -k {0|1} -n {32|64}
TEST_P(IterationCountTest, MatchesPaper)
{
	auto [meshCode, k, expectedIterations] = GetParam();

	for (size_t i = 0; i < MeshSizes.size(); i++)
	{
		ProgramResults results = RunDiffusionHHO(SquareVCycleArgs(meshCode, k, MeshSizes[i]), false);
		EXPECT_EQ(results.IterationCount, expectedIterations[i]) << "n=" << MeshSizes[i];
	}
}

INSTANTIATE_TEST_SUITE_P(Square, IterationCountTest, ::testing::Values(
	std::make_tuple(std::string("cart"), 0, std::vector<int>{ 13, 14 }),
	std::make_tuple(std::string("cart"), 1, std::vector<int>{ 16, 18 }),
	std::make_tuple(std::string("stri"), 0, std::vector<int>{ 19, 24 }),
	std::make_tuple(std::string("stri"), 1, std::vector<int>{ 24, 24 })));
#endif // ENABLE_2D

#ifdef ENABLE_3D
// Figure 4.1 documents this configuration as diverging (its CSV has no iteration count): with
// the standard coarsening, the multigrid convergence rate degrades toward 1 at k=0 on
// structured tetrahedral meshes. Observed: 122 iterations at n=16 (vs. 24 for k=1). The run
// still terminates, so check the degradation through the iteration count rather than the exit code.
//   ./bin/fhhos4 -geo cube -mesh stetra -mesher inhouse -s mg -cycle V,2,2 -k 0 -n 16
TEST(HMultigrid2020, DegradesAsDocumented_StructuredTetra_K0)
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

	ProgramResults results = RunDiffusionHHO(args, false);
	EXPECT_GT(results.IterationCount, 60) << "k=0 on stetra is documented as (nearly) diverging";
}
#endif // ENABLE_3D

#ifdef ENABLE_2D
// Figure 4.7: V(1,1)-cycle iteration counts on the heterogeneous Kellogg benchmark.
// Reference: Kellogg_scalability_V11_g0_p{1|2}_cart.csv.
class KelloggTest : public ::testing::TestWithParam<std::tuple<int, std::vector<int>>>
{
};

//   ./bin/fhhos4 -geo square4quadrants -tc kellogg -mesh cart -mesher inhouse -s mg -cycle V,1,1 -k {0|1} -n {32|64}
TEST_P(KelloggTest, MatchesPaper)
{
	auto [k, expectedIterations] = GetParam();

	for (size_t i = 0; i < MeshSizes.size(); i++)
	{
		ProgramArguments args;
		args.Problem.GeoCode = "square4quadrants";
		args.Problem.TestCaseCode = "kellogg";
		args.Discretization.Mesher = "inhouse";
		args.Discretization.MeshCode = "cart";
		args.Discretization.N = MeshSizes[i];
		args.Discretization.PolyDegree = k + 1;
		args.Solver.SolverCode = "mg";
		args.Solver.MG.CycleLetter = 'V';
		args.Solver.MG.PreSmoothingIterations = 1;
		args.Solver.MG.PostSmoothingIterations = 1;

		ProgramResults results = RunDiffusionHHO(args, false);
		EXPECT_EQ(results.IterationCount, expectedIterations[i]) << "n=" << MeshSizes[i];
	}
}

INSTANTIATE_TEST_SUITE_P(SquareFourQuadrants, KelloggTest, ::testing::Values(
	std::make_tuple(0, std::vector<int>{ 13, 13 }),
	std::make_tuple(1, std::vector<int>{ 15, 17 })));

// Figure 4.9(a): with the Galerkin operator and the heterogeneous weighting, the V(0,3)-cycle
// iteration count does not depend on the heterogeneity ratio.
// Reference: 2D_heterogeneity_chiasmus_n64_V03_g1_cart.csv.
class HeterogeneityRatioSweepTest : public ::testing::TestWithParam<std::tuple<int, int>>
{
};

//   ./bin/fhhos4 -geo square4quadrants -mesh cart -mesher inhouse -n 64 -s mg -g 1 -cycle V,0,3 -k {0|1} -heterog {1e0|1e2|1e4|1e6|1e8}
TEST_P(HeterogeneityRatioSweepTest, MatchesPaper)
{
	auto [k, expectedIterations] = GetParam();

	for (double ratio : { 1e0, 1e2, 1e4, 1e6, 1e8 })
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

		ProgramResults results = RunDiffusionHHO(args, false);
		EXPECT_EQ(results.IterationCount, expectedIterations) << "ratio=" << ratio;
	}
}

INSTANTIATE_TEST_SUITE_P(SquareFourQuadrants, HeterogeneityRatioSweepTest, ::testing::Values(
	std::make_tuple(0, 7),
	std::make_tuple(1, 9)));
#endif // ENABLE_2D
