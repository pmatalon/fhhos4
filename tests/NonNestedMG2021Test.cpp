// Validation tests for "Towards robust, fast solutions of elliptic equations on complex domains
// through hybrid high-order discretizations and non-nested multigrid methods" (Di Pietro,
// Hulsemann, Matalon, Mycek, Rude, Ruiz, Int. J. Numer. Methods Eng. 2021),
// reproducibility/2021_non_nested_MG_for_HHO.md.
//
// The meshes are built by GMSH, whose current version produces slightly different meshes from
// the paper's, so the iteration counts may differ a little from the paper's. The tests check
// the current counts; the paper's values, from the CSV files of the accepted version
// (additional_files.zip, Results/; in the file names, p = k+1), are given for comparison.
// To keep the suite fast, only the smallest mesh sizes and k = 0, 1 are tested.
#include <gtest/gtest.h>
#include <string>
#include <tuple>
#include <vector>
#include "support/RunHelper.h"

using namespace fhhos4_tests;

namespace
{
	enum class Coarsening { Remeshing, Agglomeration };

	ProgramArguments SquareArgs(Coarsening coarsening, int prolongationCode, bool barycentricSubtriangulation, int k, int n)
	{
		ProgramArguments args;
		args.Problem.GeoCode = "square";
		args.Discretization.MeshCode = "tri";
		args.Discretization.N = n;
		args.Discretization.PolyDegree = k + 1;
		args.Solver.SolverCode = "mg";
		args.Solver.MG.H_CS = coarsening == Coarsening::Remeshing ? H_CoarsStgy::IndependentRemeshing : H_CoarsStgy::AgglomerationCoarseningByFaceNeighbours;
		args.Solver.MG.ProlongationCode = prolongationCode;
		if (barycentricSubtriangulation)
			args.Solver.MG.SubtriangulationMethodForApproxL2Proj = PolygonalTriangulation::Barycentric;
		args.Actions.UseCache = false; // -no-cache
		return args;
	}

	void CheckIterations(const ProgramArguments& argsWithoutN, const std::vector<ExpectedIterations>& expected)
	{
		for (const ExpectedIterations& e : expected)
		{
			ProgramArguments args = argsWithoutN;
			args.Discretization.N = e.N;
			ProgramResults results = RunDiffusionHHO(args);
			EXPECT_EQ(results.IterationCount, e.Iterations) << "N=" << e.N << " (paper: " << e.PaperIterations << ")";
		}
	}
}

#ifdef ENABLE_2D
//----------------------------------------------------------------------------------//
// Figure 6: square, independent remeshing of the coarse levels (-cs m), V(0,3)     //
//----------------------------------------------------------------------------------//

class NonNestedSquareRemeshingTest : public ::testing::TestWithParam<std::tuple<int, int, std::vector<ExpectedIterations>>>
{
};

// (a) -prolong 7: exact L2-projection; (b) -prolong 9: approximate L2-projection by subtriangulation.
// Reference: square_V03_g0_csm_prolong{7|9}_p{1|2}.csv.
//   ./bin/fhhos4 -geo square -mesh tri -s mg -cs m -prolong {7|9} -k {0|1} -n {32|64} -no-cache
TEST_P(NonNestedSquareRemeshingTest, IterationCounts)
{
	auto [prolongationCode, k, expected] = GetParam();
	CheckIterations(SquareArgs(Coarsening::Remeshing, prolongationCode, false, k, 0), expected);
}

INSTANTIATE_TEST_SUITE_P(Square, NonNestedSquareRemeshingTest, ::testing::Values(
	std::make_tuple(7, 0, std::vector<ExpectedIterations>{ { 32,  8,  7 }, { 64,  8,  8 } }),
	std::make_tuple(7, 1, std::vector<ExpectedIterations>{ { 32, 12, 12 }, { 64, 13, 13 } }),
	std::make_tuple(9, 0, std::vector<ExpectedIterations>{ { 32,  8,  8 }, { 64,  9,  8 } }),
	std::make_tuple(9, 1, std::vector<ExpectedIterations>{ { 32, 12, 12 }, { 64, 13, 12 } })));

//----------------------------------------------------------------------------------//
// Figure 9: square, agglomeration coarsening (-cs n), V(0,3)                       //
//----------------------------------------------------------------------------------//

// The agglomeration depends on the thread scheduling: sequential execution (-threads 1).
class NonNestedSquareAgglomerationTest : public ::testing::TestWithParam<std::tuple<bool, int, std::vector<ExpectedIterations>>>
{
};

// (a) approximate L2-projection with the optimal subtriangulation (default), which gives the
// same results as the exact L2-projection; (b) -subtri-meth bary: barycentric subtriangulation.
// Reference: square_V03_g0_csn_prolong7_p{1|2}.csv (a), square_V03_g0_csn_prolong9barycentric_p{1|2}.csv (b).
//   ./bin/fhhos4 -geo square -mesh tri -s mg -cs n -prolong 9 [-subtri-meth bary] -k {0|1} -n {32|64} -no-cache -threads 1
TEST_P(NonNestedSquareAgglomerationTest, IterationCounts)
{
	auto [barycentric, k, expected] = GetParam();
	SequentialExecution sequential;
	CheckIterations(SquareArgs(Coarsening::Agglomeration, 9, barycentric, k, 0), expected);
}

INSTANTIATE_TEST_SUITE_P(Square, NonNestedSquareAgglomerationTest, ::testing::Values(
	std::make_tuple(false, 0, std::vector<ExpectedIterations>{ { 32,  9,  8 }, { 64, 10, 10 } }),
	std::make_tuple(false, 1, std::vector<ExpectedIterations>{ { 32, 13, 12 }, { 64, 13, 13 } }),
	std::make_tuple(true,  0, std::vector<ExpectedIterations>{ { 32,  9,  8 }, { 64, 10,  9 } }),
	std::make_tuple(true,  1, std::vector<ExpectedIterations>{ { 32, 13, 12 }, { 64, 13, 14 } })));
#endif // ENABLE_2D

#ifdef ENABLE_3D
//----------------------------------------------------------------------------------//
// Figure 8: cube, independent remeshing (-cs m), approx. L2-projection, V(0,6)     //
//----------------------------------------------------------------------------------//

class NonNestedCubeRemeshingTest : public ::testing::TestWithParam<std::tuple<int, std::vector<ExpectedIterations>>>
{
};

// Reference: cube_V06_g0_csm_prolong9_p{1|2}.csv.
//   ./bin/fhhos4 -geo cube -mesh tetra -s mg -cs m -prolong 9 -k {0|1} -n 8 -no-cache
TEST_P(NonNestedCubeRemeshingTest, IterationCounts)
{
	auto [k, expected] = GetParam();
	ProgramArguments args = SquareArgs(Coarsening::Remeshing, 9, false, k, 0);
	args.Problem.GeoCode = "cube";
	args.Discretization.MeshCode = "tetra";
	CheckIterations(args, expected);
}

INSTANTIATE_TEST_SUITE_P(Cube, NonNestedCubeRemeshingTest, ::testing::Values(
	std::make_tuple(0, std::vector<ExpectedIterations>{ { 8, 12, 11 } }),
	std::make_tuple(1, std::vector<ExpectedIterations>{ { 8, 16, 13 } })));
#endif // ENABLE_3D
