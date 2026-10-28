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

//------------------------------------------- GMSH's cube -------------------------------------------//

namespace
{
	ProgramArguments CubeArgs(const std::string& meshCode, const std::string& mesher, int k, const std::string& solverCode)
	{
		ProgramArguments args;
		args.Problem.GeoCode = "cube";
		args.Discretization.Mesher = mesher;
		args.Discretization.MeshCode = meshCode;
		args.Discretization.N = 8;
		args.Discretization.PolyDegree = k + 1;
		args.Solver.SolverCode = solverCode;
		args.Actions.UseCache = false;
		return args;
	}
}

// GMSH's Cartesian cube (hexahedra) is the in-house one: the same L2 error. It failed since 2022: its quadrilateral
// faces were refused, the corners of its hexahedra and faces taken in an order of GMSH's nodes that no longer held, the
// sides of the faces in their elements (their outer normals) not set, and its coordinates are exact to ~1e-12 only.
//   ./bin/fhhos4 -geo cube -mesh cart -mesher {gmsh|inhouse} -s ch -k {0|1|2} -n 8 -no-cache
class GMSHCartesianCubeTest : public ::testing::TestWithParam<int>
{
};

TEST_P(GMSHCartesianCubeTest, SameAsInHouse)
{
	int k = GetParam();
	double gmshError = RunDiffusionHHO(CubeArgs("cart", "gmsh", k, "ch")).L2Error;
	double inHouseError = RunDiffusionHHO(CubeArgs("cart", "inhouse", k, "ch")).L2Error;
	ASSERT_GT(inHouseError, 0);
	EXPECT_NEAR(gmshError, inHouseError, 1e-6 * inHouseError);
}

INSTANTIATE_TEST_SUITE_P(Cube, GMSHCartesianCubeTest, ::testing::Values(0, 1, 2));

// The mixed conditions of the cube (-bc m): Dirichlet on its bottom face (z = 0, "bottomBoundary"), Neumann on the
// others. The test case had no boundary functions (std::bad_function_call), and the faces of GMSH's cube no name: all
// Neumann, a singular system (200 iterations). No exact solution: the iterations only.
//   ./bin/fhhos4 -geo cube -mesh {tetra|cart} -mesher gmsh -bc m -s fcglibuamg -k 1 -n 8 -no-cache
class CubeMixedConditionsTest : public ::testing::TestWithParam<std::string>
{
};

TEST_P(CubeMixedConditionsTest, Converges)
{
	ProgramArguments args = CubeArgs(GetParam(), "gmsh", 1, "fcglibuamg");
	args.Problem.BCCode = "m";
	EXPECT_LT(RunDiffusionHHO(args).IterationCount, 40);
}

INSTANTIATE_TEST_SUITE_P(Cube, CubeMixedConditionsTest, ::testing::Values(std::string("tetra"), std::string("cart")));

// The in-house cube has no boundary parts: all its faces would be Dirichlet. Refused (Mesh::SetBoundaryConditions()).
//   ./bin/fhhos4 -geo cube -mesh cart -mesher inhouse -bc m -k 1 -n 8
TEST(CubeMixedConditions, RefusedWithoutBoundaryParts)
{
	ProgramArguments args = CubeArgs("cart", "inhouse", 1, "ch");
	args.Problem.BCCode = "m";
	EXPECT_THROW(RunDiffusionHHO(args), fhhos4::Error);
}
#endif // ENABLE_3D
