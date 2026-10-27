// The full Neumann problem (-tc fullneumann): a singular system, solved up to a constant (the program projects the
// right-hand side onto the image of the matrix and translates the solution to the average of the exact one).
// It crashed (segmentation fault) since 2022-07 (42f1004e): the compatibility condition was computed by the discrete
// spaces of the problem before the assembly that builds them. On the Cartesian meshes, the condition was also wrong
// (-2 instead of 0): the integrals of a function on the faces of a CartesianShape ignored their orientation.
#include <gtest/gtest.h>
#include <string>
#include <tuple>
#include <vector>
#include "support/RunHelper.h"
#include "Geometry/CartesianShape.h"

using namespace fhhos4_tests;

// The integral of a function on the faces of the Cartesian elements (src/Geometry/CartesianShape.h): along their
// orientation, at their position. f is linear: the quadrature is exact.
class CartesianFaceIntegral : public ::testing::Test
{
protected:
	void SetUp() override { GaussLegendre::Init(); }
	void TearDown() override { GaussLegendre::Free(); }
};

TEST_F(CartesianFaceIntegral, FollowsTheOrientation)
{
	DomFunction f2D = [](const DomPoint& p) { return p.X + 2 * p.Y; };
	DomPoint o1(1, 0.25);
	CartesianShape<2, 1> vertical(&o1, 0.25, CartesianShapeOrientation::Vertical); // x = 1, y in [0.25, 0.5]
	EXPECT_NEAR(vertical.Integral(f2D), 0.4375, 1e-14);
	DomPoint o2(0.5, 1);
	CartesianShape<2, 1> horizontal(&o2, 0.25, CartesianShapeOrientation::Horizontal); // y = 1, x in [0.5, 0.75]
	EXPECT_NEAR(horizontal.Integral(f2D), 0.65625, 1e-14);

	DomFunction f3D = [](const DomPoint& p) { return p.X + 2 * p.Y + 3 * p.Z; };
	DomPoint o3(0, 0, 1);
	CartesianShape<3, 2> xoy(&o3, 0.5, CartesianShapeOrientation::InXOY); // z = 1
	EXPECT_NEAR(xoy.Integral(f3D), 0.9375, 1e-14);
	DomPoint o4(0, 1, 0);
	CartesianShape<3, 2> xoz(&o4, 0.5, CartesianShapeOrientation::InXOZ); // y = 1
	EXPECT_NEAR(xoz.Integral(f3D), 0.75, 1e-14);
	DomPoint o5(1, 0, 0);
	CartesianShape<3, 2> yoz(&o5, 0.5, CartesianShapeOrientation::InYOZ); // x = 1
	EXPECT_NEAR(yoz.Integral(f3D), 0.5625, 1e-14);
}

#ifdef ENABLE_2D
namespace
{
	// mesher "" : the default (GMSH for tri)
	ProgramArguments FullNeumannArgs(const std::string& meshCode, const std::string& mesher, int k, int n, const std::string& solverCode)
	{
		ProgramArguments args;
		args.Problem.GeoCode = "square";
		args.Problem.TestCaseCode = "fullneumann";
		if (!mesher.empty())
			args.Discretization.Mesher = mesher;
		args.Discretization.MeshCode = meshCode;
		args.Discretization.N = n;
		args.Discretization.PolyDegree = k + 1;
		args.Solver.SolverCode = solverCode;
		args.Actions.UseCache = false;
		return args;
	}
}

// Mesh, mesher, k, L2 error at n = 32
class FullNeumannTest : public ::testing::TestWithParam<std::tuple<std::string, std::string, int, double>>
{
};

// The L2 error converges at the theoretical rate (h^2 for k=0, h^(k+2) for k>=1; direct solver), and the solvers give
// the same solution at n = 64, where the multigrid methods have several levels (tolerance 1e-11: at k = 2, the L2 error
// at the default tolerance is that of the solver): the geometric multigrid (which applies the zero-mean and
// compatibility conditions on its levels), U-AMG in the program and through the library (on the condensed system),
// with their kernel handling (SingularMatrix: without it, the singular coarsest matrix was factorized with a tiny pivot,
// and tri k=2 took 198 iterations through the library).
//   ./bin/fhhos4 -geo square -tc fullneumann -mesh {tri|cart -mesher inhouse} -k {0|1|2} -n {8|16|32|64} -s ch -no-cache
//   ./bin/fhhos4 -geo square -tc fullneumann -mesh {tri|cart -mesher inhouse} -k {0|1|2} -n 64 -s {mg|fcguamg|fcglibuamg} -tol 1e-11 -no-cache
TEST_P(FullNeumannTest, ConvergesAndSolversAgree)
{
	auto [meshCode, mesher, k, expectedError] = GetParam();

	std::vector<double> h;
	std::vector<double> errors;
	for (int n : { 8, 16, 32, 64 })
	{
		ProgramResults results = RunDiffusionHHO(FullNeumannArgs(meshCode, mesher, k, n, "ch"));
		ASSERT_GT(results.L2Error, 0);
		h.push_back(1.0 / n);
		errors.push_back(results.L2Error);
	}
	double order = EstimateConvergenceOrder(h, errors);
	double expectedOrder = (k == 0) ? 2.0 : static_cast<double>(k + 2);
	EXPECT_NEAR(order, expectedOrder, 0.3);
	EXPECT_NEAR(errors[2], expectedError, 0.01 * expectedError); // n = 32

	double directError = errors.back();
	for (std::string solverCode : { "mg", "fcguamg", "fcglibuamg" })
	{
		ProgramArguments args = FullNeumannArgs(meshCode, mesher, k, 64, solverCode);
		args.Solver.Tolerance = 1e-11;
		ProgramResults results = RunDiffusionHHO(args);
		EXPECT_GT(results.IterationCount, 1) << solverCode; // several levels
		EXPECT_LT(results.IterationCount, 40) << solverCode;
		EXPECT_NEAR(results.L2Error, directError, 1e-4 * directError) << solverCode;
	}
}

INSTANTIATE_TEST_SUITE_P(Square, FullNeumannTest,
	::testing::Values(
		std::make_tuple(std::string("tri"), std::string(""), 0, 5.67e-4),
		std::make_tuple(std::string("tri"), std::string(""), 1, 1.65e-6),
		std::make_tuple(std::string("tri"), std::string(""), 2, 5.35e-9),
		std::make_tuple(std::string("cart"), std::string("inhouse"), 0, 5.27e-4),
		std::make_tuple(std::string("cart"), std::string("inhouse"), 1, 3.26e-6)),
	[](const ::testing::TestParamInfo<FullNeumannTest::ParamType>& info)
	{
		return std::get<0>(info.param) + "_k" + std::to_string(std::get<2>(info.param));
	});
#endif // ENABLE_2D
