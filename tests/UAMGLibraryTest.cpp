// Tests of the library fhhos4 (library/), U-AMG for the codes that discretize the problem themselves.
// - UAMGLibraryProgram: the program runs U-AMG through the library (-s libuamg, fcglibuamg, gcrlibuamg), as an external
//   code calls it, and must get the iteration counts and solutions of the U-AMG compiled in the program (uamg, fcguamg,
//   gcruamg), including with L2-orthonormal bases, whose constant does not have coordinate 1. Through both interfaces of
//   the library: the uncondensed system (default: the library condenses it and computes the cell unknowns) and the
//   condensed one (-libuamg-condensed: the program condenses it).
// - UAMGLibrary: the public API alone (<fhhos4>), on a hybrid system built here.
#include <gtest/gtest.h>
#include <cmath>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "fhhos4"
#include "support/RunHelper.h"

using namespace fhhos4_tests;

//------------------------------------------------------------------------------------------------------------------//
//                                         Through the program (-s libuamg)                                           //
//------------------------------------------------------------------------------------------------------------------//

struct LibraryCase
{
	std::string Name;
	std::string Geo;
	std::string Mesh;
	std::string Mesher; // "" for the default
	int N;
	int K;
	int BasesOrthogonalization; // -e-ogb and -f-ogb; -1 for the defaults
};

inline std::ostream& operator<<(std::ostream& os, const LibraryCase& c)
{
	return os << c.Name;
}

// The program on the case with the solver solverCode. nonsymmetric: -uamg-nonsym; libraryCondensed: -libuamg-condensed.
ProgramResults RunCase(const LibraryCase& c, const std::string& solverCode, bool nonsymmetric = false, bool libraryCondensed = false)
{
	ProgramArguments args;
	args.Problem.GeoCode = c.Geo;
	args.Discretization.MeshCode = c.Mesh;
	if (!c.Mesher.empty())
		args.Discretization.Mesher = c.Mesher;
	args.Discretization.N = c.N;
	args.Discretization.PolyDegree = c.K + 1;
	if (c.BasesOrthogonalization != -1)
	{
		args.Discretization.OrthogonalizeElemBasesCode = c.BasesOrthogonalization;
		args.Discretization.OrthogonalizeFaceBasesCode = c.BasesOrthogonalization;
	}
	args.Solver.SolverCode = solverCode;
	args.Solver.MG.UAMGNonsymmetric = nonsymmetric;
	args.Solver.MG.UAMGLibraryCondensed = libraryCondensed;
	args.Actions.UseCache = false;
	return RunDiffusionHHO(args);
}

class UAMGLibraryProgramTest : public ::testing::TestWithParam<std::tuple<LibraryCase, std::string>>
{
};

// The library code (libuamg, fcglibuamg, gcrlibuamg) against the program's (uamg, fcguamg, gcruamg), through both
// interfaces of the library
TEST_P(UAMGLibraryProgramTest, SameAsProgram)
{
	auto [c, libraryCode] = GetParam();
	std::string programCode = libraryCode;
	programCode.replace(programCode.find("libuamg"), 7, "uamg");

	ProgramResults program = RunCase(c, programCode);
	EXPECT_GT(program.IterationCount, 0);
	EXPECT_LT(program.IterationCount, 200); // converged (default maximum: 200)
	for (bool libraryCondensed : { false, true })
	{
		ProgramResults library = RunCase(c, libraryCode, false, libraryCondensed);
		std::string label = libraryCode + (libraryCondensed ? " -libuamg-condensed" : "") + " vs " + programCode;
		EXPECT_EQ(library.IterationCount, program.IterationCount) << label;
		EXPECT_NEAR(library.L2Error, program.L2Error, 1e-6 * program.L2Error) << label;
	}
}

#ifdef ENABLE_2D
// U-AMG is a non-linear preconditioner: the library gives it with its FCG only, never to another Krylov solver
TEST(UAMGLibraryProgram, NotAPreconditionerOfTheConjugateGradient)
{
	ProgramArguments args;
	args.Problem.GeoCode = "square";
	args.Discretization.MeshCode = "tri";
	args.Discretization.N = 8;
	args.Discretization.PolyDegree = 1;
	args.Solver.SolverCode = "cglibuamg";
	args.Actions.UseCache = false;
	EXPECT_THROW(RunDiffusionHHO(args), fhhos4::Error);
}
#endif

std::vector<LibraryCase> LibraryCases()
{
	std::vector<LibraryCase> cases;
#ifdef ENABLE_2D
	cases.push_back({ "Square_tri_n16_k0", "square", "tri", "", 16, 0, -1 });
	cases.push_back({ "Square_tri_n16_k1", "square", "tri", "", 16, 1, -1 });
	cases.push_back({ "Square_tri_n16_k1_orthonormal", "square", "tri", "", 16, 1, 3 });
#endif
#ifdef ENABLE_3D
	cases.push_back({ "Cube_tetra_n8_k0", "cube", "tetra", "", 8, 0, -1 });
	cases.push_back({ "Cube_tetra_n8_k0_orthonormal", "cube", "tetra", "", 8, 0, 3 });
	cases.push_back({ "Cube_tetra_n8_k1_orthonormal", "cube", "tetra", "", 8, 1, 3 });
	cases.push_back({ "Cube_cart_n8_k2", "cube", "cart", "inhouse", 8, 2, -1 });
#endif
	return cases;
}

INSTANTIATE_TEST_SUITE_P(UAMGLibraryProgram, UAMGLibraryProgramTest,
	::testing::Combine(::testing::ValuesIn(LibraryCases()), ::testing::Values("fcglibuamg", "gcrlibuamg", "libuamg")),
	[](const ::testing::TestParamInfo<UAMGLibraryProgramTest::ParamType>& info)
	{
		return std::get<0>(info.param).Name + "_" + std::get<1>(info.param);
	});

// U-AMG for a non-symmetric matrix (-uamg-nonsym: A_FT = A_TF^T given, all the blocks in full) on the symmetric matrix:
// it coarsens its symmetric part (itself) and computes the Galerkin products of the matrix, equal to the symmetric
// ones up to rounding: the iteration count and solution of the symmetric U-AMG, with GCR. In the program (gcruamg) and
// through the library (gcrlibuamg: SetupNonsymmetric(), and SetupNonsymmetricCondensed() with -libuamg-condensed).
class UAMGNonsymmetricProgramTest : public ::testing::TestWithParam<LibraryCase>
{
};

TEST_P(UAMGNonsymmetricProgramTest, SameAsSymmetric)
{
	LibraryCase c = GetParam();
	ProgramResults symmetric = RunCase(c, "gcruamg");
	EXPECT_GT(symmetric.IterationCount, 0);
	EXPECT_LT(symmetric.IterationCount, 200);
	for (auto [code, libraryCondensed] : { std::pair{ "gcruamg", false }, { "gcrlibuamg", false }, { "gcrlibuamg", true } })
	{
		ProgramResults nonsymmetric = RunCase(c, code, true, libraryCondensed);
		std::string label = std::string(code) + " -uamg-nonsym" + (libraryCondensed ? " -libuamg-condensed" : "");
		EXPECT_EQ(nonsymmetric.IterationCount, symmetric.IterationCount) << label;
		EXPECT_NEAR(nonsymmetric.L2Error, symmetric.L2Error, 1e-6 * symmetric.L2Error) << label;
	}
}

INSTANTIATE_TEST_SUITE_P(UAMGNonsymmetricProgram, UAMGNonsymmetricProgramTest, ::testing::ValuesIn(LibraryCases()),
	[](const ::testing::TestParamInfo<LibraryCase>& info) { return info.param.Name; });

// -uamg-nonsym refuses FCG (symmetric matrices only)
#ifdef ENABLE_2D
TEST(UAMGNonsymmetricProgram, NotWithFCG)
{
	ProgramArguments args;
	args.Problem.GeoCode = "square";
	args.Discretization.MeshCode = "tri";
	args.Discretization.N = 8;
	args.Discretization.PolyDegree = 1;
	args.Solver.SolverCode = "fcguamg";
	args.Solver.MG.UAMGNonsymmetric = true;
	EXPECT_EXIT(ApplyProgramArgumentDefaults(args), ::testing::ExitedWithCode(EXIT_FAILURE), "");
}

// -libuamg-condensed requires the library
TEST(UAMGLibraryProgram, CondensedOptionRequiresTheLibrary)
{
	ProgramArguments args;
	args.Problem.GeoCode = "square";
	args.Discretization.MeshCode = "tri";
	args.Discretization.N = 8;
	args.Discretization.PolyDegree = 1;
	args.Solver.SolverCode = "fcguamg";
	args.Solver.MG.UAMGLibraryCondensed = true;
	EXPECT_EXIT(ApplyProgramArgumentDefaults(args), ::testing::ExitedWithCode(EXIT_FAILURE), "");
}
#endif

//------------------------------------------------------------------------------------------------------------------//
//                                     Through the public API only (<fhhos4>)                                    //
//------------------------------------------------------------------------------------------------------------------//

// Hybrid system of degree 0 on the n x n Cartesian mesh of the unit square: in each cell T, of coefficient a_T, the
// local matrix couples the cell DoF to each of its 4 faces with the weight a_T (the constants are in its kernel).
// The DoFs of the boundary faces (Dirichlet) are eliminated. Column-major matrices with int indices, as Eigen's
// default (and HArDCore3D's).
// Non-symmetric with asymmetry > 0, like an advection in the direction (1, 1): the cell-face coupling is
// -a_T (1 + asymmetry s) in A_TF and -a_T (1 - asymmetry s) in A_FT, s = 1 on the right and top faces of the cell (the
// outflow), -1 on the others. The symmetric part is the symmetric system: positive definite.
struct HybridSystem
{
	using Matrix = Eigen::SparseMatrix<double>;
	Matrix A_TT, A_TF, A_FT, A_FF, A; // A: the condensed matrix
	Eigen::VectorXd CellInterpOfOne, FaceInterpOfOne, b;
	bool Symmetric;

	// cellScaling, faceScaling: the cell (resp. face) basis function is multiplied by these factors
	HybridSystem(int n, const Eigen::VectorXd& cellScaling = Eigen::VectorXd(), const Eigen::VectorXd& faceScaling = Eigen::VectorXd(), double asymmetry = 0)
	{
		Symmetric = asymmetry == 0;
		int nCells = n * n;
		// Interior faces: vertical ones (between the cells (i, j) and (i+1, j)), then horizontal ones
		int nVertical = (n - 1) * n;
		int nFaces = 2 * nVertical;
		auto vertical = [n](int i, int j) { return j * (n - 1) + i; };              // between (i, j) and (i+1, j)
		auto horizontal = [n, nVertical](int i, int j) { return nVertical + j * n + i; }; // between (i, j) and (i, j+1)

		std::mt19937 generator(42);
		std::uniform_real_distribution<double> coefficient(1, 10);

		Eigen::VectorXd sT = cellScaling.size() > 0 ? cellScaling : Eigen::VectorXd::Ones(nCells);
		Eigen::VectorXd sF = faceScaling.size() > 0 ? faceScaling : Eigen::VectorXd::Ones(nFaces);

		std::vector<Eigen::Triplet<double>> TT, TF, FT, FF;
		for (int j = 0; j < n; j++)
		{
			for (int i = 0; i < n; i++)
			{
				int T = j * n + i;
				double a = coefficient(generator);
				TT.emplace_back(T, T, 4 * a * sT[T] * sT[T]);
				std::vector<std::pair<int, int>> faces; // face, s (1: outflow)
				if (i > 0)     faces.push_back({ vertical(i - 1, j), -1 });
				if (i < n - 1) faces.push_back({ vertical(i, j), 1 });
				if (j > 0)     faces.push_back({ horizontal(i, j - 1), -1 });
				if (j < n - 1) faces.push_back({ horizontal(i, j), 1 });
				for (auto [F, s] : faces)
				{
					TF.emplace_back(T, F, -a * (1 + asymmetry * s) * sT[T] * sF[F]);
					FT.emplace_back(F, T, -a * (1 - asymmetry * s) * sF[F] * sT[T]);
					FF.emplace_back(F, F, a * sF[F] * sF[F]);
				}
			}
		}
		A_TT.resize(nCells, nCells);
		A_TT.setFromTriplets(TT.begin(), TT.end());
		A_TF.resize(nCells, nFaces);
		A_TF.setFromTriplets(TF.begin(), TF.end());
		A_FT.resize(nFaces, nCells);
		A_FT.setFromTriplets(FT.begin(), FT.end());
		A_FF.resize(nFaces, nFaces);
		A_FF.setFromTriplets(FF.begin(), FF.end());

		Eigen::VectorXd invDiag = A_TT.diagonal().cwiseInverse();
		A = A_FF - Matrix(A_FT * invDiag.asDiagonal() * A_TF);
		A.makeCompressed();

		CellInterpOfOne = sT.cwiseInverse();
		FaceInterpOfOne = sF.cwiseInverse();
		// Right-hand side of the condensed system for a unit source in the cells: b = -A_FT A_TT^{-1} f
		Eigen::VectorXd f = sT;
		b = -A_FT * invDiag.cwiseProduct(f);
	}
};

// The solver for HybridSystem: degree 0 in 2D. With the sequential smoothers: the default on several threads, the
// hybrid one, depends on the number of threads (tested apart).
fhhos4::Solver ToySolver()
{
	fhhos4::Solver solver;
	solver.Dimension = 2;
	solver.FaceDegree = 0;
	solver.CellDegree = 0;
	solver.Tolerance = 1e-10;
	solver.PreSmoother = "bgs";
	solver.PostSmoother = "rbgs";
	return solver;
}

// Relative error of x against the direct solution
double RelativeError(const HybridSystem& system, const Eigen::VectorXd& x)
{
	Eigen::VectorXd exact;
	if (system.Symmetric)
		exact = Eigen::SimplicialLDLT<HybridSystem::Matrix>(system.A).solve(system.b);
	else
	{
		Eigen::SparseLU<HybridSystem::Matrix> direct(system.A);
		exact = direct.solve(system.b);
	}
	return (x - exact).norm() / exact.norm();
}

HybridSystem NonsymmetricSystem(int n, double asymmetry)
{
	return HybridSystem(n, Eigen::VectorXd(), Eigen::VectorXd(), asymmetry);
}

TEST(UAMGLibrary, SolvesTheSystem)
{
	HybridSystem system(64);
	fhhos4::Solver uamg = ToySolver();
	uamg.SetupSymmetricCondensed(system.A, system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
	EXPECT_TRUE(uamg.IsSetUp());
	EXPECT_GT(uamg.NumberOfLevels(), 2);

	Eigen::VectorXd x = Eigen::VectorXd::Zero(system.b.size());
	fhhos4::Result result = uamg.Solve(system.b, x);
	EXPECT_TRUE(result.Converged);
	EXPECT_LT(result.RelativeResidual, 1e-10);
	EXPECT_GT(result.Iterations, 0);
	EXPECT_LT(result.Iterations, 30);
	EXPECT_LT(RelativeError(system, x), 1e-8);

	// Non-zero initial guess: x is the solution already
	fhhos4::Result again = uamg.Solve(system.b, x);
	EXPECT_LE(again.Iterations, 1);
}

// Without A, the library computes it from the blocks: same solver
TEST(UAMGLibrary, CondensedMatrixComputedFromTheBlocks)
{
	HybridSystem system(64);
	fhhos4::Solver withA = ToySolver(), withoutA = ToySolver();
	withA.SetupSymmetricCondensed(system.A, system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
	withoutA.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);

	Eigen::VectorXd x1 = Eigen::VectorXd::Zero(system.b.size()), x2 = x1;
	fhhos4::Result r1 = withA.Solve(system.b, x1);
	fhhos4::Result r2 = withoutA.Solve(system.b, x2);
	EXPECT_EQ(withoutA.NumberOfLevels(), withA.NumberOfLevels());
	EXPECT_EQ(r2.Iterations, r1.Iterations);
	EXPECT_LT((x2 - x1).norm() / x1.norm(), 1e-9);

	// Lower triangular parts only of the symmetric blocks
	HybridSystem::Matrix A_TT = system.A_TT.triangularView<Eigen::Lower>(), A_FF = system.A_FF.triangularView<Eigen::Lower>();
	fhhos4::Solver lower = ToySolver();
	lower.SetupSymmetric(A_TT, system.A_TF, A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
	Eigen::VectorXd x3 = Eigen::VectorXd::Zero(system.b.size());
	EXPECT_EQ(lower.Solve(system.b, x3).Iterations, r1.Iterations);
	EXPECT_LT((x3 - x1).norm() / x1.norm(), 1e-9);
}

// The matrix of the uncondensed system [A_TT A_TF; A_FT A_FF]
HybridSystem::Matrix UncondensedMatrix(const HybridSystem& system)
{
	std::vector<Eigen::Triplet<double>> triplets;
	auto add = [&](const HybridSystem::Matrix& M, int firstRow, int firstCol)
	{
		for (int k = 0; k < M.outerSize(); k++)
			for (HybridSystem::Matrix::InnerIterator it(M, k); it; ++it)
				triplets.emplace_back(firstRow + it.row(), firstCol + it.col(), it.value());
	};
	int nT = system.A_TT.rows();
	add(system.A_TT, 0, 0);
	add(system.A_TF, 0, nT);
	add(system.A_FT, nT, 0);
	add(system.A_FF, nT, nT);
	HybridSystem::Matrix M(nT + system.A_FF.rows(), nT + system.A_FF.cols());
	M.setFromTriplets(triplets.begin(), triplets.end());
	return M;
}

// Solve(b_T, b_F, x_T, x_F) after a setup from the blocks: the solution of the uncondensed system, symmetric or not
TEST(UAMGLibrary, SolvesTheUncondensedSystem)
{
	for (double asymmetry : { 0.0, 0.3 })
	{
		HybridSystem system = NonsymmetricSystem(64, asymmetry);
		Eigen::Index nT = system.A_TT.rows(), nF = system.A_FF.rows();
		Eigen::VectorXd b_T = Eigen::VectorXd::Random(nT), b_F = Eigen::VectorXd::Random(nF);
		Eigen::VectorXd b(nT + nF);
		b << b_T, b_F;
		Eigen::VectorXd exact = Eigen::SparseLU<HybridSystem::Matrix>(UncondensedMatrix(system)).solve(b);

		fhhos4::Solver uamg = ToySolver();
		if (system.Symmetric)
			uamg.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
		else
			uamg.SetupNonsymmetric(system.A_TT, system.A_TF, system.A_FT, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
		Eigen::VectorXd x_T(nT), x_F = Eigen::VectorXd::Zero(nF);
		fhhos4::Result result = uamg.Solve(b_T, b_F, x_T, x_F);
		EXPECT_TRUE(result.Converged) << asymmetry;
		EXPECT_GT(result.Iterations, 0);
		EXPECT_LT(result.Iterations, 30);
		EXPECT_LT((x_T - exact.head(nT)).norm() / exact.head(nT).norm(), 1e-8) << asymmetry;
		EXPECT_LT((x_F - exact.tail(nF)).norm() / exact.tail(nF).norm(), 1e-8) << asymmetry;

		// The condensed system of the same right-hand side, b_F - A_FT A_TT^{-1} b_T: the same faces solution
		Eigen::VectorXd condensedRHS = b_F - system.A_FT * system.A_TT.diagonal().cwiseInverse().cwiseProduct(b_T);
		Eigen::VectorXd x = Eigen::VectorXd::Zero(nF);
		uamg.Solve(condensedRHS, x);
		EXPECT_LT((x - x_F).norm() / x_F.norm(), 1e-9) << asymmetry;

		// Non-zero initial guess: x_F is the solution already
		Eigen::VectorXd y_T(nT);
		EXPECT_LE(uamg.Solve(b_T, b_F, y_T, x_F).Iterations, 1);
		EXPECT_LT((y_T - x_T).norm() / x_T.norm(), 1e-9);
	}
}

// CondensedMatrix(): the matrix given to the setup (the caller's own if read in place, converted otherwise), or the
// one computed from the blocks
TEST(UAMGLibrary, CondensedMatrix)
{
	for (double asymmetry : { 0.0, 0.3 })
	{
		HybridSystem system = NonsymmetricSystem(16, asymmetry);
		fhhos4::SparseMatrix A = system.A, A_TT = system.A_TT, A_TF = system.A_TF, A_FT = system.A_FT, A_FF = system.A_FF;
		fhhos4::Solver inPlace = ToySolver(), converted = ToySolver(), fromBlocks = ToySolver();
		EXPECT_THROW(fromBlocks.CondensedMatrix(), fhhos4::Error); // before the setup
		if (system.Symmetric)
		{
			inPlace.SetupSymmetricCondensed(A, A_TT, A_TF, system.CellInterpOfOne, system.FaceInterpOfOne);
			converted.SetupSymmetricCondensed(system.A, system.A_TT, system.A_TF, system.CellInterpOfOne, system.FaceInterpOfOne);
			fromBlocks.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
		}
		else
		{
			inPlace.SetupNonsymmetricCondensed(A, A_TT, A_TF, A_FT, system.CellInterpOfOne, system.FaceInterpOfOne);
			converted.SetupNonsymmetricCondensed(system.A, system.A_TT, system.A_TF, system.A_FT, system.CellInterpOfOne, system.FaceInterpOfOne);
			fromBlocks.SetupNonsymmetric(system.A_TT, system.A_TF, system.A_FT, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
		}
		EXPECT_EQ(&inPlace.CondensedMatrix(), &A) << asymmetry;
		EXPECT_EQ((converted.CondensedMatrix() - A).norm(), 0) << asymmetry;
		EXPECT_LT((fromBlocks.CondensedMatrix() - A).norm(), 1e-14 * A.norm()) << asymmetry;
	}
}

// Without A_FF, not used by the algorithm: same solver
TEST(UAMGLibrary, WithoutTheFaceBlock)
{
	HybridSystem system(64);
	fhhos4::Solver withA_FF = ToySolver(), withoutA_FF = ToySolver();
	withA_FF.SetupSymmetricCondensed(system.A, system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
	withoutA_FF.SetupSymmetricCondensed(system.A, system.A_TT, system.A_TF, system.CellInterpOfOne, system.FaceInterpOfOne);

	Eigen::VectorXd x1 = Eigen::VectorXd::Zero(system.b.size()), x2 = x1;
	fhhos4::Result r1 = withA_FF.Solve(system.b, x1);
	fhhos4::Result r2 = withoutA_FF.Solve(system.b, x2);
	EXPECT_EQ(withoutA_FF.NumberOfLevels(), withA_FF.NumberOfLevels());
	EXPECT_EQ(r2.Iterations, r1.Iterations);
	EXPECT_LT((x2 - x1).norm() / x1.norm(), 1e-9);
}

// The bases scaled by arbitrary factors (the constant no longer has coordinate 1): with the interpolation of 1, U-AMG
// builds the same coarse levels, and the iterates are those of the unscaled system in the scaled bases (up to the
// rounding errors). Fixed number of iterations: the stopping criterion, on the residual of each system, differs.
TEST(UAMGLibrary, InvariantToTheScalingOfTheBases)
{
	int n = 64;
	HybridSystem reference(n);
	std::mt19937 generator(1);
	std::uniform_real_distribution<double> exponent(-2, 2);
	Eigen::VectorXd sT(n * n), sF(reference.A_TF.cols());
	for (Eigen::Index i = 0; i < sT.size(); i++)
		sT[i] = std::pow(10, exponent(generator));
	for (Eigen::Index i = 0; i < sF.size(); i++)
		sF[i] = std::pow(10, exponent(generator));
	HybridSystem scaled(n, sT, sF);

	fhhos4::Solver uamgReference = ToySolver(), uamgScaled = ToySolver();
	for (fhhos4::Solver* solver : { &uamgReference, &uamgScaled })
	{
		solver->Tolerance = 1e-30; // never reached
		solver->MaxIterations = 8;
	}
	uamgReference.SetupSymmetricCondensed(reference.A, reference.A_TT, reference.A_TF, reference.A_FF, reference.CellInterpOfOne, reference.FaceInterpOfOne);
	uamgScaled.SetupSymmetricCondensed(scaled.A, scaled.A_TT, scaled.A_TF, scaled.A_FF, scaled.CellInterpOfOne, scaled.FaceInterpOfOne);
	EXPECT_EQ(uamgScaled.NumberOfLevels(), uamgReference.NumberOfLevels());

	Eigen::VectorXd x = Eigen::VectorXd::Zero(reference.b.size()), y = x;
	fhhos4::Result r = uamgReference.Solve(reference.b, x);
	fhhos4::Result s = uamgScaled.Solve(scaled.b, y);
	EXPECT_EQ(r.Iterations, 8);
	EXPECT_EQ(s.Iterations, 8);
	EXPECT_LT(RelativeError(reference, x), 1e-3); // beyond the first iterations: the comparison is meaningful
	// y = S_F^{-1} x
	EXPECT_LT((sF.cwiseProduct(y) - x).norm() / x.norm(), 1e-8);
}

// Row-major matrices with int indices, fhhos4's format: read in place if compressed, copied otherwise. Same solver as
// with the column-major ones (converted).
TEST(UAMGLibrary, RowMajorMatrices)
{
	HybridSystem system(64);
	using RowMajorMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor, int>;
	RowMajorMatrix A = system.A, A_TT = system.A_TT, A_TF = system.A_TF, A_FF = system.A_FF;
	RowMajorMatrix uncompressedA = A, uncompressedA_TT = A_TT, uncompressedA_TF = A_TF;
	for (RowMajorMatrix* M : { &uncompressedA, &uncompressedA_TT, &uncompressedA_TF })
		M->uncompress();

	fhhos4::Solver colMajor = ToySolver(), inPlace = ToySolver(), fromBlocks = ToySolver(), uncompressed = ToySolver();
	colMajor.SetupSymmetricCondensed(system.A, system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
	inPlace.SetupSymmetricCondensed(A, A_TT, A_TF, A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
	fromBlocks.SetupSymmetric(A_TT, A_TF, A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
	uncompressed.SetupSymmetricCondensed(uncompressedA, uncompressedA_TT, uncompressedA_TF, system.CellInterpOfOne, system.FaceInterpOfOne);

	Eigen::VectorXd x = Eigen::VectorXd::Zero(system.b.size());
	fhhos4::Result reference = colMajor.Solve(system.b, x);
	for (fhhos4::Solver* solver : { &inPlace, &fromBlocks, &uncompressed })
	{
		Eigen::VectorXd y = Eigen::VectorXd::Zero(system.b.size());
		EXPECT_EQ(solver->NumberOfLevels(), colMajor.NumberOfLevels());
		EXPECT_EQ(solver->Solve(system.b, y).Iterations, reference.Iterations);
		EXPECT_LT((y - x).norm() / x.norm(), 1e-9);
	}
}

// Whether SetupSymmetricCondensed() accepts an argument A of type AType: const Matrix& (a variable), Matrix (a temporary)
template <typename AType, typename = void>
struct SetupAccepts : std::false_type {};
template <typename AType>
struct SetupAccepts<AType, std::void_t<decltype(std::declval<fhhos4::Solver&>().SetupSymmetricCondensed(std::declval<AType>(), std::declval<const HybridSystem::Matrix&>(), std::declval<const HybridSystem::Matrix&>(), Eigen::VectorXd(), Eigen::VectorXd()))>> : std::true_type {};
static_assert(SetupAccepts<const HybridSystem::Matrix&>::value);
static_assert(!SetupAccepts<HybridSystem::Matrix>::value, "the solver keeps a reference to A: a temporary A must not compile");

// BiCGSTAB preconditioned by a V-cycle (linear)
TEST(UAMGLibrary, BiCGSTAB)
{
	HybridSystem system(64);
	fhhos4::Solver uamg = ToySolver();
	uamg.Krylov = "bicgstab";
	uamg.Cycle = 'V';
	uamg.SetupSymmetricCondensed(system.A, system.A_TT, system.A_TF, system.CellInterpOfOne, system.FaceInterpOfOne);
	Eigen::VectorXd x = Eigen::VectorXd::Zero(system.b.size());
	fhhos4::Result result = uamg.Solve(system.b, x);
	EXPECT_TRUE(result.Converged);
	EXPECT_LT(result.Iterations, 30);
	EXPECT_LT(RelativeError(system, x), 1e-8);
}

// GCR preconditioned by the K-cycle, which then iterates with GCR on the coarse levels: for any matrix, here symmetric
TEST(UAMGLibrary, GCR)
{
	HybridSystem system(64);
	fhhos4::Solver uamg = ToySolver();
	uamg.Krylov = "gcr";
	uamg.SetupSymmetricCondensed(system.A, system.A_TT, system.A_TF, system.CellInterpOfOne, system.FaceInterpOfOne);
	Eigen::VectorXd x = Eigen::VectorXd::Zero(system.b.size());
	fhhos4::Result result = uamg.Solve(system.b, x);
	EXPECT_TRUE(result.Converged);
	EXPECT_LT(result.Iterations, 30);
	EXPECT_LT(RelativeError(system, x), 1e-8);
}

// The default smoothers: the hybrid Gauss-Seidel (parallel) with FCG on several threads, the sequential one otherwise
TEST(UAMGLibrary, DefaultSmoothers)
{
	HybridSystem system(64);
	for (int threads : { 4, 1 })
	{
		fhhos4::Solver uamg = ToySolver();
		uamg.PreSmoother = "";
		uamg.PostSmoother = "";
		uamg.Threads = threads;
		uamg.Verbosity = 2; // prints the smoothers
		::testing::internal::CaptureStdout();
		uamg.SetupSymmetricCondensed(system.A, system.A_TT, system.A_TF, system.CellInterpOfOne, system.FaceInterpOfOne);
		Eigen::VectorXd x = Eigen::VectorXd::Zero(system.b.size());
		fhhos4::Result result = uamg.Solve(system.b, x);
		std::string output = ::testing::internal::GetCapturedStdout();
		EXPECT_EQ(output.find("Hybrid Gauss-Seidel (forward)") != std::string::npos, threads > 1) << threads << " threads";
		EXPECT_EQ(output.find("Hybrid Gauss-Seidel (backward)") != std::string::npos, threads > 1) << threads << " threads";
		EXPECT_NE(output.find("Gauss-Seidel (backward)"), std::string::npos);
		EXPECT_TRUE(result.Converged);
		EXPECT_LT(result.Iterations, 30);
		EXPECT_LT(RelativeError(system, x), 1e-8);
	}
}

// U-AMG alone (no FCG)
TEST(UAMGLibrary, MultigridAlone)
{
	HybridSystem system(64);
	fhhos4::Solver uamg = ToySolver();
	uamg.Krylov = "none";
	uamg.SetupSymmetricCondensed(system.A, system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
	Eigen::VectorXd x = uamg.Solve(system.b);
	EXPECT_LT(RelativeError(system, x), 1e-8);
}

// Non-symmetric system (SetupNonsymmetricCondensed()): by default GCR, the K-cycle with GCR, LU on the coarsest level.
// Also BiCGSTAB with a V-cycle, and the blocks only (SetupNonsymmetric(): the library computes A).
TEST(UAMGLibrary, NonsymmetricSystem)
{
	for (double asymmetry : { 0.3, 0.9 })
	{
		HybridSystem system = NonsymmetricSystem(64, asymmetry);
		ASSERT_GT((system.A - HybridSystem::Matrix(system.A.transpose())).norm(), 0.1 * system.A.norm());

		fhhos4::Solver gcr = ToySolver(), bicgstab = ToySolver(), fromBlocks = ToySolver();
		bicgstab.Krylov = "bicgstab";
		bicgstab.Cycle = 'V';
		gcr.SetupNonsymmetricCondensed(system.A, system.A_TT, system.A_TF, system.A_FT, system.CellInterpOfOne, system.FaceInterpOfOne);
		bicgstab.SetupNonsymmetricCondensed(system.A, system.A_TT, system.A_TF, system.A_FT, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
		fromBlocks.SetupNonsymmetric(system.A_TT, system.A_TF, system.A_FT, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
		EXPECT_GT(gcr.NumberOfLevels(), 2);

		Eigen::VectorXd x = Eigen::VectorXd::Zero(system.b.size()), y = x, z = x;
		fhhos4::Result r = gcr.Solve(system.b, x);
		fhhos4::Result s = bicgstab.Solve(system.b, y);
		fhhos4::Result t = fromBlocks.Solve(system.b, z);
		for (const fhhos4::Result& result : { r, s, t })
		{
			EXPECT_TRUE(result.Converged);
			EXPECT_LT(result.Iterations, 30);
		}
		EXPECT_LT(RelativeError(system, x), 1e-8);
		EXPECT_LT(RelativeError(system, y), 1e-8);
		EXPECT_EQ(t.Iterations, r.Iterations);
		EXPECT_LT((z - x).norm() / x.norm(), 1e-9);
	}
}

// SetupNonsymmetricCondensed() on a symmetric system (A_FT = A_TF^T): the coarsening of its symmetric part, the same
// as SetupSymmetricCondensed()'s, and Galerkin coarse operators equal to the symmetric ones up to rounding: the same
// iterations with GCR
TEST(UAMGLibrary, NonsymmetricSetupOfASymmetricSystem)
{
	HybridSystem system(64);
	fhhos4::Solver symmetric = ToySolver(), nonsymmetric = ToySolver();
	symmetric.Krylov = "gcr";
	symmetric.SetupSymmetricCondensed(system.A, system.A_TT, system.A_TF, system.CellInterpOfOne, system.FaceInterpOfOne);
	nonsymmetric.SetupNonsymmetricCondensed(system.A, system.A_TT, system.A_TF, system.A_FT, system.CellInterpOfOne, system.FaceInterpOfOne);
	EXPECT_EQ(nonsymmetric.NumberOfLevels(), symmetric.NumberOfLevels());

	Eigen::VectorXd x = Eigen::VectorXd::Zero(system.b.size()), y = x;
	fhhos4::Result r = symmetric.Solve(system.b, x);
	fhhos4::Result s = nonsymmetric.Solve(system.b, y);
	EXPECT_EQ(s.Iterations, r.Iterations);
	EXPECT_LT((y - x).norm() / x.norm(), 1e-9);
}

// fhhos4::Version(), compiled in the library, starts with the version of the header: both come from project() in
// CMakeLists.txt.
TEST(UAMGLibrary, Version)
{
	std::string version = fhhos4::Version();
	EXPECT_EQ(version.rfind(FHHOS4_VERSION, 0), 0u) << version;
	EXPECT_EQ(std::string(FHHOS4_VERSION), std::to_string(FHHOS4_VERSION_MAJOR) + "." +
		std::to_string(FHHOS4_VERSION_MINOR) + "." + std::to_string(FHHOS4_VERSION_PATCH));
}

TEST(UAMGLibrary, Errors)
{
	HybridSystem system(16);
	fhhos4::Solver uamg = ToySolver();
	Eigen::VectorXd x = Eigen::VectorXd::Zero(system.b.size());
	EXPECT_THROW(uamg.Solve(system.b, x), fhhos4::Error); // before the setup

	// The uncondensed system needs a setup from its blocks, and vectors of the sizes of the cells and faces DoFs
	Eigen::VectorXd b_T = Eigen::VectorXd::Ones(system.A_TT.rows()), x_T(system.A_TT.rows()), b_F = Eigen::VectorXd::Zero(system.b.size());
	Eigen::VectorXd tooShortT = b_T.head(3), tooShortF = b_F.head(3);
	fhhos4::Solver condensed = ToySolver(), uncondensed = ToySolver();
	condensed.SetupSymmetricCondensed(system.A, system.A_TT, system.A_TF, system.CellInterpOfOne, system.FaceInterpOfOne);
	EXPECT_THROW(condensed.Solve(b_T, b_F, x_T, x), fhhos4::Error);
	uncondensed.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
	EXPECT_THROW(uncondensed.Solve(tooShortT, b_F, x_T, x), fhhos4::Error);
	EXPECT_THROW(uncondensed.Solve(b_T, b_F, tooShortT, x), fhhos4::Error);
	EXPECT_THROW(uncondensed.Solve(b_T, tooShortF, x_T, x), fhhos4::Error);
	EXPECT_THROW(uncondensed.Solve(b_T, b_F, x_T, tooShortF), fhhos4::Error);

	// Interpolations of 1 of the wrong size, or with a zero coefficient
	Eigen::VectorXd tooShort = system.FaceInterpOfOne.head(10);
	EXPECT_THROW(uamg.SetupSymmetricCondensed(system.A, system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, tooShort), fhhos4::Error);
	Eigen::VectorXd zero = system.CellInterpOfOne;
	zero[3] = 0;
	EXPECT_THROW(uamg.SetupSymmetricCondensed(system.A, system.A_TT, system.A_TF, system.A_FF, zero, system.FaceInterpOfOne), fhhos4::Error);
	EXPECT_FALSE(uamg.IsSetUp());

	// Blocks of inconsistent sizes
	HybridSystem other(8);
	EXPECT_THROW(uamg.SetupSymmetricCondensed(system.A, other.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);

	// Wrong degrees: 3 DoFs per cell, or 2 per face, in 2D
	fhhos4::Solver cellDegree1 = ToySolver();
	cellDegree1.CellDegree = 1;
	EXPECT_THROW(cellDegree1.SetupSymmetricCondensed(system.A, system.A_TT, system.A_TF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);
	fhhos4::Solver faceDegree1 = ToySolver();
	faceDegree1.FaceDegree = 1;
	EXPECT_THROW(faceDegree1.SetupSymmetricCondensed(system.A, system.A_TT, system.A_TF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);

	// Parameters not set, or invalid
	fhhos4::Solver noDimension = ToySolver();
	noDimension.Dimension = -1;
	EXPECT_THROW(noDimension.SetupSymmetricCondensed(system.A, system.A_TT, system.A_TF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);
	fhhos4::Solver noFaceDegree = ToySolver();
	noFaceDegree.FaceDegree = -1;
	EXPECT_THROW(noFaceDegree.SetupSymmetricCondensed(system.A, system.A_TT, system.A_TF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);
	fhhos4::Solver noCellDegree = ToySolver();
	noCellDegree.CellDegree = -1;
	EXPECT_THROW(noCellDegree.SetupSymmetricCondensed(system.A, system.A_TT, system.A_TF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);
	fhhos4::Solver wrongCycle = ToySolver();
	wrongCycle.Cycle = 'X';
	EXPECT_THROW(wrongCycle.SetupSymmetricCondensed(system.A, system.A_TT, system.A_TF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);
	fhhos4::Solver wrongKrylov = ToySolver();
	wrongKrylov.Krylov = "gmres";
	EXPECT_THROW(wrongKrylov.SetupSymmetricCondensed(system.A, system.A_TT, system.A_TF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);
	fhhos4::Solver bicgstabKCycle = ToySolver(); // BiCGSTAB is not flexible: the K-cycle is not linear
	bicgstabKCycle.Krylov = "bicgstab";
	EXPECT_THROW(bicgstabKCycle.SetupSymmetricCondensed(system.A, system.A_TT, system.A_TF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);

	// Non-symmetric matrix: FCG or Cholesky refused, A_FT of the wrong size
	fhhos4::Solver fcg = ToySolver();
	fcg.Krylov = "fcg";
	EXPECT_THROW(fcg.SetupNonsymmetricCondensed(system.A, system.A_TT, system.A_TF, system.A_FT, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);
	fhhos4::Solver cholesky = ToySolver();
	cholesky.CoarseSolver = "ch";
	EXPECT_THROW(cholesky.SetupNonsymmetricCondensed(system.A, system.A_TT, system.A_TF, system.A_FT, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);
	EXPECT_THROW(uamg.SetupNonsymmetricCondensed(system.A, system.A_TT, system.A_TF, system.A_TF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);

	// A_TT not block diagonal
	HybridSystem::Matrix A_TT = system.A_TT;
	A_TT.coeffRef(1, 0) = -1;
	EXPECT_THROW(uamg.SetupSymmetricCondensed(system.A, A_TT, system.A_TF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);

	// Vector of the wrong size
	uamg.SetupSymmetricCondensed(system.A, system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
	Eigen::VectorXd wrongSize = Eigen::VectorXd::Zero(3);
	EXPECT_THROW(uamg.Solve(system.b, wrongSize), fhhos4::Error);
}

// Verbosity 0: nothing printed. The number of threads and the format of cout are restored after each call.
TEST(UAMGLibrary, LeavesTheCallerUnchanged)
{
	HybridSystem system(32);
	fhhos4::Solver uamg = ToySolver();
	uamg.Threads = 1;
#ifdef _OPENMP
	int threadsBefore = omp_get_max_threads();
#endif

	::testing::internal::CaptureStdout();
	uamg.SetupSymmetricCondensed(system.A, system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
	Eigen::VectorXd x = uamg.Solve(system.b);
	std::string output = ::testing::internal::GetCapturedStdout();
	EXPECT_EQ(output, "");

#ifdef _OPENMP
	EXPECT_EQ(omp_get_max_threads(), threadsBefore);
#endif

	// Verbosity 2: fhhos4's messages, with the iterations, which change the precision of cout
	fhhos4::Solver verbose = ToySolver();
	verbose.Verbosity = 2;
	std::streamsize precision = std::cout.precision(7);
	std::ios_base::fmtflags flags = std::cout.flags();
	::testing::internal::CaptureStdout();
	verbose.SetupSymmetricCondensed(system.A, system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
	verbose.Solve(system.b);
	output = ::testing::internal::GetCapturedStdout();
	EXPECT_NE(output.find("Setup"), std::string::npos);
	EXPECT_EQ(std::cout.precision(), 7);
	EXPECT_EQ(std::cout.flags(), flags);
#ifdef _OPENMP
	EXPECT_EQ(omp_get_max_threads(), threadsBefore); // Threads = 0: at most one per physical core during the calls only
#endif
	std::cout.precision(precision);
}
