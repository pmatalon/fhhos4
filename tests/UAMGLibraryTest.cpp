// Tests of the library fhhos4 (library/), U-AMG for the codes that discretize the problem themselves.
// - UAMGLibraryProgram: the program runs U-AMG through the library (-s libuamg, fcglibuamg, gcrlibuamg), as an external
//   code calls it, and must get the iteration counts and solutions of the U-AMG compiled in the program (uamg, fcguamg,
//   gcruamg), including with L2-orthonormal bases, whose constant does not have coordinate 1. The program gives the
//   library the uncondensed system (the library condenses it and computes the cell unknowns), or, when it needs the
//   condensed system (the export of the linear system), solves it by the library on its condensed matrix (Solve(b, x)).
// - UAMGMixedPrecisionProgram: mixed precision (-uamg-mixed: the levels in single precision), in the program and through
//   the library.
// - UAMGBlockCSRProgram: the operators in block CSR in double (-uamg-bcsr), in the program and through the library.
// - UAMGLibrary: the public API alone (<fhhos4>), on a hybrid system built here.
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
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
	P_CoarsStgy PCoarsening = P_CoarsStgy::Minus2; // -p-cs
	int CoarseningProlong = 0;                     // -coarsening-prolong (0: the default)
	bool PolynomialSpaceQ = false;                 // -poly-space q
};

inline std::ostream& operator<<(std::ostream& os, const LibraryCase& c)
{
	return os << c.Name;
}

// The program on the case with the solver solverCode. nonsymmetric: -uamg-nonsym; exportSystem: -export lsys (the program
// solves the condensed system), into a temporary directory; chebyshev: the Chebyshev smoothers of the fourth kind,
// degree 2 (the library's default for a symmetric matrix), instead of the program's; mixedPrecision: -uamg-mixed;
// blockCSR: -uamg-bcsr.
ProgramResults RunCase(const LibraryCase& c, const std::string& solverCode, bool nonsymmetric = false, bool exportSystem = false,
                       bool chebyshev = false, bool mixedPrecision = false, bool blockCSR = false)
{
	ProgramArguments args;
	args.Problem.GeoCode = c.Geo;
	args.Discretization.MeshCode = c.Mesh;
	if (!c.Mesher.empty())
		args.Discretization.Mesher = c.Mesher;
	args.Discretization.N = c.N;
	args.Discretization.PolyDegree = c.K + 1;
	args.Discretization.UsePolynomialSpaceQ = c.PolynomialSpaceQ;
	if (c.BasesOrthogonalization != -1)
	{
		args.Discretization.OrthogonalizeElemBasesCode = c.BasesOrthogonalization;
		args.Discretization.OrthogonalizeFaceBasesCode = c.BasesOrthogonalization;
	}
	args.Solver.SolverCode = solverCode;
	args.Solver.MG.P_CS = c.PCoarsening;
	args.Solver.MG.CoarseningProlongationCode = c.CoarseningProlong;
	args.Solver.MG.UAMGNonsymmetric = nonsymmetric;
	if (exportSystem)
	{
		args.Actions.Export.LinearSystem = true;
		args.OutputDirectory = ::testing::TempDir();
	}
	if (chebyshev)
	{
		args.Solver.MG.PreSmoother = SmootherType::BlockChebyshevFourthKind;
		args.Solver.MG.PostSmoother = SmootherType::Default;
		args.Solver.MG.PreSmoothingIterations = 2;
		args.Solver.MG.PostSmoothingIterations = 2;
	}
	args.Solver.MG.UAMGMixedPrecision = mixedPrecision;
	args.Solver.MG.UAMGBlockCSRStorage = blockCSR;
	args.Actions.UseCache = false;
	return RunDiffusionHHO(args);
}

class UAMGLibraryProgramTest : public ::testing::TestWithParam<std::tuple<LibraryCase, std::string>>
{
};

// The library code (libuamg, fcglibuamg, gcrlibuamg) against the program's (uamg, fcguamg, gcruamg), on the uncondensed
// system and on the condensed one (exportSystem)
TEST_P(UAMGLibraryProgramTest, SameAsProgram)
{
	auto [c, libraryCode] = GetParam();
	std::string programCode = libraryCode;
	programCode.replace(programCode.find("libuamg"), 7, "uamg");

	ProgramResults program = RunCase(c, programCode);
	EXPECT_GT(program.IterationCount, 0);
	EXPECT_LT(program.IterationCount, 200); // converged (default maximum: 200)
	for (bool exportSystem : { false, true })
	{
		ProgramResults library = RunCase(c, libraryCode, false, exportSystem);
		std::string label = libraryCode + (exportSystem ? " -export lsys" : "") + " vs " + programCode;
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
	// The p-coarsening straight to degree 0 (PCoarsening::ToZero), which differs from the default from k = 3
	cases.push_back({ "Square_tri_n16_k3_pcs0", "square", "tri", "", 16, 3, -1, P_CoarsStgy::DirectToLow });
	// The coarsening prolongation HarmonicExtensionOrInject (Prolongation::HarmonicExtension in the library)
	cases.push_back({ "Square_tri_n16_k1_orthonormal_harmonic", "square", "tri", "", 16, 1, 3, P_CoarsStgy::Minus2, (int)UAMGProlongation::HarmonicExtensionOrInject });
	// The space Q (-poly-space q, PolynomialSpace::Q in the library): blocks of Q^k, and a p-level Q^1 at k = 3
	cases.push_back({ "Square_cart_n16_k3_Q", "square", "cart", "inhouse", 16, 3, -1, P_CoarsStgy::Minus2, 0, true });
#endif
#ifdef ENABLE_3D
	cases.push_back({ "Cube_tetra_n8_k0", "cube", "tetra", "", 8, 0, -1 });
	cases.push_back({ "Cube_tetra_n8_k0_orthonormal", "cube", "tetra", "", 8, 0, 3 });
	cases.push_back({ "Cube_tetra_n8_k1_orthonormal", "cube", "tetra", "", 8, 1, 3 });
	cases.push_back({ "Cube_tetra_n8_k0_orthonormal_harmonic", "cube", "tetra", "", 8, 0, 3, P_CoarsStgy::Minus2, (int)UAMGProlongation::HarmonicExtensionOrInject });
	cases.push_back({ "Cube_cart_n8_k2", "cube", "cart", "inhouse", 8, 2, -1 });
	cases.push_back({ "Cube_cart_n8_k2_Q", "cube", "cart", "inhouse", 8, 2, -1, P_CoarsStgy::Minus2, 0, true });
#endif
	return cases;
}

// Q^3 in 3D: faces of 16 DoFs (the kernels of block CSR for this size, in double and in single precision), and a p-level
// Q^1 between Q^3 and Q^0. Only in the suites of block CSR, with one Krylov method: the assembly of the program is
// long (the reconstruction in Q^3 + P^4: 2.5 s for 125 cells, the solve 0.01 s). n = 5: on the mesh n = 4, the traces
// of the solution vanish on all the faces (the condensed right-hand side is rounding noise, which the condensations of
// the program and of the library don't round alike).
std::vector<LibraryCase> HighDegreeQCases()
{
	return { { "Cube_cart_n5_k3_Q", "cube", "cart", "inhouse", 5, 3, -1, P_CoarsStgy::Minus2, 0, true } };
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
// through the library (gcrlibuamg: SetupNonsymmetric(), on the uncondensed and the condensed system).
class UAMGNonsymmetricProgramTest : public ::testing::TestWithParam<LibraryCase>
{
};

TEST_P(UAMGNonsymmetricProgramTest, SameAsSymmetric)
{
	LibraryCase c = GetParam();
	ProgramResults symmetric = RunCase(c, "gcruamg");
	EXPECT_GT(symmetric.IterationCount, 0);
	EXPECT_LT(symmetric.IterationCount, 200);
	for (auto [code, exportSystem] : { std::pair{ "gcruamg", false }, { "gcrlibuamg", false }, { "gcrlibuamg", true } })
	{
		ProgramResults nonsymmetric = RunCase(c, code, true, exportSystem);
		std::string label = std::string(code) + " -uamg-nonsym" + (exportSystem ? " -export lsys" : "");
		EXPECT_EQ(nonsymmetric.IterationCount, symmetric.IterationCount) << label;
		EXPECT_NEAR(nonsymmetric.L2Error, symmetric.L2Error, 1e-6 * symmetric.L2Error) << label;
	}
}

INSTANTIATE_TEST_SUITE_P(UAMGNonsymmetricProgram, UAMGNonsymmetricProgramTest, ::testing::ValuesIn(LibraryCases()),
	[](const ::testing::TestParamInfo<LibraryCase>& info) { return info.param.Name; });

// Mixed precision (-uamg-mixed: the levels in single precision; with the Chebyshev smoothers, which have a single
// precision version): the same iterations and solutions in the program (uamg, fcguamg, gcruamg) and through the library
// (libuamg..., on the uncondensed and the condensed system), and those in double: the same levels (the setup in double), the iterations within 1 (equal here), the
// same L2 error up to the tolerance
class UAMGMixedPrecisionProgramTest : public ::testing::TestWithParam<std::tuple<LibraryCase, std::string>>
{
};

TEST_P(UAMGMixedPrecisionProgramTest, SameInProgramAndLibrary)
{
	auto [c, programCode] = GetParam();
	std::string libraryCode = programCode;
	libraryCode.replace(libraryCode.find("uamg"), 4, "libuamg");

	ProgramResults doublePrecision = RunCase(c, programCode, false, false, true, false);
	ProgramResults program = RunCase(c, programCode, false, false, true, true);
	EXPECT_LT(program.IterationCount, 200); // converged (default maximum: 200)
	EXPECT_NEAR(program.IterationCount, doublePrecision.IterationCount, 1) << programCode;
	EXPECT_NEAR(program.L2Error, doublePrecision.L2Error, 1e-6 * doublePrecision.L2Error) << programCode;
	for (bool exportSystem : { false, true })
	{
		ProgramResults library = RunCase(c, libraryCode, false, exportSystem, true, true);
		std::string label = libraryCode + (exportSystem ? " -export lsys" : "") + " vs " + programCode + ", -uamg-mixed";
		EXPECT_EQ(library.IterationCount, program.IterationCount) << label;
		EXPECT_NEAR(library.L2Error, program.L2Error, 1e-6 * program.L2Error) << label;
	}
}

#ifdef ENABLE_3D
INSTANTIATE_TEST_SUITE_P(UAMGMixedPrecisionProgramQ, UAMGMixedPrecisionProgramTest,
	::testing::Combine(::testing::ValuesIn(HighDegreeQCases()), ::testing::Values("fcguamg")),
	[](const ::testing::TestParamInfo<UAMGMixedPrecisionProgramTest::ParamType>& info)
	{
		return std::get<0>(info.param).Name + "_" + std::get<1>(info.param);
	});
#endif

INSTANTIATE_TEST_SUITE_P(UAMGMixedPrecisionProgram, UAMGMixedPrecisionProgramTest,
	::testing::Combine(::testing::ValuesIn(LibraryCases()), ::testing::Values("fcguamg", "gcruamg", "uamg")),
	[](const ::testing::TestParamInfo<UAMGMixedPrecisionProgramTest::ParamType>& info)
	{
		return std::get<0>(info.param).Name + "_" + std::get<1>(info.param);
	});

// The operators in block CSR in double (-uamg-bcsr; with the Chebyshev smoothers, which support it): the same iterations
// and solutions as in CSR (to rounding: Lanczos' estimates of lambdaMax), in the program (fcguamg) and through the
// library (fcglibuamg, on the uncondensed and the condensed system)
class UAMGBlockCSRProgramTest : public ::testing::TestWithParam<LibraryCase>
{
};

TEST_P(UAMGBlockCSRProgramTest, SameAsCSR)
{
	LibraryCase c = GetParam();
	ProgramResults csr = RunCase(c, "fcguamg", false, false, true);
	EXPECT_LT(csr.IterationCount, 200);
	for (auto [code, exportSystem] : { std::pair{ "fcguamg", false }, { "fcglibuamg", false }, { "fcglibuamg", true } })
	{
		ProgramResults blockCSR = RunCase(c, code, false, exportSystem, true, false, true);
		std::string label = std::string(code) + " -uamg-bcsr" + (exportSystem ? " -export lsys" : "");
		EXPECT_EQ(blockCSR.IterationCount, csr.IterationCount) << label;
		EXPECT_NEAR(blockCSR.L2Error, csr.L2Error, 1e-10 * csr.L2Error) << label;
	}
}

INSTANTIATE_TEST_SUITE_P(UAMGBlockCSRProgram, UAMGBlockCSRProgramTest, ::testing::ValuesIn(LibraryCases()),
	[](const ::testing::TestParamInfo<LibraryCase>& info) { return info.param.Name; });
#ifdef ENABLE_3D
INSTANTIATE_TEST_SUITE_P(UAMGBlockCSRProgramQ, UAMGBlockCSRProgramTest, ::testing::ValuesIn(HighDegreeQCases()),
	[](const ::testing::TestParamInfo<LibraryCase>& info) { return info.param.Name; });
#endif

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

#endif

//------------------------------------------------------------------------------------------------------------------//
//                                     Through the public API only (<fhhos4>)                                    //
//------------------------------------------------------------------------------------------------------------------//

// Hybrid system of degree 0 on the n x n Cartesian mesh of the unit square: in each cell T, of coefficient a_T, the
// local matrix couples the cell DoF to each of its 4 faces with the weight a_T (the constants are in its kernel).
// The DoFs of the boundary faces (Dirichlet) are eliminated. Matrices in the format of the library
// (fhhos4::SparseMatrix: row-major, its index type), read in place.
// Non-symmetric with asymmetry > 0, like an advection in the direction (1, 1): the cell-face coupling is
// -a_T (1 + asymmetry s) in A_TF and -a_T (1 - asymmetry s) in A_FT, s = 1 on the right and top faces of the cell (the
// outflow), -1 on the others. The symmetric part is the symmetric system: positive definite.
struct HybridSystem
{
	using Matrix = fhhos4::SparseMatrix;
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

// The solver for HybridSystem: degree 0 in 2D. With the sequential smoothers: the defaults depend on the number of
// threads (tested apart).
fhhos4::Solver ToySolver()
{
	fhhos4::Solver solver;
#if FHHOS4_MPI
	solver.Comm = MPI_COMM_SELF;
#endif
	solver.Dimension = 2;
	solver.FaceDegree = 0;
	solver.CellDegree = 0;
	solver.Tolerance = 1e-10;
	solver.PreSmoother = fhhos4::Smoother::BlockGaussSeidel;
	solver.PostSmoother = fhhos4::Smoother::BackwardBlockGaussSeidel;
	return solver;
}

// Column-major matrices, for Eigen's direct solvers
using ColMajorMatrix = Eigen::SparseMatrix<double>;

// Relative error of x against the direct solution
double RelativeError(const HybridSystem& system, const Eigen::VectorXd& x)
{
	Eigen::VectorXd exact;
	ColMajorMatrix A = system.A;
	if (system.Symmetric)
		exact = Eigen::SimplicialLDLT<ColMajorMatrix>(A).solve(system.b);
	else
	{
		Eigen::SparseLU<ColMajorMatrix> direct(A);
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
	uamg.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
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

// The symmetric blocks given by their lower triangular parts only: the same solver
TEST(UAMGLibrary, LowerTriangularParts)
{
	HybridSystem system(64);
	fhhos4::Solver full = ToySolver();
	full.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
	Eigen::VectorXd x1 = Eigen::VectorXd::Zero(system.b.size());
	fhhos4::Result r1 = full.Solve(system.b, x1);

	HybridSystem::Matrix A_TT = system.A_TT.triangularView<Eigen::Lower>(), A_FF = system.A_FF.triangularView<Eigen::Lower>();
	A_TT.makeCompressed();
	A_FF.makeCompressed();
	fhhos4::Solver lower = ToySolver();
	lower.SetupSymmetric(A_TT, system.A_TF, A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
	Eigen::VectorXd x3 = Eigen::VectorXd::Zero(system.b.size());
	EXPECT_EQ(lower.Solve(system.b, x3).Iterations, r1.Iterations);
	EXPECT_LT((x3 - x1).norm() / x1.norm(), 1e-9);
}

// The matrix of the uncondensed system [A_TT A_TF; A_FT A_FF] (column-major, for Eigen's direct solvers)
ColMajorMatrix UncondensedMatrix(const HybridSystem& system)
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
	ColMajorMatrix M(nT + system.A_FF.rows(), nT + system.A_FF.cols());
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
		Eigen::VectorXd exact = Eigen::SparseLU<ColMajorMatrix>(UncondensedMatrix(system)).solve(b);

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

// CondensedMatrix(): the condensed matrix computed from the blocks
TEST(UAMGLibrary, CondensedMatrix)
{
	for (double asymmetry : { 0.0, 0.3 })
	{
		HybridSystem system = NonsymmetricSystem(16, asymmetry);
		fhhos4::Solver uamg = ToySolver();
		EXPECT_THROW(uamg.CondensedMatrix(), fhhos4::Error); // before the setup
		if (system.Symmetric)
			uamg.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
		else
			uamg.SetupNonsymmetric(system.A_TT, system.A_TF, system.A_FT, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
		EXPECT_LT((uamg.CondensedMatrix() - system.A).norm(), 1e-14 * system.A.norm()) << asymmetry;
	}
}

// The bases scaled by arbitrary factors (the constant no longer has coordinate 1): with the interpolation of 1, U-AMG
// builds the same coarse levels, and the iterates are those of the unscaled system in the scaled bases (up to the
// rounding errors). Fixed number of iterations: the stopping criterion, on the residual of each system, differs.
// The scaling of the bases given by the interpolations of 1: the same aggregates (the couplings of the rescaled
// matrices), and the same iterates (the scaled ones)
void ExpectInvariantToTheScalingOfTheBases(fhhos4::Prolongation prolongation)
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
		solver->Prolongation = prolongation;
	}
	uamgReference.SetupSymmetric(reference.A_TT, reference.A_TF, reference.A_FF, reference.CellInterpOfOne, reference.FaceInterpOfOne);
	uamgScaled.SetupSymmetric(scaled.A_TT, scaled.A_TF, scaled.A_FF, scaled.CellInterpOfOne, scaled.FaceInterpOfOne);
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

TEST(UAMGLibrary, InvariantToTheScalingOfTheBases)
{
	ExpectInvariantToTheScalingOfTheBases(fhhos4::Prolongation::Reconstruction);
}

// Prolongation::HarmonicExtension (the default) converts no cell value into face values: invariant too, and it solves
// the system
TEST(UAMGLibrary, HarmonicExtension)
{
	EXPECT_TRUE(fhhos4::Solver().Prolongation == fhhos4::Prolongation::HarmonicExtension);
	ExpectInvariantToTheScalingOfTheBases(fhhos4::Prolongation::HarmonicExtension);

	HybridSystem system(64);
	fhhos4::Solver reconstruction = ToySolver(), harmonic = ToySolver();
	harmonic.Prolongation = fhhos4::Prolongation::HarmonicExtension;
	for (fhhos4::Solver* solver : { &reconstruction, &harmonic })
		solver->SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
	EXPECT_EQ(harmonic.NumberOfLevels(), reconstruction.NumberOfLevels()); // the same aggregates
	Eigen::VectorXd x = Eigen::VectorXd::Zero(system.b.size());
	fhhos4::Result result = harmonic.Solve(system.b, x);
	EXPECT_TRUE(result.Converged);
	EXPECT_LT(result.Iterations, 30);
	EXPECT_LT(RelativeError(system, x), 1e-8);
}

// The matrices in other formats: column-major ones converted by Eigen (a copy), or viewed by TransposedView() (no copy:
// symmetric matrices), 64-bit indices converted by Eigen, Eigen::Map of the arrays: the same solver as with the matrices
// in the format of the library. Uncompressed matrices refused.
TEST(UAMGLibrary, OtherFormats)
{
	HybridSystem system(64);
	ColMajorMatrix colMajorA_TT = system.A_TT, colMajorA_TF = system.A_TF, colMajorA_FF = system.A_FF;
	fhhos4::SparseMatrix A_TT = colMajorA_TT, A_TF = colMajorA_TF;
	using LongMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor, std::int64_t>;
	LongMatrix longA_FF = system.A_FF;
	fhhos4::SparseMatrix A_FF = longA_FF;
	EXPECT_EQ((A_TF - system.A_TF).norm(), 0);
	EXPECT_EQ((A_FF - system.A_FF).norm(), 0);
	const HybridSystem::Matrix& M = system.A_FF;
	fhhos4::SparseMatrixMap mapA_FF(M.rows(), M.cols(), M.nonZeros(), M.outerIndexPtr(), M.innerIndexPtr(), M.valuePtr());

	fhhos4::Solver native = ToySolver(), converted = ToySolver(), transposedViews = ToySolver(), maps = ToySolver();
	native.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
	converted.SetupSymmetric(A_TT, A_TF, A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
	maps.SetupSymmetric(system.A_TT, system.A_TF, mapA_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
	// The symmetric ones viewed in place (A_TF is not square: converted)
	transposedViews.SetupSymmetric(fhhos4::TransposedView(colMajorA_TT), A_TF, fhhos4::TransposedView(colMajorA_FF), system.CellInterpOfOne, system.FaceInterpOfOne);

	Eigen::VectorXd x = Eigen::VectorXd::Zero(system.b.size());
	fhhos4::Result reference = native.Solve(system.b, x);
	for (fhhos4::Solver* solver : { &converted, &transposedViews, &maps })
	{
		Eigen::VectorXd y = Eigen::VectorXd::Zero(system.b.size());
		EXPECT_EQ(solver->NumberOfLevels(), native.NumberOfLevels());
		EXPECT_EQ(solver->Solve(system.b, y).Iterations, reference.Iterations);
		EXPECT_LT((y - x).norm() / x.norm(), 1e-9);
	}

	HybridSystem::Matrix uncompressed = system.A_FF;
	uncompressed.uncompress();
	ColMajorMatrix uncompressedColMajor = colMajorA_FF;
	uncompressedColMajor.uncompress();
	fhhos4::Solver solver = ToySolver();
	EXPECT_THROW(solver.SetupSymmetric(system.A_TT, system.A_TF, uncompressed, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);
	EXPECT_THROW(fhhos4::TransposedView(uncompressedColMajor), fhhos4::Error);
}

// The global matrix of the system, in its own numbering (here: the hybrid system with the DoFs of the boundary faces,
// scrambled), with 64-bit indices: ExtractSymmetricBlocks() and ExtractNonsymmetricBlocks() give the blocks of the lists
// (also the compiled ones that the Python module calls, with 64-bit row pointers and 32-bit columns), Gather() and
// Scatter() the vectors: the same solver and solution as from the blocks.
TEST(UAMGLibrary, BlocksOfTheGlobalMatrix)
{
	for (double asymmetry : { 0.0, 0.3 })
	{
		HybridSystem system = NonsymmetricSystem(32, asymmetry);
		std::int64_t nT = system.A_TT.rows(), nF = system.A_FF.rows(), nD = 50; // nD Dirichlet DoFs, coupled to everything
		std::int64_t n = nT + nF + nD;
		// The global numbering: a random permutation of (cells, faces, Dirichlet DoFs)
		std::vector<std::int64_t> global(n);
		for (std::int64_t i = 0; i < n; i++)
			global[i] = i;
		std::shuffle(global.begin(), global.end(), std::mt19937(7));
		std::vector<Eigen::Triplet<double>> triplets;
		auto add = [&](const HybridSystem::Matrix& M, std::int64_t firstRow, std::int64_t firstCol)
		{
			for (int k = 0; k < M.outerSize(); k++)
				for (HybridSystem::Matrix::InnerIterator it(M, k); it; ++it)
					triplets.emplace_back(global[firstRow + it.row()], global[firstCol + it.col()], it.value());
		};
		add(system.A_TT, 0, 0);
		add(system.A_TF, 0, nT);
		add(system.A_FT, nT, 0);
		add(system.A_FF, nT, nT);
		for (std::int64_t d = nT + nF; d < n; d++) // the Dirichlet DoFs: coupled to themselves and to a cell and a face
		{
			triplets.emplace_back(global[d], global[d], 1.0);
			triplets.emplace_back(global[d], global[d % nT], 0.5);
			triplets.emplace_back(global[d % nT], global[d], 0.5);
			triplets.emplace_back(global[d], global[nT + d % nF], 0.5);
			triplets.emplace_back(global[nT + d % nF], global[d], 0.5);
		}
		Eigen::SparseMatrix<double, Eigen::RowMajor, std::int64_t> M64(n, n);
		M64.setFromTriplets(triplets.begin(), triplets.end());
		std::vector<std::int32_t> columns(M64.innerIndexPtr(), M64.innerIndexPtr() + M64.nonZeros()); // 32-bit columns
		fhhos4::detail::BasicCSRView<std::int64_t, std::int32_t> M(n, n, M64.outerIndexPtr(), columns.data(), M64.valuePtr());

		// One DoF per cell and per face (degree 0), in the order of the blocks
		std::vector<std::int64_t> cellOffsets(nT + 1), cellDoFs(nT), faceOffsets(nF + 1), faceDoFs(nF);
		for (std::int64_t i = 0; i <= nT; i++)
			cellOffsets[i] = i;
		for (std::int64_t i = 0; i <= nF; i++)
			faceOffsets[i] = i;
		for (std::int64_t i = 0; i < nT; i++)
			cellDoFs[i] = global[i];
		for (std::int64_t i = 0; i < nF; i++)
			faceDoFs[i] = global[nT + i];
		fhhos4::DoFLists<std::int64_t> cells{ nT, cellOffsets.data(), cellDoFs.data() }, faces{ nF, faceOffsets.data(), faceDoFs.data() };

		fhhos4::HybridBlocks blocks = system.Symmetric ? fhhos4::ExtractSymmetricBlocks(M64, cells, faces) : fhhos4::ExtractNonsymmetricBlocks(M64, cells, faces);
		fhhos4::HybridBlocks mixed = fhhos4::detail::ToEigen(system.Symmetric ? fhhos4::detail::ExtractSymmetricBlocks(M, cells, faces) : fhhos4::detail::ExtractNonsymmetricBlocks(M, cells, faces));
		for (const fhhos4::HybridBlocks* b : { &blocks, &mixed })
		{
			if (system.Symmetric)
			{
				HybridSystem::Matrix lowerA_TT = system.A_TT.triangularView<Eigen::Lower>(), lowerA_FF = system.A_FF.triangularView<Eigen::Lower>();
				EXPECT_EQ((b->A_TT - lowerA_TT).norm(), 0);
				EXPECT_EQ((b->A_FF - lowerA_FF).norm(), 0);
				EXPECT_EQ(b->A_FT.rows(), 0);
			}
			else
			{
				EXPECT_EQ((b->A_TT - system.A_TT).norm(), 0);
				EXPECT_EQ((b->A_FT - system.A_FT).norm(), 0);
				EXPECT_EQ((b->A_FF - system.A_FF).norm(), 0);
			}
			EXPECT_EQ((b->A_TF - system.A_TF).norm(), 0);
		}

		// The vectors: b in the global numbering (zero on the Dirichlet DoFs: homogeneous conditions)
		Eigen::VectorXd b_T = Eigen::VectorXd::Random(nT), b_F = Eigen::VectorXd::Random(nF), b = Eigen::VectorXd::Zero(n);
		fhhos4::Scatter(b_T, b_F, cells, faces, b);
		Eigen::VectorXd gathered_T(nT), gathered_F(nF);
		fhhos4::Gather(b, cells, faces, gathered_T, gathered_F);
		EXPECT_TRUE(gathered_T == b_T && gathered_F == b_F);

		fhhos4::Solver fromGlobal = ToySolver(), fromBlocks = ToySolver();
		if (system.Symmetric)
		{
			fromGlobal.SetupSymmetric(blocks.A_TT, blocks.A_TF, blocks.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
			fromBlocks.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
		}
		else
		{
			fromGlobal.SetupNonsymmetric(blocks.A_TT, blocks.A_TF, blocks.A_FT, blocks.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
			fromBlocks.SetupNonsymmetric(system.A_TT, system.A_TF, system.A_FT, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
		}
		Eigen::VectorXd x_T(nT), x_F = Eigen::VectorXd::Zero(nF), y_T(nT), y_F = Eigen::VectorXd::Zero(nF);
		fhhos4::Result r = fromBlocks.Solve(b_T, b_F, x_T, x_F);
		fhhos4::Result s = fromGlobal.Solve(b_T, b_F, y_T, y_F);
		EXPECT_TRUE(r.Converged) << asymmetry;
		EXPECT_EQ(s.Iterations, r.Iterations) << asymmetry;
		EXPECT_TRUE(y_T == x_T && y_F == x_F) << asymmetry;

		// The solution in the global numbering: the Dirichlet DoFs unchanged
		Eigen::VectorXd x = Eigen::VectorXd::Constant(n, 2.0);
		fhhos4::Scatter(y_T, y_F, cells, faces, x);
		for (std::int64_t d = nT + nF; d < n; d++)
			EXPECT_EQ(x[global[d]], 2.0);
		EXPECT_EQ(x[global[0]], y_T[0]);
		EXPECT_EQ(x[global[nT]], y_F[0]);
	}
}

// The inputs not in the CSR format are refused, before the previous setup is replaced
TEST(UAMGLibrary, CSRFormat)
{
	HybridSystem system(16);
	fhhos4::Solver uamg = ToySolver();
	uamg.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
	Eigen::VectorXd x = uamg.Solve(system.b);

	const HybridSystem::Matrix& A_TF = system.A_TF;
	std::vector<int> rowPointers(A_TF.outerIndexPtr(), A_TF.outerIndexPtr() + A_TF.rows() + 1);
	std::vector<int> columns(A_TF.innerIndexPtr(), A_TF.innerIndexPtr() + A_TF.nonZeros());
	auto setupWith = [&](const std::vector<int>& p, const std::vector<int>& c)
	{
		fhhos4::SparseMatrixMap view(A_TF.rows(), A_TF.cols(), A_TF.nonZeros(), p.data(), c.data(), A_TF.valuePtr());
		uamg.SetupSymmetric(system.A_TT, view, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
	};
	EXPECT_NO_THROW(setupWith(rowPointers, columns));
	std::vector<int> unsorted = columns, duplicate = columns, outside = columns, shifted = rowPointers;
	std::swap(unsorted[0], unsorted[1]);
	duplicate[1] = duplicate[0];
	outside[5] = A_TF.cols();
	shifted[0] = 1;
	for (auto [p, c] : { std::pair{ &rowPointers, &unsorted }, { &rowPointers, &duplicate }, { &rowPointers, &outside }, { &shifted, &columns } })
	{
		EXPECT_THROW(setupWith(*p, *c), fhhos4::Error);
		EXPECT_TRUE(uamg.IsSetUp());
		EXPECT_LT((uamg.Solve(system.b) - x).norm(), 1e-12 * x.norm());
	}
}

// Whether SetupNonsymmetric() accepts an argument A_FT of type AType: const Matrix& (a variable), Matrix (a temporary)
template <typename AType, typename = void>
struct SetupAccepts : std::false_type {};
template <typename AType>
struct SetupAccepts<AType, std::void_t<decltype(std::declval<fhhos4::Solver&>().SetupNonsymmetric(std::declval<const HybridSystem::Matrix&>(), std::declval<const HybridSystem::Matrix&>(), std::declval<AType>(), std::declval<const HybridSystem::Matrix&>(), Eigen::VectorXd(), Eigen::VectorXd()))>> : std::true_type {};
static_assert(SetupAccepts<const HybridSystem::Matrix&>::value);
static_assert(!SetupAccepts<HybridSystem::Matrix>::value, "the solver keeps a reference to A_FT: a temporary A_FT must not compile");
static_assert(SetupAccepts<fhhos4::SparseMatrixMap>::value, "a temporary map owns no array: it may be given");

// BiCGSTAB preconditioned by a V-cycle (linear)
TEST(UAMGLibrary, BiCGSTAB)
{
	HybridSystem system(64);
	fhhos4::Solver uamg = ToySolver();
	uamg.Krylov = fhhos4::Krylov::BiCGSTAB;
	uamg.Cycle = fhhos4::Cycle::V;
	uamg.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
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
	uamg.Krylov = fhhos4::Krylov::GCR;
	uamg.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
	Eigen::VectorXd x = Eigen::VectorXd::Zero(system.b.size());
	fhhos4::Result result = uamg.Solve(system.b, x);
	EXPECT_TRUE(result.Converged);
	EXPECT_LT(result.Iterations, 30);
	EXPECT_LT(RelativeError(system, x), 1e-8);
}

// The default smoothers: the Chebyshev of the fourth kind of degree 2 for a symmetric matrix, the same results on 4 and
// 1 threads; for a non-symmetric one, the hybrid Gauss-Seidel (parallel) on several threads, the sequential one on 1
TEST(UAMGLibrary, DefaultSmoothers)
{
	for (bool symmetric : { true, false })
	{
		HybridSystem system = symmetric ? HybridSystem(64) : NonsymmetricSystem(64, 0.3);
		Eigen::VectorXd onFourThreads;
		for (int threads : { 4, 1 })
		{
			fhhos4::Solver uamg = ToySolver();
			uamg.PreSmoother = fhhos4::Smoother::Default;
			uamg.PostSmoother = fhhos4::Smoother::Default;
			uamg.Threads = threads;
			uamg.Verbosity = 2; // prints the smoothers
			::testing::internal::CaptureStdout();
			if (symmetric)
				uamg.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
			else
				uamg.SetupNonsymmetric(system.A_TT, system.A_TF, system.A_FT, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
			Eigen::VectorXd x = Eigen::VectorXd::Zero(system.b.size());
			fhhos4::Result result = uamg.Solve(system.b, x);
			std::string output = ::testing::internal::GetCapturedStdout();
			std::string description = std::string(symmetric ? "symmetric" : "non-symmetric") + ", " + std::to_string(threads) + " threads";
			bool chebyshev = symmetric, hybrid = !symmetric && threads > 1;
			EXPECT_EQ(output.find("PreSmoothing        : 2 iterations of Chebyshev (4th kind") != std::string::npos, chebyshev) << description;
			EXPECT_EQ(output.find("PostSmoothing       : 2 iterations of Chebyshev (4th kind") != std::string::npos, chebyshev) << description;
			EXPECT_EQ(output.find("1 iteration of Hybrid Gauss-Seidel (forward)") != std::string::npos, hybrid) << description;
			EXPECT_EQ(output.find("1 iteration of Hybrid Gauss-Seidel (backward)") != std::string::npos, hybrid) << description;
			EXPECT_EQ(output.find("1 iteration of Gauss-Seidel (backward)") != std::string::npos, !symmetric && threads == 1) << description;
			EXPECT_TRUE(result.Converged) << description;
			EXPECT_LT(result.Iterations, 30) << description;
			EXPECT_LT(RelativeError(system, x), 1e-8) << description;
			if (threads == 4)
				onFourThreads = x;
			else if (symmetric) // Chebyshev: the same results whatever the number of threads
				EXPECT_TRUE((x.array() == onFourThreads.array()).all()) << description;
		}
	}
}

// The Chebyshev smoothers (first and fourth kind), of degree 2 by default, the post-smoother being the pre-smoother by
// default
TEST(UAMGLibrary, ChebyshevSmoothers)
{
	HybridSystem system(64);
	for (auto [smoother, kind] : { std::pair{ fhhos4::Smoother::BlockChebyshevFirstKind, "1st kind" }, { fhhos4::Smoother::BlockChebyshevFourthKind, "4th kind" } })
	{
		fhhos4::Solver uamg = ToySolver();
		uamg.PreSmoother = smoother;
		uamg.PostSmoother = fhhos4::Smoother::Default;
		uamg.Verbosity = 2; // prints the smoothers
		::testing::internal::CaptureStdout();
		uamg.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
		Eigen::VectorXd x = Eigen::VectorXd::Zero(system.b.size());
		fhhos4::Result result = uamg.Solve(system.b, x);
		std::string output = ::testing::internal::GetCapturedStdout();
		EXPECT_NE(output.find("PreSmoothing        : 2 iterations of Chebyshev (" + std::string(kind)), std::string::npos) << kind;
		EXPECT_NE(output.find("PostSmoothing       : 2 iterations of Chebyshev (" + std::string(kind)), std::string::npos) << kind;
		EXPECT_TRUE(result.Converged) << kind;
		EXPECT_LT(result.Iterations, 30) << kind;
		EXPECT_LT(RelativeError(system, x), 1e-8) << kind;
	}
}

// The default post-smoother: the adjoint of the pre-smoother, which makes the cycle symmetric (here at k = 0: the
// scalar kernels)
TEST(UAMGLibrary, DefaultPostSmoother)
{
	HybridSystem system(32);
	std::vector<std::pair<fhhos4::Smoother, std::string>> cases = {
		{ fhhos4::Smoother::BlockGaussSeidel,               "Gauss-Seidel (backward)" },
		{ fhhos4::Smoother::BackwardBlockGaussSeidel,       "Gauss-Seidel (forward)" },
		{ fhhos4::Smoother::SymmetricBlockGaussSeidel,      "Symmetric Gauss-Seidel" },
		{ fhhos4::Smoother::HybridBlockGaussSeidel,         "Hybrid Gauss-Seidel (backward)" },
		{ fhhos4::Smoother::HybridBackwardBlockGaussSeidel, "Hybrid Gauss-Seidel (forward)" },
		{ fhhos4::Smoother::DampedBlockJacobi,              "Jacobi (omega=2/3)" }
	};
	for (const auto& [preSmoother, postSmoother] : cases)
	{
		fhhos4::Solver uamg = ToySolver();
		uamg.PreSmoother = preSmoother;
		uamg.PostSmoother = fhhos4::Smoother::Default;
		uamg.Verbosity = 2; // prints the smoothers
		::testing::internal::CaptureStdout();
		uamg.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
		Eigen::VectorXd x = Eigen::VectorXd::Zero(system.b.size());
		fhhos4::Result result = uamg.Solve(system.b, x);
		std::string output = ::testing::internal::GetCapturedStdout();
		EXPECT_NE(output.find("PostSmoothing       : 1 iteration of " + postSmoother), std::string::npos) << postSmoother;
		EXPECT_TRUE(result.Converged) << postSmoother;
		EXPECT_LT(RelativeError(system, x), 1e-8) << postSmoother;
	}
}

// U-AMG alone (no FCG)
TEST(UAMGLibrary, MultigridAlone)
{
	HybridSystem system(64);
	fhhos4::Solver uamg = ToySolver();
	uamg.Krylov = fhhos4::Krylov::None;
	uamg.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
	Eigen::VectorXd x = uamg.Solve(system.b);
	EXPECT_LT(RelativeError(system, x), 1e-8);
}

// Mixed precision (MixedPrecision, the default: the levels in single precision), with the Chebyshev smoothers (the
// default for a symmetric matrix): the Krylov methods and the cycles alone (each one on the residual of the current x,
// in double) reach the tolerance of the double precision, in the iterations of the double precision (within 1). With
// the options that the mixed precision doesn't support (here the Gauss-Seidel smoothers), the levels are in double:
// the same solution.
TEST(UAMGLibrary, MixedPrecision)
{
	HybridSystem system(64);
	for (fhhos4::Krylov krylov : { fhhos4::Krylov::FCG, fhhos4::Krylov::GCR, fhhos4::Krylov::BiCGSTAB, fhhos4::Krylov::None })
	{
		int iterations[2];
		for (bool mixed : { false, true })
		{
			std::string label = "Krylov " + std::to_string((int)krylov) + (mixed ? ", mixed precision" : ", double precision");
			fhhos4::Solver uamg = ToySolver();
			uamg.PreSmoother = fhhos4::Smoother::Default;
			uamg.PostSmoother = fhhos4::Smoother::Default;
			uamg.Krylov = krylov;
			if (krylov == fhhos4::Krylov::BiCGSTAB)
				uamg.Cycle = fhhos4::Cycle::V;
			uamg.MixedPrecision = mixed;
			uamg.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
			Eigen::VectorXd x = Eigen::VectorXd::Zero(system.b.size());
			fhhos4::Result result = uamg.Solve(system.b, x);
			EXPECT_TRUE(result.Converged) << label;
			EXPECT_LT(result.RelativeResidual, 1e-10) << label;
			EXPECT_LT(RelativeError(system, x), 1e-8) << label;
			iterations[mixed] = result.Iterations;
			// Non-zero initial guess: x is the solution already
			EXPECT_LE(uamg.Solve(system.b, x).Iterations, 1) << label;
		}
		EXPECT_NEAR(iterations[1], iterations[0], 1) << "Krylov " << (int)krylov;
	}

	Eigen::VectorXd solutions[2];
	for (bool mixed : { false, true })
	{
		fhhos4::Solver uamg = ToySolver(); // block Gauss-Seidel
		uamg.MixedPrecision = mixed;
		uamg.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
		solutions[mixed] = uamg.Solve(system.b);
	}
	EXPECT_EQ(solutions[1], solutions[0]);
}

// Non-symmetric system (SetupNonsymmetric()): by default GCR, the K-cycle with GCR, LU on the coarsest level. Also
// BiCGSTAB with a V-cycle.
TEST(UAMGLibrary, NonsymmetricSystem)
{
	for (double asymmetry : { 0.3, 0.9 })
	{
		HybridSystem system = NonsymmetricSystem(64, asymmetry);
		ASSERT_GT((system.A - HybridSystem::Matrix(system.A.transpose())).norm(), 0.1 * system.A.norm());

		fhhos4::Solver gcr = ToySolver(), bicgstab = ToySolver();
		bicgstab.Krylov = fhhos4::Krylov::BiCGSTAB;
		bicgstab.Cycle = fhhos4::Cycle::V;
		gcr.SetupNonsymmetric(system.A_TT, system.A_TF, system.A_FT, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
		bicgstab.SetupNonsymmetric(system.A_TT, system.A_TF, system.A_FT, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
		EXPECT_GT(gcr.NumberOfLevels(), 2);

		Eigen::VectorXd x = Eigen::VectorXd::Zero(system.b.size()), y = x;
		fhhos4::Result r = gcr.Solve(system.b, x);
		fhhos4::Result s = bicgstab.Solve(system.b, y);
		for (const fhhos4::Result& result : { r, s })
		{
			EXPECT_TRUE(result.Converged);
			EXPECT_LT(result.Iterations, 30);
		}
		EXPECT_LT(RelativeError(system, x), 1e-8);
		EXPECT_LT(RelativeError(system, y), 1e-8);
	}
}

// SetupNonsymmetric() on a symmetric system (A_FT = A_TF^T): the coarsening of its symmetric part, the same as
// SetupSymmetric()'s, and Galerkin coarse operators equal to the symmetric ones up to rounding: the same
// iterations with GCR
TEST(UAMGLibrary, NonsymmetricSetupOfASymmetricSystem)
{
	HybridSystem system(64);
	fhhos4::Solver symmetric = ToySolver(), nonsymmetric = ToySolver();
	symmetric.Krylov = fhhos4::Krylov::GCR;
	symmetric.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
	nonsymmetric.SetupNonsymmetric(system.A_TT, system.A_TF, system.A_FT, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
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

	// Vectors of the sizes of the cells and faces DoFs
	Eigen::VectorXd b_T = Eigen::VectorXd::Ones(system.A_TT.rows()), x_T(system.A_TT.rows()), b_F = Eigen::VectorXd::Zero(system.b.size());
	Eigen::VectorXd tooShortT = b_T.head(3), tooShortF = b_F.head(3);
	fhhos4::Solver uncondensed = ToySolver();
	uncondensed.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
	EXPECT_THROW(uncondensed.Solve(tooShortT, b_F, x_T, x), fhhos4::Error);
	EXPECT_THROW(uncondensed.Solve(b_T, b_F, tooShortT, x), fhhos4::Error);
	EXPECT_THROW(uncondensed.Solve(b_T, tooShortF, x_T, x), fhhos4::Error);
	EXPECT_THROW(uncondensed.Solve(b_T, b_F, x_T, tooShortF), fhhos4::Error);

	// Interpolations of 1 of the wrong size, or with a zero coefficient
	Eigen::VectorXd tooShort = system.FaceInterpOfOne.head(10);
	EXPECT_THROW(uamg.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, tooShort), fhhos4::Error);
	Eigen::VectorXd zero = system.CellInterpOfOne;
	zero[3] = 0;
	EXPECT_THROW(uamg.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, zero, system.FaceInterpOfOne), fhhos4::Error);
	EXPECT_FALSE(uamg.IsSetUp());

	// Blocks of inconsistent sizes
	HybridSystem other(8);
	EXPECT_THROW(uamg.SetupSymmetric(other.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);

	// Wrong degrees: 3 DoFs per cell, or 2 per face, in 2D
	fhhos4::Solver cellDegree1 = ToySolver();
	cellDegree1.CellDegree = 1;
	EXPECT_THROW(cellDegree1.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);
	fhhos4::Solver faceDegree1 = ToySolver();
	faceDegree1.FaceDegree = 1;
	EXPECT_THROW(faceDegree1.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);
	// Degree 0 in Q: 1 DoF per cell and per face, as in P
	fhhos4::Solver q = ToySolver();
	q.PolynomialSpace = fhhos4::PolynomialSpace::Q;
	EXPECT_NO_THROW(q.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne));

	// Parameters not set, or invalid
	fhhos4::Solver noDimension = ToySolver();
	noDimension.Dimension = -1;
	EXPECT_THROW(noDimension.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);
	fhhos4::Solver noFaceDegree = ToySolver();
	noFaceDegree.FaceDegree = -1;
	EXPECT_THROW(noFaceDegree.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);
	fhhos4::Solver noCellDegree = ToySolver();
	noCellDegree.CellDegree = -1;
	EXPECT_THROW(noCellDegree.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);
	fhhos4::Solver wrongCycle = ToySolver();
	wrongCycle.Cycle = static_cast<fhhos4::Cycle>(42); // out of the enum
	EXPECT_THROW(wrongCycle.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);
	fhhos4::Solver wrongPCoarsening = ToySolver();
	wrongPCoarsening.PCoarsening = static_cast<fhhos4::PCoarsening>(42);
	EXPECT_THROW(wrongPCoarsening.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);
	fhhos4::Solver wrongProlongation = ToySolver();
	wrongProlongation.Prolongation = static_cast<fhhos4::Prolongation>(42);
	EXPECT_THROW(wrongProlongation.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);
	fhhos4::Solver wrongPolynomialSpace = ToySolver();
	wrongPolynomialSpace.PolynomialSpace = static_cast<fhhos4::PolynomialSpace>(42);
	EXPECT_THROW(wrongPolynomialSpace.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);
	fhhos4::Solver negativeSmoothing = ToySolver();
	negativeSmoothing.PostSmoothingIterations = -2;
	EXPECT_THROW(negativeSmoothing.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);
	fhhos4::Solver wrongKrylov = ToySolver();
	wrongKrylov.Krylov = static_cast<fhhos4::Krylov>(42);
	EXPECT_THROW(wrongKrylov.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);
	fhhos4::Solver bicgstabKCycle = ToySolver(); // BiCGSTAB is not flexible: the K-cycle is not linear
	bicgstabKCycle.Krylov = fhhos4::Krylov::BiCGSTAB;
	EXPECT_THROW(bicgstabKCycle.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);

	// Non-symmetric matrix: FCG or Cholesky refused, A_FT of the wrong size
	fhhos4::Solver fcg = ToySolver();
	fcg.Krylov = fhhos4::Krylov::FCG;
	EXPECT_THROW(fcg.SetupNonsymmetric(system.A_TT, system.A_TF, system.A_FT, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);
	fhhos4::Solver cholesky = ToySolver();
	cholesky.CoarseSolver = fhhos4::CoarseSolver::Cholesky;
	EXPECT_THROW(cholesky.SetupNonsymmetric(system.A_TT, system.A_TF, system.A_FT, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);
	fhhos4::Solver chebyshev = ToySolver();
	chebyshev.PreSmoother = fhhos4::Smoother::BlockChebyshevFourthKind;
	chebyshev.PostSmoother = fhhos4::Smoother::Default;
	EXPECT_THROW(chebyshev.SetupNonsymmetric(system.A_TT, system.A_TF, system.A_FT, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);
	EXPECT_THROW(uamg.SetupNonsymmetric(system.A_TT, system.A_TF, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);

	// A_TT not block diagonal
	HybridSystem::Matrix A_TT = system.A_TT;
	A_TT.coeffRef(1, 0) = -1;
	A_TT.makeCompressed();
	EXPECT_THROW(uamg.SetupSymmetric(A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);

	// Vector of the wrong size
	uamg.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
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
	uamg.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
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
	verbose.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
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
