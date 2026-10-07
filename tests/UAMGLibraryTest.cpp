// Tests of the library fhhos4_AMG (library/), U-AMG for the codes that discretize the problem themselves.
// - UAMGLibraryProgram: the program runs U-AMG through the library (-s libuamg, fcglibuamg), as an external code calls
//   it, and must get the iteration counts and solutions of the U-AMG compiled in the program (uamg, fcguamg), including
//   with L2-orthonormal bases, whose constant does not have coordinate 1.
// - UAMGLibrary: the public API alone (<fhhos4>), on a hybrid system built here.
#include <gtest/gtest.h>
#include <cmath>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <Eigen/SparseCholesky>
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

class UAMGLibraryProgramTest : public ::testing::TestWithParam<std::tuple<LibraryCase, std::string>>
{
};

// The library code (libuamg, fcglibuamg) against the program's (uamg, fcguamg)
TEST_P(UAMGLibraryProgramTest, SameAsProgram)
{
	auto [c, libraryCode] = GetParam();
	std::string programCode = libraryCode;
	programCode.replace(programCode.find("libuamg"), 7, "uamg");

	auto run = [&](const std::string& solverCode)
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
		args.Actions.UseCache = false;
		return RunDiffusionHHO(args);
	};

	ProgramResults program = run(programCode);
	ProgramResults library = run(libraryCode);
	EXPECT_GT(program.IterationCount, 0);
	EXPECT_LT(program.IterationCount, 200); // converged (default maximum: 200)
	EXPECT_EQ(library.IterationCount, program.IterationCount) << libraryCode << " vs " << programCode;
	EXPECT_NEAR(library.L2Error, program.L2Error, 1e-6 * program.L2Error) << libraryCode << " vs " << programCode;
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
	::testing::Combine(::testing::ValuesIn(LibraryCases()), ::testing::Values("fcglibuamg", "libuamg")),
	[](const ::testing::TestParamInfo<UAMGLibraryProgramTest::ParamType>& info)
	{
		return std::get<0>(info.param).Name + "_" + std::get<1>(info.param);
	});

//------------------------------------------------------------------------------------------------------------------//
//                                     Through the public API only (<fhhos4>)                                    //
//------------------------------------------------------------------------------------------------------------------//

// Hybrid system of degree 0 on the n x n Cartesian mesh of the unit square: in each cell T, of coefficient a_T, the
// local matrix couples the cell DoF to each of its 4 faces with the weight a_T (the constants are in its kernel).
// The DoFs of the boundary faces (Dirichlet) are eliminated. Column-major matrices with int indices, as Eigen's
// default (and HArDCore3D's).
struct HybridSystem
{
	using Matrix = Eigen::SparseMatrix<double>;
	Matrix A_TT, A_TF, A_FF, A; // A: the condensed matrix
	Eigen::VectorXd CellInterpOfOne, FaceInterpOfOne, b;

	// cellScaling, faceScaling: the cell (resp. face) basis function is multiplied by these factors
	HybridSystem(int n, const Eigen::VectorXd& cellScaling = Eigen::VectorXd(), const Eigen::VectorXd& faceScaling = Eigen::VectorXd())
	{
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

		std::vector<Eigen::Triplet<double>> TT, TF, FF;
		for (int j = 0; j < n; j++)
		{
			for (int i = 0; i < n; i++)
			{
				int T = j * n + i;
				double a = coefficient(generator);
				TT.emplace_back(T, T, 4 * a * sT[T] * sT[T]);
				std::vector<int> faces;
				if (i > 0)     faces.push_back(vertical(i - 1, j));
				if (i < n - 1) faces.push_back(vertical(i, j));
				if (j > 0)     faces.push_back(horizontal(i, j - 1));
				if (j < n - 1) faces.push_back(horizontal(i, j));
				for (int F : faces)
				{
					TF.emplace_back(T, F, -a * sT[T] * sF[F]);
					FF.emplace_back(F, F, a * sF[F] * sF[F]);
				}
			}
		}
		A_TT.resize(nCells, nCells);
		A_TT.setFromTriplets(TT.begin(), TT.end());
		A_TF.resize(nCells, nFaces);
		A_TF.setFromTriplets(TF.begin(), TF.end());
		A_FF.resize(nFaces, nFaces);
		A_FF.setFromTriplets(FF.begin(), FF.end());

		Eigen::VectorXd invDiag = A_TT.diagonal().cwiseInverse();
		Matrix A_FT = A_TF.transpose();
		A = A_FF - Matrix(A_FT * invDiag.asDiagonal() * A_TF);
		A.makeCompressed();

		CellInterpOfOne = sT.cwiseInverse();
		FaceInterpOfOne = sF.cwiseInverse();
		// Right-hand side of the condensed system for a unit source in the cells: b = -A_FT A_TT^{-1} f
		Eigen::VectorXd f = sT;
		b = -A_FT * invDiag.cwiseProduct(f);
	}
};

// The solver for HybridSystem: degree 0 in 2D
fhhos4::Solver ToySolver()
{
	fhhos4::Solver solver;
	solver.Dimension = 2;
	solver.FaceDegree = 0;
	solver.Tolerance = 1e-10;
	return solver;
}

// Relative error of x against the direct solution
double RelativeError(const HybridSystem& system, const Eigen::VectorXd& x)
{
	Eigen::SimplicialLDLT<HybridSystem::Matrix> direct(system.A);
	Eigen::VectorXd exact = direct.solve(system.b);
	return (x - exact).norm() / exact.norm();
}

TEST(UAMGLibrary, SolvesTheSystem)
{
	HybridSystem system(64);
	fhhos4::Solver uamg = ToySolver();
	uamg.Setup(system.A, system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
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
	withA.Setup(system.A, system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
	withoutA.SetupFromBlocks(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);

	Eigen::VectorXd x1 = Eigen::VectorXd::Zero(system.b.size()), x2 = x1;
	fhhos4::Result r1 = withA.Solve(system.b, x1);
	fhhos4::Result r2 = withoutA.Solve(system.b, x2);
	EXPECT_EQ(withoutA.NumberOfLevels(), withA.NumberOfLevels());
	EXPECT_EQ(r2.Iterations, r1.Iterations);
	EXPECT_LT((x2 - x1).norm() / x1.norm(), 1e-9);

	// Lower triangular parts only of the symmetric blocks
	HybridSystem::Matrix A_TT = system.A_TT.triangularView<Eigen::Lower>(), A_FF = system.A_FF.triangularView<Eigen::Lower>();
	fhhos4::Solver lower = ToySolver();
	lower.SetupFromBlocks(A_TT, system.A_TF, A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
	Eigen::VectorXd x3 = Eigen::VectorXd::Zero(system.b.size());
	EXPECT_EQ(lower.Solve(system.b, x3).Iterations, r1.Iterations);
	EXPECT_LT((x3 - x1).norm() / x1.norm(), 1e-9);
}

// Without A_FF, not used by the algorithm: same solver
TEST(UAMGLibrary, WithoutTheFaceBlock)
{
	HybridSystem system(64);
	fhhos4::Solver withA_FF = ToySolver(), withoutA_FF = ToySolver();
	withA_FF.Setup(system.A, system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
	withoutA_FF.Setup(system.A, system.A_TT, system.A_TF, system.CellInterpOfOne, system.FaceInterpOfOne);

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
	uamgReference.Setup(reference.A, reference.A_TT, reference.A_TF, reference.A_FF, reference.CellInterpOfOne, reference.FaceInterpOfOne);
	uamgScaled.Setup(scaled.A, scaled.A_TT, scaled.A_TF, scaled.A_FF, scaled.CellInterpOfOne, scaled.FaceInterpOfOne);
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

// BiCGSTAB preconditioned by a V-cycle (linear)
TEST(UAMGLibrary, BiCGSTAB)
{
	HybridSystem system(64);
	fhhos4::Solver uamg = ToySolver();
	uamg.Krylov = "bicgstab";
	uamg.Cycle = 'V';
	uamg.Setup(system.A, system.A_TT, system.A_TF, system.CellInterpOfOne, system.FaceInterpOfOne);
	Eigen::VectorXd x = Eigen::VectorXd::Zero(system.b.size());
	fhhos4::Result result = uamg.Solve(system.b, x);
	EXPECT_TRUE(result.Converged);
	EXPECT_LT(result.Iterations, 30);
	EXPECT_LT(RelativeError(system, x), 1e-8);
}

// U-AMG alone (no FCG)
TEST(UAMGLibrary, MultigridAlone)
{
	HybridSystem system(64);
	fhhos4::Solver uamg = ToySolver();
	uamg.Krylov = "none";
	uamg.Setup(system.A, system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
	Eigen::VectorXd x = uamg.Solve(system.b);
	EXPECT_LT(RelativeError(system, x), 1e-8);
}

TEST(UAMGLibrary, Errors)
{
	HybridSystem system(16);
	fhhos4::Solver uamg = ToySolver();
	Eigen::VectorXd x = Eigen::VectorXd::Zero(system.b.size());
	EXPECT_THROW(uamg.Solve(system.b, x), fhhos4::Error); // before Setup()

	// Interpolations of 1 of the wrong size, or with a zero coefficient
	Eigen::VectorXd tooShort = system.FaceInterpOfOne.head(10);
	EXPECT_THROW(uamg.Setup(system.A, system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, tooShort), fhhos4::Error);
	Eigen::VectorXd zero = system.CellInterpOfOne;
	zero[3] = 0;
	EXPECT_THROW(uamg.Setup(system.A, system.A_TT, system.A_TF, system.A_FF, zero, system.FaceInterpOfOne), fhhos4::Error);
	EXPECT_FALSE(uamg.IsSetUp());

	// Blocks of inconsistent sizes
	HybridSystem other(8);
	EXPECT_THROW(uamg.Setup(system.A, other.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);

	// Wrong degrees: 3 DoFs per cell, or 2 per face, in 2D
	fhhos4::Solver cellDegree1 = ToySolver();
	cellDegree1.CellDegree = 1;
	EXPECT_THROW(cellDegree1.Setup(system.A, system.A_TT, system.A_TF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);
	fhhos4::Solver faceDegree1 = ToySolver();
	faceDegree1.FaceDegree = 1;
	EXPECT_THROW(faceDegree1.Setup(system.A, system.A_TT, system.A_TF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);

	// Parameters not set, or invalid
	fhhos4::Solver noDimension = ToySolver();
	noDimension.Dimension = -1;
	EXPECT_THROW(noDimension.Setup(system.A, system.A_TT, system.A_TF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);
	fhhos4::Solver noFaceDegree = ToySolver();
	noFaceDegree.FaceDegree = -1;
	EXPECT_THROW(noFaceDegree.Setup(system.A, system.A_TT, system.A_TF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);
	fhhos4::Solver wrongCycle = ToySolver();
	wrongCycle.Cycle = 'X';
	EXPECT_THROW(wrongCycle.Setup(system.A, system.A_TT, system.A_TF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);
	fhhos4::Solver wrongKrylov = ToySolver();
	wrongKrylov.Krylov = "gmres";
	EXPECT_THROW(wrongKrylov.Setup(system.A, system.A_TT, system.A_TF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);
	fhhos4::Solver bicgstabKCycle = ToySolver(); // BiCGSTAB is not flexible: the K-cycle is not linear
	bicgstabKCycle.Krylov = "bicgstab";
	EXPECT_THROW(bicgstabKCycle.Setup(system.A, system.A_TT, system.A_TF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);

	// A_TT not block diagonal
	HybridSystem::Matrix A_TT = system.A_TT;
	A_TT.coeffRef(1, 0) = -1;
	EXPECT_THROW(uamg.Setup(system.A, A_TT, system.A_TF, system.CellInterpOfOne, system.FaceInterpOfOne), fhhos4::Error);

	// Vector of the wrong size
	uamg.Setup(system.A, system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
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
	uamg.Setup(system.A, system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
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
	verbose.Setup(system.A, system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne);
	verbose.Solve(system.b);
	output = ::testing::internal::GetCapturedStdout();
	EXPECT_NE(output.find("Setup"), std::string::npos);
	EXPECT_EQ(std::cout.precision(), 7);
	EXPECT_EQ(std::cout.flags(), flags);
	std::cout.precision(precision);
}
