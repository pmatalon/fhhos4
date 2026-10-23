// Unit tests of the Chebyshev smoothers (BlockChebyshev, first and fourth kinds, scalar and block): the error after m
// iterations is p(D^-1 A) times the initial error, p given by its closed form, computed here from the eigenvalues of
// D^-1 A; the results don't depend on the number of threads.
#include <gtest/gtest.h>
#include <random>
#include <string>
#include <vector>
#include <Eigen/Eigenvalues>
#include "Solver/FixedPoint/BlockChebyshev.h"
#include "Utils/SparseMatrixOps.h"

namespace
{
	// Symmetric, strictly diagonally dominant matrix (SPD) of a 9-point stencil on an m x m grid, in blocks of blockSize
	// DoFs (dense blocks): its diagonal is the sum of the absolute values of the row plus 0.01, so that the spectrum of
	// D^-1 A is spread on (0, 2), as a Laplacian's
	SparseMatrix SPDMatrix(int blockSize, int m)
	{
		std::mt19937 gen(42);
		std::uniform_real_distribution<double> random(-1, 0);
		std::vector<Eigen::Triplet<double, SparseMatrixIndex>> triplets;
		int nb = m * m;
		Vector rowSums = Vector::Zero(nb * blockSize);
		auto add = [&](BigNumber i, BigNumber j, double value)
		{
			triplets.push_back({ (SparseMatrixIndex)i, (SparseMatrixIndex)j, value });
			triplets.push_back({ (SparseMatrixIndex)j, (SparseMatrixIndex)i, value });
			rowSums[i] += abs(value);
			rowSums[j] += abs(value);
		};
		for (int y = 0; y < m; y++)
		{
			for (int x = 0; x < m; x++)
			{
				int I = y * m + x;
				for (int k = 0; k < blockSize; k++) // diagonal block
				{
					for (int l = k + 1; l < blockSize; l++)
						add(I * blockSize + k, I * blockSize + l, 0.5 * random(gen));
				}
				for (int dy = 0; dy <= 1; dy++)
				{
					for (int dx = -1; dx <= 1; dx++)
					{
						int J = (y + dy) * m + x + dx;
						if ((dy == 0 && dx <= 0) || y + dy >= m || x + dx < 0 || x + dx >= m)
							continue;
						for (int k = 0; k < blockSize; k++)
						{
							for (int l = 0; l < blockSize; l++)
								add(I * blockSize + k, J * blockSize + l, random(gen));
						}
					}
				}
			}
		}
		for (BigNumber i = 0; i < nb * blockSize; i++)
			triplets.push_back({ (SparseMatrixIndex)i, (SparseMatrixIndex)i, rowSums[i] + 0.01 });
		SparseMatrix A(nb * blockSize, nb * blockSize);
		A.setFromTriplets(triplets.begin(), triplets.end());
		return A;
	}

	std::unique_ptr<BlockChebyshev> CreateChebyshev(const SparseMatrix& A, int blockSize, ChebyshevKind kind, int degree)
	{
		auto chebyshev = std::make_unique<BlockChebyshev>(blockSize, kind, 10.0);
		chebyshev->StoppingCrit = StoppingCriteria::MaxIterations;
		chebyshev->MaxIterations = degree;
		chebyshev->PrintIterationResults = false;
		chebyshev->Setup(A);
		return chebyshev;
	}

	// The error polynomial of the kind at lambda: T_m((theta - lambda) / delta) / T_m(theta / delta) on [lambdaMax/ratio,
	// lambdaMax] (first kind), W_m(1 - 2 lambda / lambdaMax) / (2m + 1) (fourth kind)
	double ErrorPolynomial(ChebyshevKind kind, int m, double lambdaMax, double ratio, double lambda)
	{
		auto chebyshev = [m](double t, double p1) // T_m(t) (p1 = t) or W_m(t) (p1 = 2t + 1)
		{
			double p0 = 1;
			for (int k = 1; k < m; k++)
			{
				double p2 = 2 * t * p1 - p0;
				p0 = p1;
				p1 = p2;
			}
			return m == 0 ? p0 : p1;
		};
		if (kind == ChebyshevKind::First)
		{
			double lambdaMin = lambdaMax / ratio;
			double theta = (lambdaMax + lambdaMin) / 2, delta = (lambdaMax - lambdaMin) / 2;
			double t = (theta - lambda) / delta, sigma = theta / delta;
			return chebyshev(t, t) / chebyshev(sigma, sigma);
		}
		double t = 1 - 2 * lambda / lambdaMax;
		return chebyshev(t, 2 * t + 1) / (2 * m + 1);
	}

	// D: the block diagonal of A
	DenseMatrix BlockDiagonal(const DenseMatrix& A, int blockSize)
	{
		DenseMatrix D = DenseMatrix::Zero(A.rows(), A.cols());
		for (BigNumber I = 0; I < A.rows() / blockSize; I++)
			D.block(I * blockSize, I * blockSize, blockSize, blockSize) = A.block(I * blockSize, I * blockSize, blockSize, blockSize);
		return D;
	}

	double MaxDiff(const Vector& x, const Vector& y)
	{
		return (x - y).lpNorm<Eigen::Infinity>();
	}

	// Exactly equal (+0 and -0 being equal)
	bool Identical(const Vector& x, const Vector& y)
	{
		return x.size() == y.size() && (x.array() == y.array()).all();
	}

	// Number of threads of the solvers while in scope
	class Threads
	{
	public:
		Threads(int n) { Parallelism::SetNThreads(n); }
		~Threads() { Parallelism::SetNThreads(0); }
	};
}

// x_m = x* - p(D^-1 A) (x* - x_0), with the generalized eigenvectors of (A, D): D^-1 A = V Lambda V^T D. Also the
// residual and Ax given by the smoother, and the estimate of lambdaMax (from below, close).
TEST(ChebyshevSmootherTest, ErrorPolynomial)
{
	for (int blockSize : { 1, 3 })
	{
		SparseMatrix A = SPDMatrix(blockSize, 6);
		DenseMatrix denseA = A;
		Eigen::GeneralizedSelfAdjointEigenSolver<DenseMatrix> eigen(denseA, BlockDiagonal(denseA, blockSize));
		const Vector& lambdas = eigen.eigenvalues();
		const DenseMatrix& V = eigen.eigenvectors();
		DenseMatrix VtD = V.transpose() * BlockDiagonal(denseA, blockSize);
		Vector b = Vector::LinSpaced(A.rows(), -1, 1);
		Vector exact = denseA.ldlt().solve(b);

		for (ChebyshevKind kind : { ChebyshevKind::First, ChebyshevKind::Fourth })
		{
			for (int degree : { 1, 2, 3, 5 })
			{
				auto chebyshev = CreateChebyshev(A, blockSize, kind, degree);
				std::string description = "block size " + std::to_string(blockSize) + (kind == ChebyshevKind::First ? ", 1st kind" : ", 4th kind") + ", degree " + std::to_string(degree);
				double estimate = chebyshev->LambdaMaxEstimate();
				EXPECT_LE(estimate, lambdas.maxCoeff() * (1 + 1e-12)) << description;
				EXPECT_GE(estimate, 0.95 * lambdas.maxCoeff()) << description;
				double lambdaMax = BlockChebyshev::EigenvalueSafetyFactor * estimate;

				for (bool xEquals0 : { true, false })
				{
					Vector x0 = xEquals0 ? Vector::Zero(A.rows()) : Vector(Vector::LinSpaced(A.rows(), 2, -3).array().sin().matrix());
					Vector p(lambdas.size());
					for (BigNumber i = 0; i < lambdas.size(); i++)
						p[i] = ErrorPolynomial(kind, degree, lambdaMax, 10.0, lambdas[i]);
					Vector expected = exact - V * p.asDiagonal() * VtD * (exact - x0);

					for (int outputs = 0; outputs < 4; outputs++) // x alone, the residual, Ax, both
					{
						bool computeResidual = outputs & 1, computeAx = outputs & 2;
						Vector x = x0;
						chebyshev->Solve(b, x, xEquals0, computeResidual, computeAx);
						std::string what = description + (xEquals0 ? ", from 0" : "") + ", outputs " + std::to_string(outputs);
						EXPECT_LT(MaxDiff(x, expected), 1e-11 * expected.lpNorm<Eigen::Infinity>()) << what;
						if (computeResidual)
							EXPECT_LT(MaxDiff(chebyshev->Residual, b - A * x), 1e-12 * b.lpNorm<Eigen::Infinity>()) << what;
						if (computeAx)
							EXPECT_LT(MaxDiff(chebyshev->Ax, A * x), 1e-12 * b.lpNorm<Eigen::Infinity>()) << what;
					}
				}
			}
		}
	}
}

// The same results on 1 and 4 threads (setup and smoothing), on matrices large enough for the parallel loops
TEST(ChebyshevSmootherTest, IndependentOfTheNumberOfThreads)
{
	for (int blockSize : { 1, 3 })
	{
		SparseMatrix A = SPDMatrix(blockSize, 60);
		ASSERT_GT(A.nonZeros(), 20000);
		Vector b = Vector::LinSpaced(A.rows(), -1, 1);
		Vector x0 = Vector::LinSpaced(A.rows(), 2, -3).array().sin().matrix();
		for (ChebyshevKind kind : { ChebyshevKind::First, ChebyshevKind::Fourth })
		{
			Vector x[2], r[2];
			double estimate[2];
			for (int t = 0; t < 2; t++)
			{
				Threads threads(t == 0 ? 1 : 4);
				auto chebyshev = CreateChebyshev(A, blockSize, kind, 3);
				estimate[t] = chebyshev->LambdaMaxEstimate();
				x[t] = x0;
				chebyshev->Solve(b, x[t], false, true, false);
				r[t] = chebyshev->Residual;
			}
			std::string description = "block size " + std::to_string(blockSize) + (kind == ChebyshevKind::First ? ", 1st kind" : ", 4th kind");
			EXPECT_EQ(estimate[0], estimate[1]) << description;
			EXPECT_TRUE(Identical(x[0], x[1])) << description;
			EXPECT_TRUE(Identical(r[0], r[1])) << description;
		}
	}
}

// A matrix that is not positive definite is refused in the setup
TEST(ChebyshevSmootherTest, NotPositiveDefinite)
{
	SparseMatrix A = SPDMatrix(1, 6);
	SparseMatrix negative = -A;
	auto chebyshev = std::make_unique<BlockChebyshev>(1, ChebyshevKind::Fourth);
	chebyshev->MaxIterations = 2;
	EXPECT_THROW(chebyshev->Setup(negative), fhhos4::Error);
	SparseMatrix negativeBlocks = -SPDMatrix(3, 6);
	auto block = std::make_unique<BlockChebyshev>(3, ChebyshevKind::Fourth);
	block->MaxIterations = 2;
	EXPECT_THROW(block->Setup(negativeBlocks), fhhos4::Error);
}

// The matrix stored by blocks (block CSR, SetupBlockCSR()): the same lambdaMax, iterate, residual and A x as in CSR, to
// rounding (block sizes compiled for: 2, 3, 6; at run time: 7). The passes sum in the same order, but Lanczos' products
// differ: Eigen 5's sparse product sums with two accumulators, and in single precision in float.
TEST(ChebyshevSmootherTest, BlockCSR)
{
	for (int blockSize : { 2, 3, 6, 7 })
	{
		SparseMatrix A = SPDMatrix(blockSize, 30);
		Vector b = Vector::LinSpaced(A.rows(), -1, 1);
		Vector x0 = Vector::LinSpaced(A.rows(), 2, -3).array().sin().matrix();
		std::string description = "block size " + std::to_string(blockSize);

		auto csr = CreateChebyshev(A, blockSize, ChebyshevKind::Fourth, 3);
		BlockSparseMatrix<double> blocks = BlockSparseMatrix<double>::FromCSR(A, blockSize);
		auto bcsr = std::make_unique<BlockChebyshev>(blockSize, ChebyshevKind::Fourth);
		bcsr->StoppingCrit = StoppingCriteria::MaxIterations;
		bcsr->MaxIterations = 3;
		bcsr->PrintIterationResults = false;
		ASSERT_TRUE(bcsr->SupportsBlockCSR());
		bcsr->SetupBlockCSR(blocks);
		EXPECT_NEAR(bcsr->LambdaMaxEstimate(), csr->LambdaMaxEstimate(), 1e-13 * csr->LambdaMaxEstimate()) << description;
		Vector x[2] = { x0, x0 };
		csr->Solve(b, x[0], false, true, true);
		bcsr->Solve(b, x[1], false, true, true);
		EXPECT_LT((x[1] - x[0]).norm() / x[0].norm(), 1e-12) << description;
		EXPECT_LT((bcsr->Residual - csr->Residual).norm() / csr->Residual.norm(), 1e-12) << description;
		EXPECT_LT((bcsr->Ax - csr->Ax).norm() / csr->Ax.norm(), 1e-12) << description;

		// Single precision
		FloatSparseMatrix singleA = SparseMatrixOps::Cast<float>(A);
		BlockSparseMatrix<float> singleBlocks = BlockSparseMatrix<float>::FromCSR(A, blockSize);
		csr->SetupSinglePrecision(singleA);
		bcsr->SetupBlockCSR(singleBlocks);
		FloatVector singleB = b.cast<float>(), y[2] = { x0.cast<float>(), x0.cast<float>() }, residual[2];
		csr->SolveInSinglePrecision(singleB, y[0], false, &residual[0], nullptr);
		bcsr->SolveInSinglePrecision(singleB, y[1], false, &residual[1], nullptr);
		EXPECT_LT((y[1] - y[0]).norm() / y[0].norm(), 1e-5) << description;
		EXPECT_LT((residual[1] - residual[0]).norm() / residual[0].norm(), 1e-5) << description;
	}
}
