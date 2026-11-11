#include <gtest/gtest.h>
#include "Solver/FixedPoint/BlockChebyshev.h"
#include "support/CartesianHybridSystem.h"
#include "support/Kronecker.h"
#include "support/PartitionedSystem.h"

using namespace fhhos4_tests;

namespace
{
	std::unique_ptr<BlockChebyshev> Smoother(ChebyshevKind kind, int degree, int blockSize = 1)
	{
		auto chebyshev = std::make_unique<BlockChebyshev>(blockSize, kind);
		chebyshev->StoppingCrit = StoppingCriteria::MaxIterations;
		chebyshev->MaxIterations = degree;
		chebyshev->PrintIterationResults = false;
		return chebyshev;
	}
}

// Mixed precision: the Chebyshev smoother on the distributed matrix in single precision, its estimate of lambdaMax and its
// smoothing those in double to the rounding of float; the post-smoother reuses the setup of the pre-smoother
TEST(DistributedChebyshev, SinglePrecision)
{
	int n = 30;
	CartesianHybridSystem global(n, 0);
	PartitionedSystem system(global, n, PartitionShape::Blocks, Communicator::Current());
	DistributedMatrixOf<float> singleS = system.S.Cast<float>();
	BigNumber nRows = system.S.LocalRows();
	Vector b = Vector::LinSpaced(nRows, -1, 1), x0 = Vector::LinSpaced(nRows, 2, -3).array().sin().matrix();
	for (ChebyshevKind kind : { ChebyshevKind::First, ChebyshevKind::Fourth })
	{
		std::unique_ptr<BlockChebyshev> d = Smoother(kind, 3), s = Smoother(kind, 3), r = Smoother(kind, 3);
		BlockChebyshev &inDouble = *d, &inSingle = *s, &reusing = *r;
		inDouble.Setup(system.S);
		inSingle.SetupSinglePrecision(singleS);
		ASSERT_TRUE(reusing.SetupSinglePrecisionFrom(singleS, inSingle));
		EXPECT_FALSE(reusing.SetupSinglePrecisionFrom(singleS, inDouble));
		EXPECT_NEAR(inSingle.LambdaMaxEstimate(), inDouble.LambdaMaxEstimate(), 1e-6 * inDouble.LambdaMaxEstimate());
		EXPECT_EQ(reusing.LambdaMaxEstimate(), inSingle.LambdaMaxEstimate());

		for (bool xEquals0 : { true, false })
		{
			Vector x = xEquals0 ? Vector::Zero(nRows) : x0;
			inDouble.Solve(b, x, xEquals0, true, false);
			FloatVector singleX = xEquals0 ? FloatVector::Zero(nRows) : FloatVector(x0.cast<float>()), residual;
			reusing.SolveInSinglePrecision(b.cast<float>(), singleX, xEquals0, &residual, nullptr);
			EXPECT_LE(Parallel::Norm(Vector(singleX.cast<double>() - x)), 1e-5 * Parallel::Norm(x)) << (int)kind << xEquals0;
			EXPECT_LE(Parallel::Norm(Vector(residual.cast<double>() - inDouble.Residual)), 1e-5 * Parallel::Norm(b)) << (int)kind << xEquals0;
		}
	}
}

// Block CSR on several processes: Diag stored by blocks (S ⊗ C, C SPD: dense blocks of 2), Offd in CSR; the same lambdaMax
// (Lanczos' products sum Diag in another order), iterate, residual and A x as in CSR, in double and in single precision;
// the post-smoother reuses the setup of the pre-smoother of the same distributed matrix only
TEST(DistributedChebyshev, BlockCSR)
{
	int n = 30, bs = 2;
	CartesianHybridSystem global(n, 0);
	PartitionedSystem system(global, n, PartitionShape::Blocks, Communicator::Current());
	Eigen::MatrixXd C(bs, bs);
	C << 2, 1, 1, 3;
	DistributedMatrix S = system.S;
	S.RowBlockSize = S.ColBlockSize = bs;
	S.Diag = Kron(system.S.Diag, C);
	S.Offd = Kron(system.S.Offd, C);
	DistributedMatrixOf<float> singleS = S.Cast<float>(false);
	BlockSparseMatrix<double> blocks = BlockSparseMatrix<double>::FromCSR(S.Diag, bs);
	BlockSparseMatrix<float> singleBlocks = BlockSparseMatrix<float>::FromCSR(S.Diag, bs);
	BigNumber nRows = S.LocalRows();
	Vector b = Vector::LinSpaced(nRows, -1, 1), x0 = Vector::LinSpaced(nRows, 2, -3).array().sin().matrix();

	std::unique_ptr<BlockChebyshev> c = Smoother(ChebyshevKind::Fourth, 3, bs), d = Smoother(ChebyshevKind::Fourth, 3, bs), r = Smoother(ChebyshevKind::Fourth, 3, bs);
	BlockChebyshev &csr = *c, &bcsr = *d, &reusing = *r;
	ASSERT_TRUE(bcsr.SupportsBlockCSR());
	csr.Setup(S);
	bcsr.SetupBlockCSR(blocks, &S);
	EXPECT_FALSE(reusing.SetupBlockCSRFrom(blocks, nullptr, bcsr));
	ASSERT_TRUE(reusing.SetupBlockCSRFrom(blocks, &S, bcsr));
	EXPECT_NEAR(bcsr.LambdaMaxEstimate(), csr.LambdaMaxEstimate(), 1e-13 * csr.LambdaMaxEstimate());
	Vector x[2] = { x0, x0 };
	csr.Solve(b, x[0], false, true, true);
	reusing.Solve(b, x[1], false, true, true);
	EXPECT_LE(Parallel::Norm(Vector(x[1] - x[0])), 1e-12 * Parallel::Norm(x[0]));
	EXPECT_LE(Parallel::Norm(Vector(reusing.Residual - csr.Residual)), 1e-12 * Parallel::Norm(csr.Residual));
	EXPECT_LE(Parallel::Norm(Vector(reusing.Ax - csr.Ax)), 1e-12 * Parallel::Norm(csr.Ax));
	// A x by the blocks
	EXPECT_LE(Parallel::Norm(Vector(S.Multiply(x0, &blocks) - S * x0)), 1e-14 * Parallel::Norm(Vector(S * x0)));

	// Single precision: Diag only in block CSR
	DistributedMatrixOf<float> singleCSR = S.Cast<float>();
	csr.SetupSinglePrecision(singleCSR);
	bcsr.SetupBlockCSR(singleBlocks, &singleS);
	ASSERT_TRUE(reusing.SetupBlockCSRFrom(singleBlocks, &singleS, bcsr));
	FloatVector singleB = b.cast<float>(), y[2] = { x0.cast<float>(), x0.cast<float>() }, residual[2];
	csr.SolveInSinglePrecision(singleB, y[0], false, &residual[0], nullptr);
	reusing.SolveInSinglePrecision(singleB, y[1], false, &residual[1], nullptr);
	EXPECT_LE(Parallel::Norm(Vector((y[1] - y[0]).cast<double>())), 1e-5 * Parallel::Norm(Vector(y[0].cast<double>())));
	EXPECT_LE(Parallel::Norm(Vector((residual[1] - residual[0]).cast<double>())), 1e-5 * Parallel::Norm(Vector(residual[0].cast<double>())));

	// Blocks of another size on the last process only: the error raised on all the processes
	const Communicator& comm = Communicator::Current();
	BlockSparseMatrix<double> wrongBlocks = comm.Rank() == comm.Size() - 1 ? BlockSparseMatrix<double>::FromCSR(S.Diag, 1) : blocks;
	EXPECT_THROW(bcsr.SetupBlockCSR(wrongBlocks, &S), fhhos4::Error);
}
