#include <gtest/gtest.h>
#include "Solver/FixedPoint/BlockChebyshev.h"
#include "support/CartesianHybridSystem.h"
#include "support/PartitionedSystem.h"

using namespace fhhos4_tests;

namespace
{
	std::unique_ptr<BlockChebyshev> Smoother(ChebyshevKind kind, int degree)
	{
		auto chebyshev = std::make_unique<BlockChebyshev>(1, kind);
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
