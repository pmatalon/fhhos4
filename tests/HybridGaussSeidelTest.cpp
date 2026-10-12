// Unit tests of the hybrid Gauss-Seidel sweeps (GaussSeidel and BlockGaussSeidel with hybrid = true, HybridChunks):
// Gauss-Seidel in the chunk of each thread, Jacobi between the chunks, hence deterministic; the sequential sweeps if the
// chunks are poorly separated.
#include <gtest/gtest.h>
#include <algorithm>
#include <memory>
#include <random>
#include <vector>
#include "Solver/FixedPoint/BlockGaussSeidel.h"
#include "Solver/FixedPoint/GaussSeidel.h"

namespace
{
	// Non-symmetric, diagonally dominant matrix of a 9-point stencil on an m x m grid, in blocks of blockSize DoFs (dense
	// blocks), with a coupling to a distant grid point in one block row out of 5: some block rows have non-zeros in other
	// chunks, others not. m = 50: more than 20000 non-zeros, so that the scalar sweeps run in parallel.
	SparseMatrix TestMatrix(int blockSize, int m = 50)
	{
		std::mt19937 gen(42);
		std::uniform_real_distribution<double> random(-1, 0);
		std::vector<Eigen::Triplet<double, SparseMatrixIndex>> triplets;
		int nb = m * m;
		auto addBlock = [&](int I, int J)
		{
			for (int k = 0; k < blockSize; k++)
			{
				for (int l = 0; l < blockSize; l++)
				{
					double value = I != J ? random(gen) : (k == l ? 20.0 * blockSize : 0.1 * random(gen));
					triplets.push_back({ (SparseMatrixIndex)(I * blockSize + k), (SparseMatrixIndex)(J * blockSize + l), value });
				}
			}
		};
		for (int y = 0; y < m; y++)
		{
			for (int x = 0; x < m; x++)
			{
				int I = y * m + x;
				for (int dy = -1; dy <= 1; dy++)
				{
					for (int dx = -1; dx <= 1; dx++)
					{
						if (y + dy >= 0 && y + dy < m && x + dx >= 0 && x + dx < m)
							addBlock(I, (y + dy) * m + x + dx);
					}
				}
				if (I % 5 == 0)
					addBlock(I, (I * 7919 + nb / 2) % nb);
			}
		}
		SparseMatrix A(nb * blockSize, nb * blockSize);
		A.setFromTriplets(triplets.begin(), triplets.end()); // sums the duplicates
		return A;
	}

	// The block rows in a random order (window = nb), or shuffled in windows of consecutive block rows: the chunks of the
	// order are then not ranges of columns, but still well separated
	std::vector<BigNumber> ShuffledOrder(BigNumber nb, BigNumber window)
	{
		std::vector<BigNumber> order(nb);
		for (BigNumber I = 0; I < nb; I++)
			order[I] = I;
		std::mt19937 gen(7);
		for (BigNumber begin = 0; begin < nb; begin += window)
			std::shuffle(order.begin() + begin, order.begin() + std::min(begin + window, nb), gen);
		return order;
	}

	// The block rows (and columns) of A renumbered in this order
	SparseMatrix Renumbered(const SparseMatrix& A, int blockSize, const std::vector<BigNumber>& order)
	{
		std::vector<BigNumber> newIndex(A.rows());
		for (BigNumber k = 0; k < (BigNumber)order.size(); k++)
		{
			for (int l = 0; l < blockSize; l++)
				newIndex[order[k] * blockSize + l] = k * blockSize + l;
		}
		std::vector<Eigen::Triplet<double, SparseMatrixIndex>> triplets;
		for (BigNumber i = 0; i < A.rows(); i++)
		{
			for (SparseMatrix::InnerIterator it(A, i); it; ++it)
				triplets.push_back({ (SparseMatrixIndex)newIndex[i], (SparseMatrixIndex)newIndex[it.col()], it.value() });
		}
		SparseMatrix B(A.rows(), A.cols());
		B.setFromTriplets(triplets.begin(), triplets.end());
		return B;
	}

	// Reference: the sweep sequence (the natural order, or order) split into nChunks contiguous chunks; in each chunk, a
	// block Gauss-Seidel sweep (backward: decreasing positions; never with an order, which all the sweeps follow) that
	// reads the values of the other chunks in x0
	Vector ReferenceSweep(const SparseMatrix& A, int blockSize, const Vector& b, const Vector& x0, const std::vector<BigNumber>* order, int nChunks, bool backward)
	{
		BigNumber nb = A.rows() / blockSize;
		backward = backward && !order;
		Vector x = x0;
		for (int c = 0; c < nChunks; c++)
		{
			BigNumber begin = nb * c / nChunks, end = nb * (c + 1) / nChunks;
			Vector y = x0;
			for (BigNumber s = 0; s < end - begin; s++)
			{
				BigNumber position = backward ? end - 1 - s : begin + s;
				BigNumber I = order ? (*order)[position] : position;
				// D_I y_I = b_I - sum_{J != I} A_IJ y_J
				Vector rhs = b.segment(I * blockSize, blockSize);
				DenseMatrix D = DenseMatrix::Zero(blockSize, blockSize);
				for (int k = 0; k < blockSize; k++)
				{
					for (SparseMatrix::InnerIterator it(A, I * blockSize + k); it; ++it)
					{
						if ((BigNumber)it.col() / blockSize == I)
							D(k, it.col() % blockSize) = it.value();
						else
							rhs[k] -= it.value() * y[it.col()];
					}
				}
				y.segment(I * blockSize, blockSize) = D.fullPivLu().solve(rhs);
			}
			for (BigNumber position = begin; position < end; position++)
			{
				BigNumber I = order ? (*order)[position] : position;
				x.segment(I * blockSize, blockSize) = y.segment(I * blockSize, blockSize);
			}
		}
		return x;
	}

	std::unique_ptr<IterativeSolver> CreateGaussSeidel(const SparseMatrix& A, int blockSize, bool hybrid, Direction direction, const std::vector<BigNumber>* order)
	{
		std::unique_ptr<IterativeSolver> gs;
		if (blockSize == 1)
		{
			GaussSeidel* scalar = new GaussSeidel(direction, hybrid);
			scalar->SweepOrder = order;
			gs.reset(scalar);
		}
		else
		{
			BlockGaussSeidel* block = new BlockGaussSeidel(blockSize, direction, hybrid);
			block->SweepOrder = order;
			gs.reset(block);
		}
		gs->StoppingCrit = StoppingCriteria::MaxIterations;
		gs->MaxIterations = 1;
		gs->PrintIterationResults = false;
		gs->Setup(A);
		return gs;
	}

	// One sweep from x0
	Vector Sweep(IterativeSolver& gs, const Vector& b, const Vector& x0, bool xEquals0)
	{
		Vector x = x0;
		gs.Solve(b, x, xEquals0);
		return x;
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

	// Runs test(A, blockSize, order, direction, b, x0, xEquals0, description) on every case: scalar and blocks, natural
	// order and SweepOrder, forward and backward, from 0 or not
	template <class Test>
	void ForEachCase(Test test)
	{
		for (int blockSize : { 1, 3 })
		{
			SparseMatrix A = TestMatrix(blockSize);
			std::vector<BigNumber> shuffled = ShuffledOrder(A.rows() / blockSize, 100);
			Vector b = Vector::LinSpaced(A.rows(), -1, 1);
			const std::vector<BigNumber>* orders[] = { nullptr, &shuffled };
			for (const std::vector<BigNumber>* order : orders)
			{
				for (Direction direction : { Direction::Forward, Direction::Backward })
				{
					for (bool xEquals0 : { true, false })
					{
						Vector x0 = xEquals0 ? Vector::Zero(A.rows()) : Vector(Vector::LinSpaced(A.rows(), 2, -3).array().sin().matrix());
						std::string description = "block size " + std::to_string(blockSize) + (order ? ", SweepOrder" : "") +
							(direction == Direction::Forward ? ", forward" : ", backward") + (xEquals0 ? ", from 0" : "");
						test(A, blockSize, order, direction, b, x0, xEquals0, description);
					}
				}
			}
		}
	}
}

// Gauss-Seidel in the chunks, Jacobi between them (a straightforward implementation), on 4 threads
TEST(HybridGaussSeidelTest, MatchesTheReference)
{
	Threads threads(4);
	ForEachCase([](const SparseMatrix& A, int blockSize, const std::vector<BigNumber>* order, Direction direction, const Vector& b, const Vector& x0, bool xEquals0, const std::string& description)
	{
		bool backward = direction == Direction::Backward;
		Vector reference = ReferenceSweep(A, blockSize, b, x0, order, 4, backward);
		Vector x = Sweep(*CreateGaussSeidel(A, blockSize, true, direction, order), b, x0, xEquals0);
		EXPECT_LT(MaxDiff(x, reference), 1e-13 * reference.lpNorm<Eigen::Infinity>()) << description;
		// The test discriminates: the sequential sweep differs from the reference
		Vector sequential = ReferenceSweep(A, blockSize, b, x0, order, 1, backward);
		EXPECT_GT(MaxDiff(sequential, reference), 1e-6 * reference.lpNorm<Eigen::Infinity>()) << description;
	});
}

// The same results from one run to the next
TEST(HybridGaussSeidelTest, Deterministic)
{
	Threads threads(4);
	ForEachCase([](const SparseMatrix& A, int blockSize, const std::vector<BigNumber>* order, Direction direction, const Vector& b, const Vector& x0, bool xEquals0, const std::string& description)
	{
		Vector first = Sweep(*CreateGaussSeidel(A, blockSize, true, direction, order), b, x0, xEquals0);
		for (int run = 0; run < 5; run++)
			EXPECT_TRUE(Identical(Sweep(*CreateGaussSeidel(A, blockSize, true, direction, order), b, x0, xEquals0), first)) << description << ", run " << run;
	});
}

// On 1 thread, the sequential sweep, and in the scalar case the residual deduced from the sweep as in the sequential
// smoother
TEST(HybridGaussSeidelTest, OneThreadIsSequential)
{
	Threads threads(1);
	ForEachCase([](const SparseMatrix& A, int blockSize, const std::vector<BigNumber>* order, Direction direction, const Vector& b, const Vector& x0, bool xEquals0, const std::string& description)
	{
		auto hybrid = CreateGaussSeidel(A, blockSize, true, direction, order);
		auto sequential = CreateGaussSeidel(A, blockSize, false, direction, order);
		EXPECT_TRUE(Identical(Sweep(*hybrid, b, x0, xEquals0), Sweep(*sequential, b, x0, xEquals0))) << description;
		EXPECT_EQ(hybrid->CanOptimizeResidualComputation(), sequential->CanOptimizeResidualComputation()) << description;
		if (sequential->CanOptimizeResidualComputation())
		{
			Vector x1 = x0, x2 = x0;
			hybrid->Solve(b, x1, xEquals0, true, false);
			sequential->Solve(b, x2, xEquals0, true, false);
			EXPECT_TRUE(Identical(x1, x2)) << description;
			EXPECT_TRUE(Identical(hybrid->Residual, sequential->Residual)) << description;
		}
	});
}

// The chunks follow the number of threads, which can change between two solves
TEST(HybridGaussSeidelTest, NumberOfThreadsChangesBetweenSolves)
{
	for (int blockSize : { 1, 3 })
	{
		SparseMatrix A = TestMatrix(blockSize);
		Vector b = Vector::LinSpaced(A.rows(), -1, 1);
		Vector x0 = Vector::LinSpaced(A.rows(), 2, -3).array().sin().matrix();
		auto gs = CreateGaussSeidel(A, blockSize, true, Direction::Forward, nullptr);
		for (int nThreads : { 4, 2, 4, 1, 3 })
		{
			Threads threads(nThreads);
			Vector reference = ReferenceSweep(A, blockSize, b, x0, nullptr, nThreads, false);
			EXPECT_LT(MaxDiff(Sweep(*gs, b, x0, false), reference), 1e-13 * reference.lpNorm<Eigen::Infinity>()) << "block size " << blockSize << ", " << nThreads << " threads";
		}
	}
}

// Poorly separated chunks (a random order, or a random numbering: most of the couplings between chunks): on 4 threads,
// the sequential sweep, bit for bit
TEST(HybridGaussSeidelTest, SequentialIfPoorlySeparated)
{
	Threads threads(4);
	for (int blockSize : { 1, 3 })
	{
		SparseMatrix A = TestMatrix(blockSize);
		std::vector<BigNumber> random = ShuffledOrder(A.rows() / blockSize, A.rows() / blockSize);
		SparseMatrix renumbered = Renumbered(A, blockSize, random);
		Vector b = Vector::LinSpaced(A.rows(), -1, 1);
		for (bool reorder : { true, false }) // SweepOrder, or the natural order of the renumbered matrix
		{
			const SparseMatrix& M = reorder ? A : renumbered;
			const std::vector<BigNumber>* order = reorder ? &random : nullptr;
			for (Direction direction : { Direction::Forward, Direction::Backward })
			{
				for (bool xEquals0 : { true, false })
				{
					Vector x0 = xEquals0 ? Vector::Zero(A.rows()) : Vector(Vector::LinSpaced(A.rows(), 2, -3).array().sin().matrix());
					Vector hybrid = Sweep(*CreateGaussSeidel(M, blockSize, true, direction, order), b, x0, xEquals0);
					Vector sequential = Sweep(*CreateGaussSeidel(M, blockSize, false, direction, order), b, x0, xEquals0);
					EXPECT_TRUE(Identical(hybrid, sequential)) << "block size " << blockSize << (reorder ? ", SweepOrder" : ", renumbered")
						<< (direction == Direction::Forward ? ", forward" : ", backward") << (xEquals0 ? ", from 0" : "");
				}
			}
		}
	}
}
