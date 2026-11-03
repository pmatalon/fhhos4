#pragma once
#include "BlockSOR.h"
#include "BlockDiagonalSolver.h"
#include "DiagonalBlockPositions.h"
#include "HybridChunks.h"
using namespace std;

// The matrix of the setup in CSR in double (Setup()) or in single precision (SetupSinglePrecision(), the sweeps of float
// vectors, their sums in double), or stored by blocks in either precision (SetupBlockCSR(): the same sums in the same
// order, the same results as in CSR)
class BlockGaussSeidel : public IterativeSolver
{
protected:
	int _blockSize;
	Direction _direction;
	bool _hybrid = false;
	bool _l1 = false; // l1 Gauss-Seidel (see IterativeSolver::L1Theta)

	BlockDiagonalSolver diagBlockSolver; // solves the systems with the diagonal blocks
	DiagonalBlockPositions _diagBlock;   // positions of the diagonal block in each row (in CSR)
	// For each block row, whether its rows have the same non-zero pattern (always the case for the HHO matrices; in CSR)
	vector<char> _samePattern;
	HybridChunks _chunks; // of the hybrid sweeps (see Chunks())
	// The l1 term of each row (see L1Theta) in the factorized diagonal blocks; empty: none
	Vector _l1Terms;
public:
	// If set, every sweep (forward or backward, pre- or post-smoothing) processes the block rows in this order: the
	// downwind order of a convection-dominated matrix (DownwindOrder()), which all the sweeps then follow
	const vector<BigNumber>* SweepOrder = nullptr;

	BlockGaussSeidel(int blockSize, Direction direction, bool hybrid, bool l1 = false)
	{
		this->_blockSize = blockSize;
		this->_direction = direction;
		this->_hybrid = hybrid;
		this->_l1 = l1;
	}

	virtual void Serialize(ostream& os) const override
	{
		if (_l1)
			os << "l1 ";
		if (_hybrid)
			os << "Hybrid ";
		if (_direction == Direction::Symmetric)
			os << "Symmetric ";
		os << "Block Gauss-Seidel";

		if (_direction != Direction::Symmetric || _blockSize != 1)
			os << " (";
		if (_blockSize != 1)
		{
			os << "blockSize=" << _blockSize;
			if (_direction != Direction::Symmetric)
				os << ", ";
		}
		if (_direction == Direction::Forward)
			os << "direction=forward";
		else if (_direction == Direction::Backward)
			os << "direction=backward";
		else if (_direction == Direction::AlternatingForwardFirst)
			os << "alternating directions, forward first";
		else if (_direction == Direction::AlternatingBackwardFirst)
			os << "alternating directions, backward first";

		if (_direction != Direction::Symmetric || _blockSize != 1)
			os << ")";
	}

	//------------------------------------------------//
	// Big assumption: the matrix must be symmetric!  //
	//------------------------------------------------//

	void Setup(const SparseMatrix& A) override
	{
		IterativeSolver::Setup(A);
		SetStorage(&A, nullptr, nullptr, nullptr, nullptr, nullptr);
		SetupLocal(A);
	}

	// On several processes: Gauss-Seidel in the rows of the process, Jacobi between the processes (see SweepRhs())
	void Setup(const DistributedMatrix& A) override
	{
		SetDistributedMatrix(A);
		SetStorage(&A.Diag, nullptr, nullptr, nullptr, &A, nullptr);
		SetupLocal(A.Diag);
	}

	// Mixed precision: the same sweeps of float vectors (SolveInSinglePrecision()), their sums in double
	bool SupportsSinglePrecision() const override
	{
		return true;
	}

	void SetupSinglePrecision(const FloatSparseMatrix& A) override
	{
		SetStorage(nullptr, &A, nullptr, nullptr, nullptr, nullptr);
		SetupLocal(A);
	}

	void SetupSinglePrecision(const DistributedMatrixOf<float>& A) override
	{
		SetStorage(nullptr, &A.Diag, nullptr, nullptr, nullptr, &A);
		SetupLocal(A.Diag);
	}

	// The matrix stored by blocks of the block size of the smoother (block CSR); on several processes, the block Diag of
	// distributed
	bool SupportsBlockCSR() const override
	{
		return true;
	}

	void SetupBlockCSR(const BlockSparseMatrix<double>& A, const DistributedMatrix* distributed) override
	{
		SetStorage(nullptr, nullptr, &A, nullptr, distributed, nullptr);
		SetupLocal(A);
	}

	void SetupBlockCSR(const BlockSparseMatrix<float>& A, const DistributedMatrixOf<float>* distributed) override
	{
		SetStorage(nullptr, nullptr, nullptr, &A, nullptr, distributed);
		SetupLocal(A);
	}

private:
	template <typename T>
	void SetupLocal(const SparseMatrixOf<T>& A)
	{
		this->diagBlockSolver.Setup(A, _blockSize, BlockDiagonalSolver::SmootherMethod);
		_diagBlock.Setup(A, _blockSize);
		SetupSamePattern(A);
		SetupChunksAndCost();
	}

	template <typename T>
	void SetupLocal(const BlockSparseMatrix<T>& A)
	{
		if (A.BlockSize != _blockSize)
			Utils::FatalError("The block CSR storage of a block Gauss-Seidel smoother requires blocks of its block size (" + to_string(_blockSize) + "), got " + to_string(A.BlockSize) + ".");
		this->diagBlockSolver.Setup(A, _blockSize, BlockDiagonalSolver::SmootherMethod);
		_diagBlock = DiagonalBlockPositions();
		_samePattern.clear();
		SetupChunksAndCost();
	}

	void SetupChunksAndCost()
	{
		_l1Terms = Vector();
		_chunks.Clear();
		if (_hybrid)
			Chunks();
		else
			SetupL1Terms();

		this->SetupComputationalWork = MatrixRows() / _blockSize * 2.0/3.0*pow(_blockSize, 3)*1e-6;
	}

public:
	// The setup of another block Gauss-Seidel of the same block size set up for the same matrix in the same storage (the
	// pre-smoother of the level): its factorized diagonal blocks, the positions of the blocks and, with the same order,
	// its hybrid chunks
	bool SetupFrom(const SparseMatrix& A, const IterativeSolver& other) override
	{
		return SetupFromOther(other, [&] { SetStorage(&A, nullptr, nullptr, nullptr, nullptr, nullptr); });
	}

	bool SetupFrom(const DistributedMatrix& A, const IterativeSolver& other) override
	{
		return SetupFromOther(other, [&] { SetStorage(&A.Diag, nullptr, nullptr, nullptr, &A, nullptr); });
	}

	bool SetupSinglePrecisionFrom(const FloatSparseMatrix& A, const IterativeSolver& other) override
	{
		return SetupFromOther(other, [&] { SetStorage(nullptr, &A, nullptr, nullptr, nullptr, nullptr); });
	}

	bool SetupSinglePrecisionFrom(const DistributedMatrixOf<float>& A, const IterativeSolver& other) override
	{
		return SetupFromOther(other, [&] { SetStorage(nullptr, &A.Diag, nullptr, nullptr, nullptr, &A); });
	}

	bool SetupBlockCSRFrom(const BlockSparseMatrix<double>& A, const DistributedMatrix* distributed, const IterativeSolver& other) override
	{
		return SetupFromOther(other, [&] { SetStorage(nullptr, nullptr, &A, nullptr, distributed, nullptr); });
	}

	bool SetupBlockCSRFrom(const BlockSparseMatrix<float>& A, const DistributedMatrixOf<float>* distributed, const IterativeSolver& other) override
	{
		return SetupFromOther(other, [&] { SetStorage(nullptr, nullptr, nullptr, &A, nullptr, distributed); });
	}

private:
	// setStorage(): the matrix of this setup. False (the caller then sets this smoother up) if other is not set up for it.
	template <class SetStorageFunction>
	bool SetupFromOther(const IterativeSolver& other, SetStorageFunction setStorage)
	{
		const BlockGaussSeidel* o = dynamic_cast<const BlockGaussSeidel*>(&other);
		if (!o || o->_blockSize != _blockSize || this->ComputeExactSolution)
			return false;
		setStorage();
		if (!SameStorage(*o))
			return false;
		CopySetup(o);
		return true;
	}

	void CopySetup(const BlockGaussSeidel* o)
	{
		diagBlockSolver = o->diagBlockSolver;
		_l1Terms = o->_l1Terms;
		_diagBlock = o->_diagBlock;
		_samePattern = o->_samePattern;
		_chunks.Clear();
		if (_hybrid && o->_hybrid && o->SweepOrder == SweepOrder)
			_chunks = o->_chunks;
		if (_hybrid)
			Chunks();
		SetupL1Terms(); // those of this smoother (o's: l1 or not, its chunks or not)
		this->SetupComputationalWork = o->SetupComputationalWork;
	}

	IterationResult ExecuteOneIteration(const Vector& b, Vector& x, bool& xEquals0, bool computeResidual, bool computeAx, const IterationResult& oldResult) override
	{
		IterationResult result(oldResult);
		assert(!computeResidual && !computeAx);
		Iterate(b, x, xEquals0);
		result.SetX(x);
		result.AddWorkInFlops(IterationFlops());
		return result;
	}

public:
	// The iterations in single precision (the residual and A x computed by the caller: CanOptimizeResidualComputation()
	// is false)
	void SolveInSinglePrecision(const FloatVector& b, FloatVector& x, bool xEquals0, FloatVector* residual, FloatVector* Ax) override
	{
		assert(this->MaxIterations > 0 && InSinglePrecision() && !residual && !Ax);
		for (this->IterationCount = 0; this->IterationCount < this->MaxIterations; this->IterationCount++)
			Iterate(b, x, xEquals0);
		this->SolvingComputationalWork = this->MaxIterations * IterationFlops() * 1e-6;
	}

private:
	// One iteration (one sweep, two if symmetric) in precision T
	template <typename T>
	void Iterate(const VectorOf<T>& b, VectorOf<T>& x, bool& xEquals0)
	{
		if (!_hybrid)
		{
			if (_direction == Direction::Forward)
				Sweep(b, x, false, xEquals0);
			else if (_direction == Direction::Backward)
				Sweep(b, x, true, xEquals0);
			else if (_direction == Direction::Symmetric)
			{
				Sweep(b, x, false, xEquals0); // Forward
				Sweep(b, x, true, false);     // Backward
			}
			else if (_direction == Direction::AlternatingForwardFirst || _direction == Direction::AlternatingBackwardFirst)
			{
				int modulo = _direction == Direction::AlternatingForwardFirst ? 0 : 1;
				if (this->IterationCount % 2 == modulo)
					Sweep(b, x, false, xEquals0); // Forward
				else
					Sweep(b, x, true, xEquals0);  // Backward
			}
			else
				Utils::FatalError("direction not managed");
		}
		else // Hybrid
		{
			if (_direction == Direction::Forward)
				HybridSweep(b, x, false, xEquals0);
			else if (_direction == Direction::Backward)
				HybridSweep(b, x, true, xEquals0);
			else if (_direction == Direction::Symmetric)
			{
				HybridSweep(b, x, false, xEquals0); // Forward
				HybridSweep(b, x, true, false);     // Backward
			}
			else if (_direction == Direction::AlternatingForwardFirst || _direction == Direction::AlternatingBackwardFirst)
			{
				int modulo = _direction == Direction::AlternatingForwardFirst ? 0 : 1;
				if (this->IterationCount % 2 == modulo)
					HybridSweep(b, x, false, xEquals0); // Forward
				else
					HybridSweep(b, x, true, xEquals0);  // Backward
			}
			else
				Utils::FatalError("direction not managed");
		}

		xEquals0 = false;
	}

	Flops IterationFlops() const
	{
		double sweepWork = 2 * MatrixNonZeros() + MatrixRows() / _blockSize * pow(_blockSize, 2);
		return _direction == Direction::Symmetric ? 2*sweepWork : sweepWork;
	}

	// Which off-diagonal blocks of a block row are multiplied by x: in the first sweep from x = 0, those of the block
	// rows not processed yet multiply zeros (the upper ones in a forward sweep, the lower ones in a backward sweep)
	enum class OffDiagonalBlocks { Both, Lower, Upper };

	// The work vectors of the sweeps, allocated once per sweep (per thread in the hybrid sweeps): the right-hand side of a
	// block row (its sums in double) and the solution of its diagonal block; their copies in float for the LU solves in
	// single precision
	struct SweepWork
	{
		Vector Rhs, Solution;
		FloatVector SingleRhs, SingleSolution;
		SweepWork(int blockSize) : Rhs(blockSize), Solution(blockSize), SingleRhs(blockSize), SingleSolution(blockSize) {}
	};

	template <typename T>
	void Sweep(const VectorOf<T>& givenB, VectorOf<T>& x, bool backward, bool xEquals0)
	{
		VectorOf<T> distributedRhs;
		const VectorOf<T>& b = SweepRhs(givenB, x, xEquals0, distributedRhs, DistributedOf<T>(), &_l1Terms);
		BigNumber n = Rows<T>();
		SweepWork work(_blockSize);
		if (xEquals0 && (SweepOrder || x.size() != n))
			x = VectorOf<T>::Zero(n);
		SweepRange(b, x, backward, xEquals0, 0, n / _blockSize, work);
	}

	// The block rows at the positions [begin, end) of the sweep sequence: the natural order (processed backward if
	// backward), or SweepOrder (all the sweeps follow it)
	template <typename T>
	void SweepRange(const VectorOf<T>& b, VectorOf<T>& x, bool backward, bool xEquals0, BigNumber begin, BigNumber end, SweepWork& work)
	{
		if (SweepOrder)
		{
			// The block rows not processed yet are not the upper ones: from x = 0, all the off-diagonal blocks are
			// multiplied (by zeros for those)
			for (BigNumber k = begin; k < end; ++k)
				ProcessBlockRow((*SweepOrder)[k], b, x, work, OffDiagonalBlocks::Both, true);
			return;
		}
		OffDiagonalBlocks blocks = !xEquals0 ? OffDiagonalBlocks::Both : (backward ? OffDiagonalBlocks::Upper : OffDiagonalBlocks::Lower);
		if (backward)
		{
			for (BigNumber i = end; i-- > begin; )
				ProcessBlockRow(i, b, x, work, blocks, true);
		}
		else
		{
			for (BigNumber i = begin; i < end; ++i)
				ProcessBlockRow(i, b, x, work, blocks, false);
		}
	}

	// Hybrid sweep (see HybridChunks): each thread sweeps its own contiguous chunk of the block rows (or of SweepOrder) by
	// Gauss-Seidel, and reads the values of the other chunks from the start of the sweep (Jacobi between the chunks). The
	// interior block rows are processed as in Sweep(). The result depends on the number of threads, not on their
	// timing; on one thread, or if the chunks are poorly separated (HybridChunks::Parallel()), it is Sweep(). From x = 0,
	// the products by the values of the other chunks and by those not computed yet are skipped, as in Sweep() (except in
	// SweepOrder).
	template <typename T>
	void HybridSweep(const VectorOf<T>& givenB, VectorOf<T>& x, bool backward, bool xEquals0)
	{
		if (!Chunks().Parallel())
		{
			Sweep(givenB, x, backward, xEquals0);
			return;
		}
		VectorOf<T> distributedRhs;
		const VectorOf<T>& b = SweepRhs(givenB, x, xEquals0, distributedRhs, DistributedOf<T>(), &_l1Terms);
		BigNumber n = Rows<T>();
		if (xEquals0 && (SweepOrder || x.size() != n))
			x = VectorOf<T>::Zero(n); // as Sweep()
		bool interleaveRows = backward || SweepOrder; // as in Sweep()
		bool backwardChunks = backward && !SweepOrder;
		auto products = !xEquals0 || SweepOrder ? HybridChunks::InternalProducts::Both : (backward ? HybridChunks::InternalProducts::Upper : HybridChunks::InternalProducts::Lower);
		const BlockSparseMatrix<T>* blocks = BlockMatrixOf<T>();
		#pragma omp parallel num_threads(_chunks.Count())
		{
			SweepWork work(_blockSize); // allocated once per thread
			if (!xEquals0)
			{
				_chunks.SaveExported(x);
				#pragma omp barrier
			}
			_chunks.ForEachChunk([&](int c)
			{
				_chunks.SweepChunk(c, backwardChunks,
					[&](BigNumber begin, BigNumber end) { SweepRange(b, x, backward, xEquals0, begin, end, work); },
					[&](BigNumber q, BigNumber blockRow)
					{
						if (blocks)
							DispatchBlockSize(_blockSize, [&](auto size) { _chunks.BoundaryRowRhs<T, decltype(size)::value>(*blocks, c, q, blockRow, b, x, xEquals0, products, work.Rhs.data()); });
						else
							_chunks.BoundaryRowRhs(MatrixOf<T>(), c, q, blockRow, b, x, xEquals0, products, work.Rhs.data(), interleaveRows && _samePattern[blockRow]);
						SolveDiagonalBlock(blockRow, work, x.data() + blockRow * _blockSize);
					});
			});
		}
	}

	// The chunks of the hybrid sweeps: one per thread. Set up with the matrix, again if the number of threads or
	// SweepOrder changes.
	const HybridChunks& Chunks()
	{
		int nChunks = Parallelism::NThreads();
		if (!_chunks.IsSetUp(nChunks, SweepOrder))
		{
			if (_blockMatrix)
				_chunks.Setup(*_blockMatrix, SweepOrder, nChunks);
			else if (_singleBlockMatrix)
				_chunks.Setup(*_singleBlockMatrix, SweepOrder, nChunks);
			else if (_singleMatrix)
				_chunks.Setup(*_singleMatrix, _blockSize, SweepOrder, nChunks);
			else
				_chunks.Setup(*this->Matrix, _blockSize, SweepOrder, nChunks);
			SetupL1Terms();
		}
		return _chunks;
	}

	// The l1 terms (L1Theta, if l1) of the rows coupled to other processes, and to other chunks if the hybrid sweeps run
	// in parallel; the diagonal blocks whose l1 terms change factorized again
	void SetupL1Terms()
	{
		if (InSinglePrecision())
			SetupL1Terms<float>();
		else
			SetupL1Terms<double>();
	}

	template <typename T>
	void SetupL1Terms()
	{
		double theta = _l1 ? L1Theta : 0;
		if (theta == 0 && _l1Terms.size() == 0)
			return;
		BigNumber n = Rows<T>();
		Vector terms = L1Terms<T>(theta, [&](Vector& sums)
		{
			if (!_hybrid || !_chunks.Parallel())
				return;
			if (const BlockSparseMatrix<T>* blocks = BlockMatrixOf<T>())
				_chunks.AddExternalAbsoluteSums(*blocks, sums);
			else
				_chunks.AddExternalAbsoluteSums(MatrixOf<T>(), sums);
		});
		Vector old = _l1Terms.size() > 0 ? _l1Terms : Vector::Zero(n);
		vector<BigNumber> changed;
		for (BigNumber I = 0; I < n / _blockSize; I++)
		{
			if ((terms.segment(I * _blockSize, _blockSize).array() != old.segment(I * _blockSize, _blockSize).array()).any())
				changed.push_back(I);
		}
		diagBlockSolver.Refactorize<T>(changed, [&](BigNumber I)
		{
			DenseMatrix D = DiagonalBlock<T>(I);
			D.diagonal() += terms.segment(I * _blockSize, _blockSize);
			return D;
		});
		_l1Terms = (terms.array() != 0).any() ? terms : Vector();
	}

	// The diagonal block I of the matrix in precision T, in double
	template <typename T>
	DenseMatrix DiagonalBlock(BigNumber I) const
	{
		if (const BlockSparseMatrix<T>* blocks = BlockMatrixOf<T>())
			return blocks->DiagonalBlock(I);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> D = MatrixOf<T>().block(I * _blockSize, I * _blockSize, _blockSize, _blockSize);
		return D.template cast<double>();
	}

	template <typename T>
	void SetupSamePattern(const SparseMatrixOf<T>& A)
	{
		const SparseMatrixIndex* outer = A.outerIndexPtr();
		const SparseMatrixIndex* col = A.innerIndexPtr();
		BigNumber nb = A.rows() / _blockSize;
		_samePattern.resize(nb);
		#pragma omp parallel for
		for (BigNumber iBlock = 0; iBlock < nb; iBlock++)
		{
			BigNumber i0 = iBlock * _blockSize;
			SparseMatrixIndex length = outer[i0 + 1] - outer[i0];
			bool same = true;
			for (int k = 1; k < _blockSize && same; k++)
			{
				BigNumber i = i0 + k;
				same = outer[i + 1] - outer[i] == length && equal(col + outer[i], col + outer[i + 1], col + outer[i0]);
			}
			_samePattern[iBlock] = same;
		}
	}

	// interleaveRows: whether to use BlockRowRhsSamePattern() when possible (in CSR). Measured on the HHO matrices
	// (U-AMG, k=2): in a forward sweep, BlockRowRhs() is faster (it reads the matrix as one sequential stream, which the
	// hardware prefetcher follows); in a backward sweep, BlockRowRhsSamePattern() is (BlockRowRhs() jumps back at each
	// block row, and the prefetcher doesn't follow: reading the rows together starts their memory accesses together).
	// In block CSR, the rows of a block row are read together.
	template <typename T>
	void ProcessBlockRow(BigNumber currentBlockRow, const VectorOf<T>& b, VectorOf<T>& x, SweepWork& work, OffDiagonalBlocks blocks, bool interleaveRows)
	{
		// BlockRow i: [ --- Li --- | Di | --- Ui --- ]

		/*
		auto nb = A.rows() / _blockSize;
		// x_new                                               x_new                                  x_old                                               x_old
		x.segment(i * _blockSize, _blockSize) = -_omega * Li * x.head(i * _blockSize) - _omega * Ui * x.tail((nb - i - 1) * _blockSize) + (1 - _omega)*Di*x.segment(i * _blockSize, _blockSize) + _omega * bi;
		x.segment(i * _blockSize, _blockSize) = this->invD.block(i * _blockSize, 0, _blockSize, _blockSize) * x.segment(i * _blockSize, _blockSize);*/

		// rhs = b_i - (Li | Ui) * x
		double* rhs = work.Rhs.data();
		if (const BlockSparseMatrix<T>* M = BlockMatrixOf<T>())
			DispatchBlockSize(_blockSize, [&](auto size) { BlockRowRhs<T, decltype(size)::value>(*M, currentBlockRow, b, x, rhs, blocks); });
		else if (interleaveRows && _samePattern[currentBlockRow])
		{
			switch (_blockSize)
			{
				case 2: BlockRowRhsSamePattern<T, 2>(currentBlockRow, b, x, rhs, blocks); break;
				case 3: BlockRowRhsSamePattern<T, 3>(currentBlockRow, b, x, rhs, blocks); break;
				case 4: BlockRowRhsSamePattern<T, 4>(currentBlockRow, b, x, rhs, blocks); break;
				case 5: BlockRowRhsSamePattern<T, 5>(currentBlockRow, b, x, rhs, blocks); break;
				case 6: BlockRowRhsSamePattern<T, 6>(currentBlockRow, b, x, rhs, blocks); break;
				default: BlockRowRhs(MatrixOf<T>(), currentBlockRow, b, x, rhs, blocks);
			}
		}
		else
			BlockRowRhs(MatrixOf<T>(), currentBlockRow, b, x, rhs, blocks);

		SolveDiagonalBlock(currentBlockRow, work, x.data() + currentBlockRow * _blockSize);
	}

	// x_I = D_I^-1 work.Rhs (destroyed): in double, the solve of BlockDiagonalSolver; in single precision, the inverse
	// applied with the sums in double (ApplyInverse()), or the LU solve in float
	template <typename T>
	void SolveDiagonalBlock(BigNumber I, SweepWork& work, T* x_I) const
	{
		if constexpr (std::is_same_v<T, double>)
		{
			this->diagBlockSolver.Solve(I, work.Rhs, work.Solution);
			for (int k = 0; k < _blockSize; k++)
				x_I[k] = work.Solution[k];
		}
		else if (diagBlockSolver.Method() == DiagBlockSolveMethod::Inverse)
		{
			DispatchBlockSize(_blockSize, [&](auto size) { diagBlockSolver.ApplyInverse<T, decltype(size)::value>(I, work.Rhs.data(), work.Solution.data()); });
			for (int k = 0; k < _blockSize; k++)
				x_I[k] = (T)work.Solution[k];
		}
		else
		{
			work.SingleRhs = work.Rhs.template cast<T>();
			this->diagBlockSolver.Solve(I, work.SingleRhs, work.SingleSolution);
			for (int k = 0; k < _blockSize; k++)
				x_I[k] = work.SingleSolution[k];
		}
	}

	// rhs = b_i - (Li | Ui) * x for the block row i: for each row, b_k - a_kj*x_j - ..., the products subtracted one by one
	// in the order of the row. The sum of a row is a chain of dependent operations: it is kept in a register.
	// The blocks excluded by 'blocks' multiply zeros: skipped (s - 0 = s, exactly, up to the sign of a zero s).
	// In precision T: the sums in double.
	template <typename T>
	void BlockRowRhs(const SparseMatrixOf<T>& A, BigNumber iBlock, const VectorOf<T>& b, const VectorOf<T>& x, double* rhs, OffDiagonalBlocks blocks) const
	{
		const SparseMatrixIndex* outer = A.outerIndexPtr();
		const SparseMatrixIndex* col = A.innerIndexPtr();
		const T* val = A.valuePtr();
		const T* xp = x.data();
		for (int k = 0; k < _blockSize; k++)
		{
			BigNumber i = iBlock * _blockSize + k;
			double s = b[i];
			if (blocks != OffDiagonalBlocks::Upper) // Li
			{
				for (SparseMatrixIndex p = outer[i]; p < _diagBlock.Begin[i]; p++)
					s -= (double)val[p] * xp[col[p]];
			}
			if (blocks != OffDiagonalBlocks::Lower) // Ui
			{
				for (SparseMatrixIndex p = _diagBlock.End[i]; p < outer[i + 1]; p++)
					s -= (double)val[p] * xp[col[p]];
			}
			rhs[k] = s;
		}
	}

	// Same as BlockRowRhs (same operations in the same order), when the rows of the block row have the same pattern:
	// they are then stored one after the other with the same length, and processed together, column by column.
	// Their sums are independent chains of operations, which the CPU interleaves.
	template <typename T, int BlockSize>
	void BlockRowRhsSamePattern(BigNumber iBlock, const VectorOf<T>& b, const VectorOf<T>& x, double* rhs, OffDiagonalBlocks blocks) const
	{
		const SparseMatrixOf<T>& A = MatrixOf<T>();
		BigNumber i0 = iBlock * BlockSize;
		SparseMatrixIndex start = A.outerIndexPtr()[i0];
		SparseMatrixIndex length = A.outerIndexPtr()[i0 + 1] - start;
		SparseMatrixIndex diagBegin = _diagBlock.Begin[i0] - start;
		SparseMatrixIndex diagEnd = _diagBlock.End[i0] - start;
		const SparseMatrixIndex* col = A.innerIndexPtr() + start;
		const T* val[BlockSize]; // coefficients of each row
		double s[BlockSize];
		for (int k = 0; k < BlockSize; k++)
		{
			val[k] = A.valuePtr() + start + k * length;
			s[k] = b[i0 + k];
		}
		const T* xp = x.data();

		if (blocks != OffDiagonalBlocks::Upper) // Li
		{
			for (SparseMatrixIndex p = 0; p < diagBegin; p++)
			{
				double x_j = xp[col[p]];
				for (int k = 0; k < BlockSize; k++)
					s[k] -= (double)val[k][p] * x_j;
			}
		}
		if (blocks != OffDiagonalBlocks::Lower) // Ui
		{
			for (SparseMatrixIndex p = diagEnd; p < length; p++)
			{
				double x_j = xp[col[p]];
				for (int k = 0; k < BlockSize; k++)
					s[k] -= (double)val[k][p] * x_j;
			}
		}
		for (int k = 0; k < BlockSize; k++)
			rhs[k] = s[k];
	}

	// The same on the matrix stored by blocks: the blocks of the block row but the diagonal one (and those excluded by
	// 'blocks'), each row summed in the order of its columns, as in CSR; B: the block size if it is known at compile
	// time, 0 otherwise
	template <typename T, int B>
	void BlockRowRhs(const BlockSparseMatrix<T>& A, BigNumber I, const VectorOf<T>& b, const VectorOf<T>& x, double* rhs, OffDiagonalBlocks blocks) const
	{
		const int n = B > 0 ? B : _blockSize;
		for (int r = 0; r < n; r++)
			rhs[r] = b[I * n + r];
		const T* block = A.Values.data() + A.ValueStart[I];
		const T* xp = x.data();
		for (SparseMatrixIndex p = A.RowStart[I]; p < A.RowStart[I + 1]; p++, block += n * n)
		{
			BigNumber J = A.Cols[p];
			if (J == I || (J < I && blocks == OffDiagonalBlocks::Upper) || (J > I && blocks == OffDiagonalBlocks::Lower))
				continue;
			const T* x_J = xp + J * n;
			for (int l = 0; l < n; l++)
			{
				double x_l = x_J[l];
				const T* column = block + l * n;
				for (int r = 0; r < n; r++)
					rhs[r] -= (double)column[r] * x_l;
			}
		}
	}
};
