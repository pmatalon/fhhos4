#pragma once
#include "../IterativeSolver.h"
#include "BlockDiagonalSolver.h"
using namespace std;

enum class ChebyshevKind : unsigned
{
	// Chebyshev polynomial of the first kind, minimal on [lambdaMax / Ratio, lambdaMax]: the classical Chebyshev smoother
	// (Adams, Brezina, Hu, Tuminaro, J. Comput. Phys. 188, 2003). The low end of the spectrum is left to the coarse levels.
	First,
	// Chebyshev polynomial of the fourth kind, which minimizes max sqrt(lambda / lambdaMax) |p(lambda)| on [0, lambdaMax],
	// the smoothing factor of the convergence bound of the V-cycle: no lower bound to choose (Lottes, Numer. Linear
	// Algebra Appl. 30, 2023; Phillips, Fischer, arXiv:2210.03179). Degree 1: Jacobi with omega = 4 / (3 lambdaMax).
	Fourth
};

// Chebyshev smoother preconditioned by the block Jacobi iteration: after m iterations (MaxIterations: the degree of the
// polynomial), the error is p(D^-1 A) times the initial error, where D is the block diagonal of A and p the Chebyshev
// polynomial of the kind, with p(0) = 1. Each iteration costs a product by A and the solves of the diagonal blocks, as
// a block Jacobi iteration, but its block rows are independent: it runs in parallel, and its results don't depend on
// the number of threads. Each pass over the matrix fuses the product with the vector updates of the iteration.
// lambdaMax: the largest eigenvalue of D^-1 A, estimated in the setup by Lanczos, times EigenvalueSafetyFactor (the
// estimate is from below; above lambdaMax, |p| grows quickly).
// Requires a symmetric positive (semi-)definite matrix (real spectrum of D^-1 A); the smoother is then symmetric in the
// A-inner product: the same in pre- and post-smoothing.
class BlockChebyshev : public IterativeSolver
{
public:
	// Ratio lambdaMax / lambdaMin of the interval of the first kind (-cheb-ratio). 10: PETSc's interval, [0.1, 1.1] x the
	// estimate; better than Adams et al.'s 30 with U-AMG (PERFORMANCE.md)
	inline static double DefaultRatio = 10;
	// Lanczos iterations of the estimation of lambdaMax, and the factor applied to the estimate
	static constexpr int EigenvalueEstimateIterations = 10;
	static constexpr double EigenvalueSafetyFactor = 1.1;

private:
	int _blockSize;
	ChebyshevKind _kind;
	double _ratio;
	double _lambdaMaxEstimate = 0; // Lanczos estimate of the largest eigenvalue of D^-1 A (0 before the setup)

	// Blocks of size > 1: their inverses (DiagBlockSolveMethod::Inverse, computed in double), applied by a small product in
	// the passes (ApplyInverseBlock()): in block CSR, passes 1.7x faster than with the LU solves in 2D (3x3 blocks), 1.25x
	// in 3D (6x6; PERFORMANCE.md). Their rounding doesn't limit the accuracy of the solution: the smoother is a part of
	// the preconditioner (the outer Krylov method checks the true residual); the blocks of an SPD matrix are invertible.
	BlockDiagonalSolver _diagBlockSolver;

	// The inverses of the diagonal coefficients (blocks of size 1) and the work vectors of the iterations (kept from one
	// call to the next), in precision T: double, or float (a matrix in single precision, SetupSinglePrecision())
	template <typename T>
	struct Precision
	{
		VectorOf<T> InvDiag;
		VectorOf<T> d, dNext, r;
	};
	Precision<double> _double;
	Precision<float> _single;
	const FloatSparseMatrix* _singleMatrix = nullptr; // the matrix in single precision (SetupSinglePrecision())
	// The matrix stored by blocks (block CSR, SetupBlockCSR()), in double or in single precision; null: in CSR
	const BlockSparseMatrix<double>* _blockMatrix = nullptr;
	const BlockSparseMatrix<float>* _singleBlockMatrix = nullptr;

	template <typename T>
	Precision<T>& In()
	{
		if constexpr (std::is_same_v<T, float>)
			return _single;
		else
			return _double;
	}

	template <typename T>
	const Precision<T>& In() const
	{
		if constexpr (std::is_same_v<T, float>)
			return _single;
		else
			return _double;
	}

	template <typename T>
	const SparseMatrixOf<T>& MatrixOf() const
	{
		if constexpr (std::is_same_v<T, float>)
			return *_singleMatrix;
		else
			return *this->Matrix;
	}

	template <typename T>
	const BlockSparseMatrix<T>* BlockMatrixOf() const
	{
		if constexpr (std::is_same_v<T, float>)
			return _singleBlockMatrix;
		else
			return _blockMatrix;
	}

	// The numbers of rows and of non-zeros (stored coefficients) of the matrix in precision T, in CSR or block CSR
	template <typename T>
	BigNumber Rows() const
	{
		return BlockMatrixOf<T>() ? BlockMatrixOf<T>()->Rows() : MatrixOf<T>().rows();
	}

	template <typename T>
	BigNumber NonZeros() const
	{
		return BlockMatrixOf<T>() ? BlockMatrixOf<T>()->NonZeros() : MatrixOf<T>().nonZeros();
	}
public:

	BlockChebyshev(int blockSize, ChebyshevKind kind, double ratio = DefaultRatio)
	{
		_blockSize = blockSize;
		_kind = kind;
		_ratio = ratio;
		if (_kind == ChebyshevKind::First && !(_ratio > 1))
			Utils::FatalError("Chebyshev smoother: the ratio lambdaMax / lambdaMin of its interval must be > 1, got " + to_string(_ratio) + ".");
	}

	virtual void Serialize(ostream& os) const override
	{
		os << (_blockSize == 1 ? "Chebyshev" : "Block Chebyshev") << " (";
		if (_blockSize != 1)
			os << "blockSize=" << _blockSize << ", ";
		if (_kind == ChebyshevKind::First)
			os << "1st kind on [lambdaMax/" << _ratio << ", lambdaMax]";
		else
			os << "4th kind on [0, lambdaMax]";
		if (_lambdaMaxEstimate > 0)
		{
			auto precision = os.precision(3);
			os << ", lambdaMax = " << EigenvalueSafetyFactor << " x " << _lambdaMaxEstimate;
			os.precision(precision);
		}
		os << ")";
	}

	// Lanczos estimate of the largest eigenvalue of D^-1 A (without the safety factor)
	double LambdaMaxEstimate() const
	{
		return _lambdaMaxEstimate;
	}

	void Setup(const SparseMatrix& A) override
	{
		IterativeSolver::Setup(A);
		_blockMatrix = nullptr;
		SetupBlocks(A);
	}

	// On several processes: the diagonal blocks of the rows of the process (in Diag), and the estimate of lambdaMax of
	// the whole matrix
	void Setup(const DistributedMatrix& A) override
	{
		SetDistributedMatrix(A);
		SetupBlocks(A.Diag);
	}

	// Mixed precision (one process): the same with the matrix in single precision, for the iterations in single
	// precision (SolveInSinglePrecision()), the sums of the passes in double
	bool SupportsSinglePrecision() const override
	{
		return !this->Distributed;
	}

	void SetupSinglePrecision(const FloatSparseMatrix& A) override
	{
		_singleMatrix = &A;
		_singleBlockMatrix = nullptr;
		SetupBlocks(A);
	}

	// The matrix stored by blocks (block CSR, the blocks of the block size of the smoother: of size > 1), in double or
	// in single precision (one process): the same iterations, faster passes (PassBlocks())
	bool SupportsBlockCSR() const override
	{
		return !this->Distributed;
	}

	void SetupBlockCSR(const BlockSparseMatrix<double>& A) override
	{
		this->Matrix = nullptr;
		_blockMatrix = &A;
		SetupBlocks(A);
	}

	void SetupBlockCSR(const BlockSparseMatrix<float>& A) override
	{
		_singleMatrix = nullptr;
		_singleBlockMatrix = &A;
		SetupBlocks(A);
	}

private:
	// The setup with A in precision T: the inverse diagonal (blocks of size 1) or the factorizations of the blocks
	// (computed in double, stored in T), the estimate of lambdaMax (Lanczos in T, the dot products in double)
	template <typename T>
	void SetupBlocks(const SparseMatrixOf<T>& A)
	{
		if (!A.isCompressed())
			Utils::FatalError("The Chebyshev smoothers require a compressed matrix (see Eigen's makeCompressed()).");

		if (_blockSize == 1)
		{
			VectorOf<T>& invDiag = In<T>().InvDiag;
			invDiag = VectorOf<T>(A.rows());
			bool positiveDiagonal = true;
			#pragma omp parallel for reduction(&&:positiveDiagonal)
			for (BigNumber i = 0; i < A.rows(); i++)
			{
				double a_ii = 0;
				for (typename SparseMatrixOf<T>::InnerIterator it(A, i); it; ++it)
				{
					if (it.col() == i)
						a_ii = it.value();
				}
				positiveDiagonal = positiveDiagonal && a_ii > 0;
				invDiag[i] = (T)(1 / a_ii);
			}
			if (!Communicator::Current().All(positiveDiagonal))
				Utils::FatalError("The Chebyshev smoothers require a symmetric positive definite matrix: non-positive diagonal coefficient.");
		}
		else
			_diagBlockSolver.Setup(A, _blockSize, DiagBlockSolveMethod::Inverse);
		EstimateLambdaMaxAndCost<T>();
	}

	template <typename T>
	void SetupBlocks(const BlockSparseMatrix<T>& A)
	{
		if (_blockSize == 1 || A.BlockSize != _blockSize)
			Utils::FatalError("The block CSR storage of a Chebyshev smoother requires blocks of its block size (" + to_string(_blockSize) + ", > 1), got " + to_string(A.BlockSize) + ".");
		_diagBlockSolver.Setup(A, _blockSize, DiagBlockSolveMethod::Inverse);
		EstimateLambdaMaxAndCost<T>();
	}

	template <typename T>
	void EstimateLambdaMaxAndCost()
	{
		_lambdaMaxEstimate = EstimateLambdaMax<T>();
		BigNumber n = Rows<T>();
		Flops lanczosIteration = Cost::SpMatVec(NonZeros<T>()) + 2 * _blockSize * n + 10 * n;
		this->SetupComputationalWork = (n / _blockSize * 2.0/3.0*pow(_blockSize, 3) + EigenvalueEstimateIterations * lanczosIteration) * 1e-6;
	}

public:

	// The diagonal blocks and the estimate of lambdaMax of another Chebyshev smoother set up for A with the same block
	// size, whatever its kind
	bool SetupFrom(const SparseMatrix& A, const IterativeSolver& other) override
	{
		const BlockChebyshev* o = dynamic_cast<const BlockChebyshev*>(&other);
		if (!o || o->Matrix != &A || o->Distributed || o->_blockSize != _blockSize || this->ComputeExactSolution)
			return false;
		IterativeSolver::Setup(A);
		_blockMatrix = nullptr;
		CopySetup(*o);
		return true;
	}

	bool SetupSinglePrecisionFrom(const FloatSparseMatrix& A, const IterativeSolver& other) override
	{
		const BlockChebyshev* o = dynamic_cast<const BlockChebyshev*>(&other);
		if (!o || o->_singleMatrix != &A || o->_blockSize != _blockSize || this->ComputeExactSolution)
			return false;
		_singleMatrix = &A;
		_singleBlockMatrix = nullptr;
		CopySetup(*o);
		return true;
	}

	bool SetupBlockCSRFrom(const BlockSparseMatrix<double>& A, const IterativeSolver& other) override
	{
		const BlockChebyshev* o = dynamic_cast<const BlockChebyshev*>(&other);
		if (!o || o->_blockMatrix != &A || o->_blockSize != _blockSize || this->ComputeExactSolution)
			return false;
		this->Matrix = nullptr;
		_blockMatrix = &A;
		CopySetup(*o);
		return true;
	}

	bool SetupBlockCSRFrom(const BlockSparseMatrix<float>& A, const IterativeSolver& other) override
	{
		const BlockChebyshev* o = dynamic_cast<const BlockChebyshev*>(&other);
		if (!o || o->_singleBlockMatrix != &A || o->_blockSize != _blockSize || this->ComputeExactSolution)
			return false;
		_singleMatrix = nullptr;
		_singleBlockMatrix = &A;
		CopySetup(*o);
		return true;
	}

	bool SetupFrom(const DistributedMatrix& A, const IterativeSolver& other) override
	{
		const BlockChebyshev* o = dynamic_cast<const BlockChebyshev*>(&other);
		if (!o || o->Distributed != &A || o->_blockSize != _blockSize || this->ComputeExactSolution)
			return false;
		SetDistributedMatrix(A);
		CopySetup(*o);
		return true;
	}

private:
	void CopySetup(const BlockChebyshev& o)
	{
		_diagBlockSolver = o._diagBlockSolver;
		_double.InvDiag = o._double.InvDiag;
		_single.InvDiag = o._single.InvDiag;
		_lambdaMaxEstimate = o._lambdaMaxEstimate;
		this->SetupComputationalWork = o.SetupComputationalWork;
	}

public:

	bool CanOptimizeResidualComputation() override
	{
		return true;
	}

	// m = MaxIterations iterations (the degree of the polynomial), with
	//   r_0 = b - A x_0,  d_0 = c_0 D^-1 r_0,
	//   x_k = x_{k-1} + d_{k-1},  r_k = r_{k-1} - A d_{k-1},  d_k = a_k d_{k-1} + c_k D^-1 r_k  (k = 1..m-1),
	//   x_m = x_{m-1} + d_{m-1}  (and r_m = r_{m-1} - A d_{m-1} for the residual or Ax).
	// Each line is one pass over the block rows. The residual is that of the recurrence: b - A x up to rounding.
	void Solve(const Vector& b, Vector& x, bool xEquals0, bool computeResidual, bool computeAx) override
	{
		if (this->MaxIterations == 0)
		{
			IterativeSolver::Solve(b, x, xEquals0, computeResidual, computeAx);
			return;
		}
		assert(this->StoppingCrit == StoppingCriteria::MaxIterations); // a smoother only: a fixed polynomial

		this->Residual = Vector();
		Flops flops = Iterate(b, x, xEquals0, computeResidual ? &this->Residual : nullptr, computeAx ? &this->Ax : nullptr);
		this->IterationCount = this->MaxIterations;
		this->SolvingComputationalWork = flops * 1e-6;
	}

	void SolveInSinglePrecision(const FloatVector& b, FloatVector& x, bool xEquals0, FloatVector* residual, FloatVector* Ax) override
	{
		assert(this->MaxIterations > 0 && (_singleMatrix || _singleBlockMatrix));
		Flops flops = Iterate(b, x, xEquals0, residual, Ax);
		this->IterationCount = this->MaxIterations;
		this->SolvingComputationalWork = flops * 1e-6;
	}

private:
	// The m = MaxIterations iterations in precision T (Solve()); the residual (that of the recurrence) and A x in residual
	// and Ax if not null. Returns the flops.
	template <typename T>
	Flops Iterate(const VectorOf<T>& b, VectorOf<T>& x, bool xEquals0, VectorOf<T>* residual, VectorOf<T>* Ax)
	{
		Precision<T>& W = In<T>();
		BigNumber n = Rows<T>();
		int m = this->MaxIterations;
		if (xEquals0 && x.rows() != n)
			x.resize(n);

		vector<double> a, c;
		Coefficients(m, a, c);

		W.d.resize(n);
		W.dNext.resize(n);
		W.r.resize(n);
		Flops matVec = Cost::SpMatVec(NonZeros<T>());
		Flops update = (2 * _blockSize + 3) * n; // D^-1 r and the linear combination of d
		Flops flops = 0;

		const T* r = b.data(); // r_k (b while x = 0)
		if (xEquals0)
		{
			Pass<T>(nullptr, b.data(), nullptr, nullptr, W.d.data(), 0, c[0], nullptr, false);                         flops += update;
		}
		else
		{
			Pass<T>(x.data(), b.data(), W.r.data(), nullptr, W.d.data(), 0, c[0], nullptr, false);                     flops += matVec + update;
			r = W.r.data();
		}
		for (int k = 1; k < m; k++)
		{
			Pass<T>(W.d.data(), r, W.r.data(), W.d.data(), W.dNext.data(), a[k], c[k], x.data(), xEquals0);            flops += matVec + update + n;
			r = W.r.data();
			xEquals0 = false;
			W.d.swap(W.dNext);
		}
		bool computeResidual = residual || Ax;
		Pass<T>(computeResidual ? W.d.data() : nullptr, r, W.r.data(), W.d.data(), nullptr, 0, 0, x.data(), xEquals0); flops += (computeResidual ? matVec : 0) + n;

		if (Ax)
		{
			*Ax = b - W.r;                                                                                        flops += n;
		}
		if (residual)
			*residual = std::move(W.r);
		return flops;
	}

private:
	// The coefficients a_k (a_0 = 0, unused) and c_k of the m iterations, from lambdaMax = EigenvalueSafetyFactor x the
	// estimate
	void Coefficients(int m, vector<double>& a, vector<double>& c) const
	{
		double lambdaMax = EigenvalueSafetyFactor * _lambdaMaxEstimate;
		a.assign(m, 0);
		c.assign(m, 0);
		if (_kind == ChebyshevKind::First)
		{
			// Saad, Iterative Methods for Sparse Linear Systems, 2nd ed., 2003, Algorithm 12.1, on [lambdaMin, lambdaMax]:
			// error polynomial T_m((theta - lambda) / delta) / T_m(theta / delta)
			double lambdaMin = lambdaMax / _ratio;
			double theta = (lambdaMax + lambdaMin) / 2;
			double delta = (lambdaMax - lambdaMin) / 2;
			double sigma = theta / delta;
			double rho = 1 / sigma;
			c[0] = 1 / theta;
			for (int k = 1; k < m; k++)
			{
				double rhoNext = 1 / (2 * sigma - rho);
				a[k] = rhoNext * rho;
				c[k] = 2 * rhoNext / delta;
				rho = rhoNext;
			}
		}
		else
		{
			// Lottes 2023: error polynomial W_m(1 - 2 lambda / lambdaMax) / (2m + 1), W_m the Chebyshev polynomial of the
			// fourth kind
			c[0] = 4 / (3 * lambdaMax);
			for (int k = 1; k < m; k++)
			{
				a[k] = (2.0 * k - 1) / (2.0 * k + 3);
				c[k] = (8.0 * k + 4) / ((2.0 * k + 3) * lambdaMax);
			}
		}
	}

	// One pass over the block rows I, independent of each other:
	//   r_I = rIn_I - (A v)_I               if v (otherwise r = rIn, and r is not written),
	//   dOut_I = a dIn_I + c D_I^-1 r_I     if dOut (a dIn_I omitted if dIn is null),
	//   x_I = dIn_I (xEquals0) or x_I += dIn_I   if x.
	// rIn may be r (in place). v must be neither r nor dOut nor x.
	// On several processes, A v is Diag v + Offd vGhosts, with the values of v of the ghosts of the column map (halo
	// exchange).
	// In precision T: the matrix and the vectors in T, the sums in double (in single precision: rounded when stored).
	template <typename T>
	void Pass(const T* v, const T* rIn, T* r, const T* dIn, T* dOut, double a, double c, T* x, bool xEquals0)
	{
		if (const BlockSparseMatrix<T>* blocks = BlockMatrixOf<T>())
		{
			DispatchBlockSize(blocks->BlockSize, [&](auto size)
			{
				PassBlocks<T, decltype(size)::value>(*blocks, v, rIn, r, dIn, dOut, a, c, x, xEquals0);
			});
			return;
		}
		const SparseMatrixOf<T>& A = MatrixOf<T>();
		const SparseMatrixIndex* outer = A.outerIndexPtr();
		const SparseMatrixIndex* col = A.innerIndexPtr();
		const T* val = A.valuePtr();
		BigNumber nb = A.rows() / _blockSize;
		// Off-process columns (none on one process; the single precision is on one process only)
		Vector vGhosts;
		const SparseMatrixIndex* offdOuter = nullptr;
		const SparseMatrixIndex* offdCol = nullptr;
		const double* offdVal = nullptr;
		if constexpr (std::is_same_v<T, double>)
		{
			if (this->Distributed && v)
			{
				vGhosts = this->Distributed->GhostValues(v);
				const SparseMatrix& offd = this->Distributed->Offd;
				if (offd.nonZeros() > 0)
				{
					offdOuter = offd.outerIndexPtr();
					offdCol = offd.innerIndexPtr();
					offdVal = offd.valuePtr();
				}
			}
		}
		const double* vg = vGhosts.data();

		if (_blockSize == 1)
		{
			const T* invDiag = In<T>().InvDiag.data();
			#pragma omp parallel for if (A.nonZeros() > 20000)
			for (BigNumber i = 0; i < nb; i++)
			{
				if (v || dOut)
				{
					double r_i = rIn[i];
					if (v)
					{
						for (SparseMatrixIndex p = outer[i]; p < outer[i + 1]; p++)
							r_i -= (double)val[p] * v[col[p]];
						if (offdOuter)
							for (SparseMatrixIndex p = offdOuter[i]; p < offdOuter[i + 1]; p++)
								r_i -= offdVal[p] * vg[offdCol[p]];
						r[i] = (T)r_i;
					}
					if (dOut)
						dOut[i] = (T)((dIn ? a * dIn[i] : 0) + c * (invDiag[i] * r_i));
				}
				if (x)
					x[i] = xEquals0 ? dIn[i] : x[i] + dIn[i];
			}
		}
		else
		{
			#pragma omp parallel if (A.nonZeros() > 20000)
			{
				vector<T> r_I(_blockSize); // work vectors, allocated once per thread
				vector<double> z_I(_blockSize);
				#pragma omp for
				for (BigNumber I = 0; I < nb; I++)
				{
					BigNumber i0 = I * _blockSize;
					for (int k = 0; k < _blockSize && (v || dOut); k++)
					{
						BigNumber i = i0 + k;
						double r_i = rIn[i];
						if (v)
						{
							for (SparseMatrixIndex p = outer[i]; p < outer[i + 1]; p++)
								r_i -= (double)val[p] * v[col[p]];
							if (offdOuter)
								for (SparseMatrixIndex p = offdOuter[i]; p < offdOuter[i + 1]; p++)
									r_i -= offdVal[p] * vg[offdCol[p]];
							r[i] = (T)r_i;
						}
						r_I[k] = (T)r_i;
					}
					if (dOut)
					{
						ApplyInverseBlock<T, 0>(I, r_I.data(), z_I.data());
						for (int k = 0; k < _blockSize; k++)
							dOut[i0 + k] = (T)((dIn ? a * dIn[i0 + k] : 0) + c * z_I[k]);
					}
					if (x)
					{
						for (int k = 0; k < _blockSize; k++)
							x[i0 + k] = xEquals0 ? dIn[i0 + k] : x[i0 + k] + dIn[i0 + k];
					}
				}
			}
		}
	}

	// Pass() on the matrix stored by blocks: the same operations by block rows, the same results; B the block size if
	// it is known at compile time (DispatchBlockSize()), 0 otherwise
	template <typename T, int B>
	void PassBlocks(const BlockSparseMatrix<T>& A, const T* v, const T* rIn, T* r, const T* dIn, T* dOut, double a, double c, T* x, bool xEquals0)
	{
		const int bs = B > 0 ? B : _blockSize;
		BigNumber nb = A.NBlocks();
		#pragma omp parallel if (A.NonZeros() > 20000)
		{
			T r_I[B > 0 ? B : BlockSparseMatrix<T>::MaxBlockSize];
			double acc[B > 0 ? B : BlockSparseMatrix<T>::MaxBlockSize], z_I[B > 0 ? B : BlockSparseMatrix<T>::MaxBlockSize];
			#pragma omp for
			for (BigNumber I = 0; I < nb; I++)
			{
				BigNumber i0 = I * bs;
				if (v)
				{
					A.template RowResidual<B>(I, v, rIn, acc);
					for (int k = 0; k < bs; k++)
					{
						r[i0 + k] = (T)acc[k];
						r_I[k] = (T)acc[k];
					}
				}
				else if (dOut)
				{
					for (int k = 0; k < bs; k++)
						r_I[k] = rIn[i0 + k];
				}
				if (dOut)
				{
					ApplyInverseBlock<T, B>(I, r_I, z_I);
					for (int k = 0; k < bs; k++)
						dOut[i0 + k] = (T)((dIn ? a * dIn[i0 + k] : 0) + c * z_I[k]);
				}
				if (x)
				{
					for (int k = 0; k < bs; k++)
						x[i0 + k] = xEquals0 ? dIn[i0 + k] : x[i0 + k] + dIn[i0 + k];
				}
			}
		}
	}

	// z_I = D_I^-1 r_I by the inverse of the diagonal block I (blocks of size > 1), summed in double in the order of the
	// columns of the inverse; B: the block size if it is known at compile time, 0 otherwise
	template <typename T, int B>
	inline void ApplyInverseBlock(BigNumber I, const T* r_I, double* z_I) const
	{
		const int bs = B > 0 ? B : _blockSize;
		const T* inverse = _diagBlockSolver.template Block<T>(I);
		for (int k = 0; k < bs; k++)
			z_I[k] = 0;
		for (int l = 0; l < bs; l++)
		{
			double r_l = r_I[l];
			const T* column = inverse + l * bs;
			for (int k = 0; k < bs; k++)
				z_I[k] += (double)column[k] * r_l;
		}
	}

	// z = D^-1 r, in precision T
	template <typename T>
	void ApplyInverseDiagonal(const VectorOf<T>& r, VectorOf<T>& z) const
	{
		if (_blockSize == 1)
		{
			z = In<T>().InvDiag.cwiseProduct(r);
			return;
		}
		BigNumber nb = Rows<T>() / _blockSize;
		#pragma omp parallel if (NonZeros<T>() > 20000)
		{
			vector<double> z_I(_blockSize);
			#pragma omp for
			for (BigNumber I = 0; I < nb; I++)
			{
				ApplyInverseBlock<T, 0>(I, r.data() + I * _blockSize, z_I.data());
				for (int k = 0; k < _blockSize; k++)
					z[I * _blockSize + k] = (T)z_I[k];
			}
		}
	}

	// Largest eigenvalue of D^-1 A, by Lanczos: D^-1 A is similar to the symmetric D^-1/2 A D^-1/2, whose Lanczos
	// tridiagonal matrix is given by the coefficients of the conjugate gradient preconditioned by D (Saad, Iterative
	// Methods for Sparse Linear Systems, 2nd ed., 2003, §6.7.3). Its largest eigenvalue approaches lambdaMax from below.
	// Pseudo-random start vector, a function of the index only: the estimate is the same in every run, whatever the
	// number of threads. In precision T (the dot products in double).
	template <typename T>
	double EstimateLambdaMax() const
	{
		BigNumber n = Rows<T>();
		// On several processes: the global numbers of the rows, the products with the distributed matrix
		BigNumber firstRow = this->Distributed ? this->Distributed->RowMap->Begin * this->Distributed->RowBlockSize : 0;
		VectorOf<T> r(n), z(n), p, Ap(n);
		for (BigNumber i = 0; i < n; i++)
			r[i] = PseudoRandom(firstRow + i);
		ApplyInverseDiagonal(r, z);
		double rz = Parallel::Dot(r, z);
		if (!(rz > 0))
			Utils::FatalError("The Chebyshev smoothers require a symmetric positive definite matrix: its diagonal blocks are not positive definite.");
		double rz0 = rz;
		p = z;
		vector<double> alpha, beta;
		for (int j = 0; j < EigenvalueEstimateIterations; j++)
		{
			if (const BlockSparseMatrix<T>* blocks = BlockMatrixOf<T>())
				blocks->Multiply(p.data(), Ap.data());
			else if constexpr (std::is_same_v<T, double>)
			{
				if (this->Distributed)
					Ap = *this->Distributed * p;
				else
					Ap.noalias() = *this->Matrix * p;
			}
			else
				Ap.noalias() = *_singleMatrix * p;
			double pAp = Parallel::Dot(p, Ap);
			if (!(pAp > 0))
				Utils::FatalError("The Chebyshev smoothers require a symmetric positive definite matrix: p^T A p = " + to_string(pAp) + " in the estimation of its largest eigenvalue.");
			alpha.push_back(rz / pAp);
			if (j == EigenvalueEstimateIterations - 1)
				break;
			r -= T(alpha.back()) * Ap;
			ApplyInverseDiagonal(r, z);
			double rzNext = Parallel::Dot(r, z);
			if (!(rzNext > 1e-28 * rz0)) // invariant Krylov subspace: the eigenvalues of the tridiagonal matrix are exact
				break;
			beta.push_back(rzNext / rz);
			p = z + T(beta.back()) * p;
			rz = rzNext;
		}

		int m = alpha.size();
		DenseMatrix tridiagonal = DenseMatrix::Zero(m, m);
		for (int j = 0; j < m; j++)
		{
			tridiagonal(j, j) = 1 / alpha[j] + (j > 0 ? beta[j - 1] / alpha[j - 1] : 0);
			if (j > 0)
				tridiagonal(j, j - 1) = tridiagonal(j - 1, j) = sqrt(beta[j - 1]) / alpha[j - 1];
		}
		Eigen::SelfAdjointEigenSolver<DenseMatrix> eigenSolver(tridiagonal, Eigen::EigenvaluesOnly);
		return eigenSolver.eigenvalues().maxCoeff();
	}

	// In [-1, 1), from the index (SplitMix64)
	static double PseudoRandom(uint64_t i)
	{
		uint64_t z = i + 0x9E3779B97F4A7C15ull;
		z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
		z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
		z ^= z >> 31;
		return (z >> 11) * 0x1.0p-52 - 1;
	}
};
