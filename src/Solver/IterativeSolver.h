#pragma once
#include "Solver.h"
#include "IterationResult.h"
#include "../Utils/BlockSparseMatrix.h"
using namespace std;

class IterativeSolver : public Solver
{
protected:
	// Computes ExactSolution, the exact solution of the linear system, against which the algebraic
	// error of the iterates is computed (see ComputeExactSolution).
	// Column-major type: see EigenSparseLU.
	Eigen::SparseLU<ColMajorSparseMatrix> _directSolver;

	// The matrix of the smoothers that support the single precision and the block CSR storage (SetupSinglePrecision(),
	// SetupBlockCSR()), besides Matrix (CSR in double): in CSR in single precision, in block CSR in double or in single
	// precision; on several processes, the distributed matrix in single precision whose block Diag is _singleMatrix or
	// _singleBlockMatrix (null on one process)
	const FloatSparseMatrix* _singleMatrix = nullptr;
	const DistributedMatrixOf<float>* _singleDistributed = nullptr;
	const BlockSparseMatrix<double>* _blockMatrix = nullptr;
	const BlockSparseMatrix<float>* _singleBlockMatrix = nullptr;

	// The matrix in CSR, the distributed matrix (null on one process) and the matrix in block CSR (null in CSR) in
	// precision T
	template <typename T>
	const SparseMatrixOf<T>& MatrixOf() const
	{
		if constexpr (std::is_same_v<T, float>)
			return *_singleMatrix;
		else
			return *this->Matrix;
	}

	template <typename T>
	const DistributedMatrixOf<T>* DistributedOf() const
	{
		if constexpr (std::is_same_v<T, float>)
			return _singleDistributed;
		else
			return this->Distributed;
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

	// The matrix of the setup in one storage, the others null (the Gauss-Seidel smoothers)
	void SetStorage(const SparseMatrix* A, const FloatSparseMatrix* singleA, const BlockSparseMatrix<double>* blocks, const BlockSparseMatrix<float>* singleBlocks, const DistributedMatrix* distributed, const DistributedMatrixOf<float>* singleDistributed)
	{
		this->Matrix = A;
		_singleMatrix = singleA;
		_blockMatrix = blocks;
		_singleBlockMatrix = singleBlocks;
		this->Distributed = distributed;
		_singleDistributed = singleDistributed;
	}

	bool SameStorage(const IterativeSolver& o) const
	{
		return o.Matrix == this->Matrix && o._singleMatrix == _singleMatrix && o._blockMatrix == _blockMatrix && o._singleBlockMatrix == _singleBlockMatrix
			&& o.Distributed == this->Distributed && o._singleDistributed == _singleDistributed;
	}

	// Whether the matrix of the setup is in single precision
	bool InSinglePrecision() const
	{
		return _singleMatrix || _singleBlockMatrix;
	}

	// The numbers of rows and of non-zeros of the matrix of the setup (SetStorage())
	BigNumber MatrixRows() const
	{
		return InSinglePrecision() ? Rows<float>() : Rows<double>();
	}

	BigNumber MatrixNonZeros() const
	{
		return InSinglePrecision() ? NonZeros<float>() : NonZeros<double>();
	}
public:
	Vector ExactSolution;
	double Tolerance = 1e-3;
	double StagnationConvRate = 0.90;
	int MaxIterations = 200;
	int IterationCount = 0;
	bool PrintIterationResults = true;
	bool ComputeExactSolution = false;
	function<void(IterationResult&, const Vector&)> OnNewSolution = nullptr;
	function<void(const IterationResult&)> OnIterationEnd = nullptr;
	StoppingCriteria StoppingCrit = StoppingCriteria::NormalizedResidual;
	Vector Residual;
	Vector Ax;
	IterationResult LastIterationResult;

	IterativeSolver() : Solver() {}

	using Solver::Setup; // the distributed one, among the overloads below

	// Setup for A from that of other, if other is a solver of the same kind set up for A (e.g. the pre- and
	// post-smoothers of a multigrid level): the result of Setup(A), at the cost of a copy. Returns false (and does
	// nothing) if it can't.
	virtual bool SetupFrom(const SparseMatrix& A, const IterativeSolver& other)
	{
		return false;
	}

	void Setup(const SparseMatrix& A) override
	{
		Solver::Setup(A);
		if (this->ComputeExactSolution)
		{
			// Copy: A is row-major, SparseLU needs column-major. Freed at the end of this block.
			ColMajorSparseMatrix colMajorA = A;
			this->_directSolver.compute(colMajorA);
			// A singular matrix: no exact solution (solving from a failed factorization reads uninitialized memory)
			if (this->_directSolver.info() != Eigen::Success)
			{
				Utils::Warning("The direct solver could not factorize the matrix: the algebraic error is not computed.");
				this->ComputeExactSolution = false;
			}
		}
	}

	// Same as SetupFrom() on several processes
	virtual bool SetupFrom(const DistributedMatrix& A, const IterativeSolver& other)
	{
		return false;
	}

	// Setup on several processes from the blocks of the uncondensed system of a hybrid discretization (see the
	// sequential one above), distributed by rows: A (faces x faces) the rows of the faces of the process;
	// A_T_T and A_T_F the rows of its cells, the columns of A_T_F numbered in faceMap (its faces, then the other faces of
	// its cells, the ghosts); A_F_T empty if the matrix is symmetric, otherwise the rows of all the faces of faceMap and
	// the columns of the cells of the process (DistributedMatrix::RowsWithGhosts()); the interpolations of 1 on its cells
	// and faces.
	virtual void Setup(const DistributedMatrix& A, const SparseMatrix& A_T_T, const SparseMatrix& A_T_F, const SparseMatrix& A_F_T, shared_ptr<const IndexMap> faceMap, const Vector& cellInterpOfOne, const Vector& faceInterpOfOne)
	{
		Utils::FatalError("This solver is not compatible with matrix blocks arising from a hybrid discretization on several processes.");
	}

	// b - A x, with the matrix of the setup (distributed or not)
	Vector ResidualOf(const Vector& b, const Vector& x) const
	{
		if (this->Distributed)
			return b - *this->Distributed * x;
		return b - *this->Matrix * x;
	}

	// A x, with the matrix of the setup (distributed or not)
	Vector MatVec(const Vector& x) const
	{
		if (this->Distributed)
			return *this->Distributed * x;
		return *this->Matrix * x;
	}

protected:
	// On several processes, the right-hand side of a sweep of a relaxation method (Gauss-Seidel, Jacobi) on the rows of
	// the process (Diag): b - Offd x, with the values of x of the other processes from the start of the sweep (Jacobi
	// between the processes, Gauss-Seidel within each one, as the hybrid smoothers between threads); b on one process, or
	// from x = 0. With l1 terms (l1 Gauss-Seidel, BlockGaussSeidel::L1Theta: not empty), + l1Terms .* x. buffer: its
	// storage. distributed: the distributed matrix (null on one process), in precision T (float: the l1 terms added in
	// double).
	// Collective on the halo.
	template <typename T>
	const VectorOf<T>& SweepRhs(const VectorOf<T>& b, const VectorOf<T>& x, bool xEquals0, VectorOf<T>& buffer, const DistributedMatrixOf<T>* distributed, const Vector* l1Terms = nullptr) const
	{
		bool l1 = l1Terms && l1Terms->size() > 0;
		if ((!distributed && !l1) || xEquals0)
			return b;
		buffer = b;
		if (distributed)
		{
			VectorOf<T> ghosts = distributed->GhostValues(x);
			if (distributed->Offd.nonZeros() > 0)
				buffer.noalias() -= distributed->Offd * ghosts;
		}
		if (l1)
		{
			if constexpr (std::is_same_v<T, double>)
				buffer += l1Terms->cwiseProduct(x);
			else
				buffer = (buffer.template cast<double>() + l1Terms->cwiseProduct(x.template cast<double>())).template cast<T>();
		}
		return buffer;
	}

	// The part of the setup on several processes common to the iterative solvers: A distributed by rows (Matrix: its
	// block Diag). The vectors of the solves are then the rows of the process. Each solver that runs on several
	// processes overrides Setup(const DistributedMatrix&) and calls it.
	void SetDistributedMatrix(const DistributedMatrix& A)
	{
		if (this->ComputeExactSolution)
			Utils::FatalError("The exact solution of the linear system (algebraic error) is not computed on several processes.");
		Solver::Setup(A.Diag);
		this->Distributed = &A;
	}

public:

	// Setup from the blocks of the uncondensed system of a hybrid discretization (A being the condensed one), and the
	// interpolation of the function 1 on the polynomial bases of the cells (cellInterpOfOne, one coefficient per cell,
	// numbered as the blocks of rows of A_T_F) and of the faces of the system (faceInterpOfOne, one per face, as the
	// blocks of columns of A_T_F). The bases must be hierarchical with a constant first function phi_0: 1 = c phi_0 (see
	// UncondensedLevel::CoarsenMesh()).
	// A_F_T: empty (0 x 0) if the matrix is symmetric (A_F_T = A_T_F^T); otherwise, the matrix is non-symmetric, and all
	// the blocks are read in full.
	virtual void Setup(const SparseMatrix& A, const SparseMatrix& A_T_T, const SparseMatrix& A_T_F, const SparseMatrix& A_F_T, const SparseMatrix& A_F_F, const Vector& cellInterpOfOne, const Vector& faceInterpOfOne)
	{
		Utils::FatalError("This solver is not compatible with matrix blocks arising from a hybrid discretization.");
	}

	// Mixed precision: a solver usable in single precision (the smoothers of the levels of U-AMG in single precision,
	// Multigrid::MixedPrecision), set up with a matrix in single precision (instead of Setup()).
	virtual bool SupportsSinglePrecision() const { return false; }
	virtual void SetupSinglePrecision(const FloatSparseMatrix& A)
	{
		Utils::FatalError("This solver has no single precision version.");
	}
	// As SetupFrom()
	virtual bool SetupSinglePrecisionFrom(const FloatSparseMatrix& A, const IterativeSolver& other)
	{
		return false;
	}
	// On several processes: the matrix in single precision distributed by rows
	virtual void SetupSinglePrecision(const DistributedMatrixOf<float>& A)
	{
		Utils::FatalError("This solver has no single precision version on several processes.");
	}
	virtual bool SetupSinglePrecisionFrom(const DistributedMatrixOf<float>& A, const IterativeSolver& other)
	{
		return false;
	}

	// The matrix stored by blocks (block CSR, see BlockSparseMatrix), in double or in single precision (as
	// SetupSinglePrecision()): the solvers that support it (the smoothers of U-AMG), instead of Setup(). On several
	// processes, distributed: the matrix distributed by rows whose block Diag A stores (its Offd and ghosts; its Diag is
	// not read); null on one process.
	virtual bool SupportsBlockCSR() const { return false; }
	virtual void SetupBlockCSR(const BlockSparseMatrix<double>& A, const DistributedMatrix* distributed)
	{
		Utils::FatalError("This solver does not support the block CSR storage.");
	}
	virtual void SetupBlockCSR(const BlockSparseMatrix<float>& A, const DistributedMatrixOf<float>* distributed)
	{
		Utils::FatalError("This solver does not support the block CSR storage.");
	}
	// As SetupFrom()
	virtual bool SetupBlockCSRFrom(const BlockSparseMatrix<double>& A, const DistributedMatrix* distributed, const IterativeSolver& other)
	{
		return false;
	}
	virtual bool SetupBlockCSRFrom(const BlockSparseMatrix<float>& A, const DistributedMatrixOf<float>* distributed, const IterativeSolver& other)
	{
		return false;
	}
	// As Solve(b, x, xEquals0, computeResidual, computeAx), in single precision: the residual and A x in residual and Ax
	// if not null
	virtual void SolveInSinglePrecision(const FloatVector& b, FloatVector& x, bool xEquals0, FloatVector* residual, FloatVector* Ax)
	{
		Utils::FatalError("This solver has no single precision version.");
	}

	virtual bool CanOptimizeResidualComputation()
	{
		return false;
	}

	Vector Solve() override
	{
		Utils::FatalError("Solve() not implemented in this solver.");
		return Vector::Zero(0);
	}

	Vector Solve(const Vector& b) override
	{
		return Solve(b, "0");
	}

	virtual Vector Solve(const Vector& b, string initialGuessCode)
	{
		Vector x = InitialGuess(initialGuessCode, b.rows());
		Solve(b, x, initialGuessCode.compare("0") == 0);
		return x;
	}

	// The initial guess of the code (-initial-guess): "0", "1" or "rand"
	static Vector InitialGuess(const string& initialGuessCode, BigNumber size)
	{
		if (initialGuessCode.compare("0") == 0)
			return Vector::Zero(size);
		if (initialGuessCode.compare("1") == 0)
			return Vector::Ones(size);
		if (initialGuessCode.compare("rand") == 0)
			return Vector::Random(size);
		Utils::FatalError("Unknown initial guess code");
		return Vector();
	}

	virtual void Solve(const Vector& b, Vector& x, bool xEquals0)
	{
		bool computeResidual = false;
		bool computeAx = false;
		Solve(b, x, xEquals0, computeResidual, computeAx);
	}

	virtual void Solve(const Vector& b, Vector& x, bool xEquals0, bool computeResidual, bool computeAx)
	{
		DoBeforeSolving();

		const SparseMatrix& A = *this->Matrix;

		this->SolvingComputationalWork = 0;
		this->Residual = Vector();

		if (this->ComputeExactSolution)
			this->ExactSolution = this->_directSolver.solve(b);

		this->IterationCount = 0;
		
		IterationResult result = CreateFirstIterationResult(b, x);
		if (MaxIterations == 0)
		{
			if (computeResidual && !computeAx)
			{
				if (xEquals0)
					this->Residual = b;
				else
				{
					this->Residual = ResidualOf(b, x);                       result.AddWorkInFlops(Cost::DAXPY(A));
				}
			}
			else if (computeResidual && computeAx)
			{
				if (xEquals0)
				{
					this->Residual = b;
					this->Ax = Vector::Zero(b.rows());
				}
				else
				{
					this->Ax = MatVec(x);                                    result.AddWorkInFlops(Cost::MatVec(A));
					this->Residual = b - this->Ax;                           result.AddWorkInFlops(Cost::AddVec(b));
				}
			}
			else if (computeAx)
			{
				if (xEquals0)
					this->Ax = Vector::Zero(b.rows());
				else
				{
					this->Ax = MatVec(x);                                    result.AddWorkInFlops(Cost::MatVec(A));
				}
			}
			return;
		}

		if (StoppingCrit == StoppingCriteria::NormalizedResidual || StoppingCrit == StoppingCriteria::Stagnation)
		{
			if (xEquals0)
				result.SetResidualAsB();
			else
			{
				this->Residual = ResidualOf(b, x);                                   result.AddWorkInFlops(Cost::DAXPY(A));
				result.SetResidualNorm(Parallel::Norm(this->Residual));              result.AddWorkInFlops(Cost::Norm(b));
			}
		}

		if (this->PrintIterationResults)
			Utils::Log() << result << endl;


		while (!StoppingCriteriaReached(result))
		{
			if (!CanOptimizeResidualComputation())
			{
				result = ExecuteOneIteration(b, x, xEquals0, false, false, result);
				if (!result.IsResidualSet() && (StoppingCrit == StoppingCriteria::NormalizedResidual || StoppingCrit == StoppingCriteria::Stagnation))
				{
					this->Residual = ResidualOf(b, x);                                   result.AddWorkInFlops(Cost::DAXPY(A));
					result.SetResidualNorm(Parallel::Norm(this->Residual));              result.AddWorkInFlops(Cost::Norm(b));
				}                                      
			}
			else
			{
				if (StoppingCrit == StoppingCriteria::MaxIterations)
				{
					if ((!computeResidual && !computeAx) || this->IterationCount < this->MaxIterations - 1)
						result = ExecuteOneIteration(b, x, xEquals0, false, false, result); // No need to compute the residual until the last iteration
					else
						result = ExecuteOneIteration(b, x, xEquals0, computeResidual, computeAx, result);
				}
				else
				{
					result = ExecuteOneIteration(b, x, xEquals0, true, computeAx, result);
					assert(result.Residual.rows() > 0);
					result.SetResidualNorm(Parallel::Norm(result.Residual));                                       result.AddWorkInFlops(Cost::Norm(b));
				}
			}
			this->IterationCount++;

			if (this->PrintIterationResults)
				Utils::Log() << result << endl;
		}

		result.CopyInfoInto(LastIterationResult);

		if (this->PrintIterationResults)
			Utils::Log() << endl;

		this->SolvingComputationalWork = result.SolvingComputationalWork();
		if (computeResidual)
			this->Residual = std::move(result.Residual);
		if (computeAx)
			this->Ax = std::move(result.Ax);
	}

	virtual ~IterativeSolver() {}

protected:
	virtual void DoBeforeSolving()
	{}

	IterationResult CreateFirstIterationResult(const Vector& b, const Vector& x)
	{
		IterationResult result;
		if (this->MaxIterations == 0 || this->StoppingCrit == StoppingCriteria::MaxIterations)
			return result;
		if (this->Matrix)
			result.SetA(*this->Matrix);
		result.SetB(b);
		if (this->ComputeExactSolution)
			result.SetExactSolution(this->ExactSolution);
		if (this->OnNewSolution)
			result.OnNewSolution = this->OnNewSolution;
		result.SetX(x);
		if (this->StoppingCrit == StoppingCriteria::NormalizedResidual)
			result.SetTolerance(this->Tolerance);
		return result;
	}

	virtual IterationResult ExecuteOneIteration(const Vector& b, Vector& x, bool& xEquals0, bool computeResidual, bool computeAx, const IterationResult& oldResult) 
	{ 
		Utils::FatalError("The function Solve() or ExecuteOneIteration() must be overriden in the subclass.");
		return IterationResult(); // to avoid warning
	};


	bool StoppingCriteriaReached(const IterationResult& result)
	{
		if (this->OnIterationEnd)
			this->OnIterationEnd(result);

		return StoppingCriteriaMet(result);
	}

	// The stopping criteria alone, without calling OnIterationEnd
	bool StoppingCriteriaMet(const IterationResult& result) const
	{
		if (MaxIterations == 0)
			return true;
		if (IterationCount == 0)
			return false;
		if (IterationCount >= MaxIterations)
			return true;
		if (StoppingCrit == StoppingCriteria::NormalizedResidual || StoppingCrit == StoppingCriteria::Stagnation)
		{
			if (std::isinf(result.NormalizedResidualNorm)) // do not remove the prefix std:: or it can become ambiguous according to the compiler
				Utils::FatalError("divergence of the solver");
			if (std::isnan(result.NormalizedResidualNorm))
				Utils::FatalError("the residual is NaN.");
			if (StoppingCrit == StoppingCriteria::NormalizedResidual && result.NormalizedResidualNorm < this->Tolerance)
				return true;
			if (StoppingCrit == StoppingCriteria::Stagnation && result.AsymptoticConvRate >= this->StagnationConvRate && result.NormalizedResidualNorm <= result.PreviousNormalizedResidualNorm)
				return true;
		}
		return false;
	}
};