#pragma once
#include "../Solver.h"

// The solver of a singular symmetric matrix whose kernel is spanned by the interpolation of 1 (the coarsest level of
// U-AMG for a pure Neumann problem, see UncondensedAMG::SingularMatrix): the kernel z is c on the first DoF of each
// block and 0 on the others (1 = c phi_0). The matrix is pinned on the DoF of largest |z_i| (its row and column replaced
// by those of the identity): regular, factorized by the solver given. Each solve projects the right-hand side onto the
// image of the matrix (orthogonal to z), zeroes the pinned DoF, solves, and returns the solution orthogonal to z.
// Without the pin, the factorization meets a pivot of the order of the rounding errors, and the solutions have a
// kernel component that ruins the Krylov method. The kernel is known: no threshold to detect the null pivot (as MUMPS's
// ICNTL(24) would need), and the same for all the factorizations, on one process or several.
class SingularSolver : public Solver
{
private:
	unique_ptr<Solver> _pinnedSolver;
	Vector _interpOfOne; // one coefficient per block; empty: 1
	int _blockSize = 1;
	Vector _z;           // the kernel, normalized
	Eigen::Index _pin = -1; // the pinned DoF, local; -1 on the processes that don't own it
	SparseMatrix _pinnedA;
	DistributedMatrix _pinnedDistributedA;
public:
	// pinnedSolver: owned
	SingularSolver(Solver* pinnedSolver)
		: _pinnedSolver(pinnedSolver)
	{}

	// The kernel on the n DoFs of the process, normalized (on all the processes of Communicator::Current()): the
	// interpolation of 1 c on the first DoF of each block of blockSize DoFs (1 if interpOfOne is empty; its first n /
	// blockSize coefficients: those of the owned blocks), 0 on the others
	static Vector Kernel(const Vector& interpOfOne, int blockSize, Eigen::Index n)
	{
		Vector z = Vector::Zero(n);
		for (Eigen::Index i = 0; i < n; i += blockSize)
			z[i] = interpOfOne.rows() > 0 ? interpOfOne[i / blockSize] : 1;
		return z / Parallel::Norm(z);
	}

	// The interpolation of 1 on the blocks of blockSize DoFs of the matrix of the next setup (empty: 1 everywhere)
	void SetKernel(const Vector& interpOfOne, int blockSize)
	{
		_interpOfOne = interpOfOne;
		_blockSize = blockSize;
	}

	void Serialize(ostream& os) const override
	{
		os << *_pinnedSolver << " (singular matrix: one DoF pinned, projections onto the image)";
	}

	void Setup(const SparseMatrix& A) override
	{
		Solver::Setup(A);
		_z = Kernel(_interpOfOne, _blockSize, A.rows());
		_z.cwiseAbs().maxCoeff(&_pin);
		_pinnedA = A;
		Pin(_pinnedA, _pin, _pin);
		_pinnedSolver->Setup(_pinnedA);
		this->SetupComputationalWork = _pinnedSolver->SetupComputationalWork;
	}

	// On several processes: the pinned DoF is that of largest |z_i| on the smallest rank that has it. Its owner pins its
	// row and its column of Diag; the processes that have it as a ghost zero its column of Offd. Collective.
	void Setup(const DistributedMatrix& A) override
	{
		this->Matrix = nullptr;
		this->Distributed = &A;
		const Communicator& comm = A.Comm();
		_z = Kernel(_interpOfOne, _blockSize, A.LocalRows());
		Eigen::Index localPin = -1;
		double localMax = _z.rows() > 0 ? _z.cwiseAbs().maxCoeff(&localPin) : -1;
		int owner = comm.Min(localMax == comm.Max(localMax) ? comm.Rank() : comm.Size());
		_pin = comm.Rank() == owner ? localPin : -1;
		BigNumber pinGlobal = comm.Sum<BigNumber>(_pin >= 0 ? A.RowMap->Begin * A.RowBlockSize + _pin : 0);

		_pinnedDistributedA = A;
		BigNumber ghostColumn = -1;
		const vector<BigNumber>& ghosts = A.ColMap->Ghosts;
		auto ghost = find(ghosts.begin(), ghosts.end(), pinGlobal / A.ColBlockSize);
		if (ghost != ghosts.end())
			ghostColumn = (ghost - ghosts.begin()) * A.ColBlockSize + pinGlobal % A.ColBlockSize;
		Pin(_pinnedDistributedA.Diag, _pin, _pin);
		Pin(_pinnedDistributedA.Offd, _pin, ghostColumn);
		_pinnedSolver->Setup(_pinnedDistributedA);
		this->SetupComputationalWork = _pinnedSolver->SetupComputationalWork;
	}

	Vector Solve(const Vector& b) override
	{
		Vector rhs = b - Parallel::Dot(b, _z) * _z;
		if (_pin >= 0)
			rhs[_pin] = 0;
		Vector x = _pinnedSolver->Solve(rhs);
		x -= Parallel::Dot(x, _z) * _z;
		this->SolvingComputationalWork = _pinnedSolver->SolvingComputationalWork; // + the projections, negligible
		return x;
	}

private:
	// Zeroes the row and the column of M (-1: none), but their common coefficient, set to 1 (the diagonal of the pinned DoF)
	static void Pin(SparseMatrix& M, Eigen::Index row, BigNumber col)
	{
		if (row < 0 && col < 0)
			return;
		for (Eigen::Index i = 0; i < M.rows(); ++i)
		{
			for (SparseMatrix::InnerIterator it(M, i); it; ++it)
			{
				if (i == row || it.col() == col)
					it.valueRef() = i == row && it.col() == col ? 1 : 0;
			}
		}
	}
};
