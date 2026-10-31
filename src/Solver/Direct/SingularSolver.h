#pragma once
#include "../Solver.h"

// The solver of a singular symmetric matrix whose kernel is spanned by the interpolation of 1 (the coarsest level of
// U-AMG for a pure Neumann problem, see UncondensedAMG::SingularMatrix): the kernel z is c on the first DoF of each
// block and 0 on the others (1 = c phi_0); with n components (UncondensedAMG::VectorComponents), n vectors z_k, c on the
// DoF k of each block (the constant of the component k), orthogonal (disjoint supports). The matrix is pinned on the DoF
// of largest |z_i| of each (its row and column replaced by those of the identity): regular (the kernel vectors restricted
// to the pinned DoFs are independent), factorized by the solver given. Each solve projects the right-hand side onto the
// image of the matrix (orthogonal to the kernel), zeroes the pinned DoFs, solves, and returns the solution orthogonal to
// the kernel. Without the pins, the factorization meets pivots of the order of the rounding errors, and the solutions
// have a kernel component that ruins the Krylov method. The kernel is known: no threshold to detect the null pivots (as
// MUMPS's ICNTL(24) would need), and the same for all the factorizations, on one process or several.
class SingularSolver : public Solver
{
private:
	unique_ptr<Solver> _pinnedSolver;
	Vector _interpOfOne; // one coefficient per block; empty: 1
	int _blockSize = 1;
	int _components = 1;
	vector<Vector> _z;            // the kernel, normalized vectors
	vector<Eigen::Index> _pins;   // the pinned DoFs (one per vector of the kernel), local: those of the process
	SparseMatrix _pinnedA;
	DistributedMatrix _pinnedDistributedA;
public:
	// pinnedSolver: owned
	SingularSolver(Solver* pinnedSolver)
		: _pinnedSolver(pinnedSolver)
	{}

	// The kernel on the n DoFs of the process, normalized vectors (on all the processes of Communicator::Current()): for
	// each of the components k, the interpolation of 1 c on the DoF k of each block of blockSize DoFs (1 if interpOfOne
	// is empty; its first n / blockSize coefficients: those of the owned blocks), 0 on the others
	static vector<Vector> Kernel(const Vector& interpOfOne, int blockSize, Eigen::Index n, int components)
	{
		vector<Vector> kernel;
		for (int k = 0; k < components; ++k)
		{
			Vector z = Vector::Zero(n);
			for (Eigen::Index i = k; i < n; i += blockSize)
				z[i] = interpOfOne.rows() > 0 ? interpOfOne[i / blockSize] : 1;
			kernel.push_back(z / Parallel::Norm(z));
		}
		return kernel;
	}

	// x projected orthogonally to the kernel (its vectors orthogonal: one after the other); returns the norm of the
	// component removed
	static double Project(Vector& x, const vector<Vector>& kernel)
	{
		double removed = 0;
		for (const Vector& z : kernel)
		{
			double component = Parallel::Dot(x, z);
			x -= component * z;
			removed = hypot(removed, component);
		}
		return removed;
	}

	// The interpolation of 1 on the blocks of blockSize DoFs of the matrix of the next setup (empty: 1 everywhere), and
	// the number of components, whose constants are the first DoFs of each block
	void SetKernel(const Vector& interpOfOne, int blockSize, int components)
	{
		_interpOfOne = interpOfOne;
		_blockSize = blockSize;
		_components = components;
	}

	void Serialize(ostream& os) const override
	{
		os << *_pinnedSolver << " (singular matrix: " << (_components == 1 ? "one DoF" : to_string(_components) + " DoFs") << " pinned, projections onto the image)";
	}

	void Setup(const SparseMatrix& A) override
	{
		Solver::Setup(A);
		_z = Kernel(_interpOfOne, _blockSize, A.rows(), _components);
		_pinnedA = A;
		_pins.clear();
		for (const Vector& z : _z)
		{
			Eigen::Index pin;
			z.cwiseAbs().maxCoeff(&pin);
			_pins.push_back(pin);
			Pin(_pinnedA, pin, pin);
		}
		_pinnedSolver->Setup(_pinnedA);
		this->SetupComputationalWork = _pinnedSolver->SetupComputationalWork;
	}

	// On several processes: the pinned DoF of each vector of the kernel is that of largest |z_i| on the smallest rank that
	// has it. Its owner pins its row and its column of Diag; the processes that have it as a ghost zero its column of
	// Offd. Collective.
	void Setup(const DistributedMatrix& A) override
	{
		this->Matrix = nullptr;
		this->Distributed = &A;
		const Communicator& comm = A.Comm();
		_z = Kernel(_interpOfOne, _blockSize, A.LocalRows(), _components);
		_pinnedDistributedA = A;
		_pins.clear();
		for (const Vector& z : _z)
		{
			Eigen::Index localPin = -1;
			double localMax = z.rows() > 0 ? z.cwiseAbs().maxCoeff(&localPin) : -1;
			int owner = comm.Min(localMax == comm.Max(localMax) ? comm.Rank() : comm.Size());
			Eigen::Index pin = comm.Rank() == owner ? localPin : -1;
			if (pin >= 0)
				_pins.push_back(pin);
			BigNumber pinGlobal = comm.Sum<BigNumber>(pin >= 0 ? A.RowMap->Begin * A.RowBlockSize + pin : 0);

			BigNumber ghostColumn = -1;
			const vector<BigNumber>& ghosts = A.ColMap->Ghosts;
			auto ghost = find(ghosts.begin(), ghosts.end(), pinGlobal / A.ColBlockSize);
			if (ghost != ghosts.end())
				ghostColumn = (ghost - ghosts.begin()) * A.ColBlockSize + pinGlobal % A.ColBlockSize;
			Pin(_pinnedDistributedA.Diag, pin, pin);
			Pin(_pinnedDistributedA.Offd, pin, ghostColumn);
		}
		_pinnedSolver->Setup(_pinnedDistributedA);
		this->SetupComputationalWork = _pinnedSolver->SetupComputationalWork;
	}

	Vector Solve(const Vector& b) override
	{
		Vector rhs = b;
		Project(rhs, _z);
		for (Eigen::Index pin : _pins)
			rhs[pin] = 0;
		Vector x = _pinnedSolver->Solve(rhs);
		Project(x, _z);
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
