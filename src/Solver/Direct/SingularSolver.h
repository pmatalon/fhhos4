#pragma once
#include "../Solver.h"

// The solver of a singular symmetric matrix whose kernel is spanned by the interpolation of 1 (the coarsest level of
// U-AMG for a pure Neumann problem, see UncondensedAMG::SingularMatrix): the kernel z is c on the first DoF of each
// block and 0 on the others (1 = c phi_0); with n components (UncondensedAMG::VectorComponents), n vectors z_k, c on the
// DoF k of each block (the constant of the component k), orthogonal (disjoint supports). Or given functions (e.g. the
// rigid body modes) on the first DoFs of each block, orthonormalized (see Kernel()). The matrix is pinned on one DoF per
// vector of the kernel (its row and column replaced by those of the identity), chosen so that the kernel restricted to
// the pinned DoFs is regular (see Setup()): the pinned matrix is regular, factorized by the solver given. Each solve
// projects the right-hand side onto the image of the matrix (orthogonal to the kernel), zeroes the pinned DoFs, solves,
// and returns the solution orthogonal to the kernel. Without the pins, the factorization meets pivots of the order of
// the rounding errors, and the solutions have a kernel component that ruins the Krylov method. The kernel is known: no
// threshold to detect the null pivots (as MUMPS's ICNTL(24) would need), and the same for all the factorizations, on
// one process or several.
class SingularSolver : public Solver
{
private:
	unique_ptr<Solver> _pinnedSolver;
	Vector _interpOfOne; // one coefficient per block; empty: 1
	int _blockSize = 1;
	int _components = 1;
	const DenseMatrix* _functions = nullptr; // see Kernel()
	int _functionDoFs = 0;
	vector<Vector> _z;            // the kernel, orthonormal vectors
	vector<Eigen::Index> _pins;   // the pinned DoFs (one per vector of the kernel), local: those of the process
	SparseMatrix _pinnedA;
	DistributedMatrix _pinnedDistributedA;
public:
	// pinnedSolver: owned
	SingularSolver(Solver* pinnedSolver)
		: _pinnedSolver(pinnedSolver)
	{}

	// The kernel on the n DoFs of the process, orthonormal vectors (on all the processes of Communicator::Current()): for
	// each of the components k, the interpolation of 1 c on the DoF k of each block of blockSize DoFs (1 if interpOfOne
	// is empty; its first n / blockSize coefficients: those of the owned blocks), 0 on the others. functions (if not
	// null): instead, its columns on the first functionDoFs DoFs of each block (functionDoFs rows per block), 0 on the
	// others, orthonormalized (modified Gram-Schmidt, twice)
	static vector<Vector> Kernel(const Vector& interpOfOne, int blockSize, Eigen::Index n, int components, const DenseMatrix* functions = nullptr, int functionDoFs = 0)
	{
		vector<Vector> kernel;
		if (functions)
		{
			for (Eigen::Index j = 0; j < functions->cols(); ++j)
			{
				Vector z = Vector::Zero(n);
				for (Eigen::Index b = 0; b < n / blockSize; ++b)
					z.segment(b * blockSize, functionDoFs) = functions->col(j).segment(b * functionDoFs, functionDoFs);
				kernel.push_back(z);
			}
			for (int pass = 0; pass < 2; ++pass)
			{
				for (size_t j = 0; j < kernel.size(); ++j)
				{
					for (size_t i = 0; i < j; ++i)
						kernel[j] -= Parallel::Dot(kernel[i], kernel[j]) * kernel[i];
					kernel[j] /= Parallel::Norm(kernel[j]);
				}
			}
			return kernel;
		}
		for (int k = 0; k < components; ++k)
		{
			Vector z = Vector::Zero(n);
			for (Eigen::Index i = k; i < n; i += blockSize)
				z[i] = interpOfOne.rows() > 0 ? interpOfOne[i / blockSize] : 1;
			kernel.push_back(z / Parallel::Norm(z));
		}
		return kernel;
	}

	// x projected orthogonally to the kernel (its vectors orthonormal: one after the other); returns the norm of the
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

	// The kernel of the matrix of the next setup, see Kernel(): the interpolation of 1 on the blocks of blockSize DoFs
	// (empty: 1 everywhere) and the number of components, whose constants are the first DoFs of each block; or the
	// functions (not copied: kept until the setup)
	void SetKernel(const Vector& interpOfOne, int blockSize, int components, const DenseMatrix* functions = nullptr, int functionDoFs = 0)
	{
		_interpOfOne = interpOfOne;
		_blockSize = blockSize;
		_components = components;
		_functions = functions;
		_functionDoFs = functionDoFs;
	}

	void Serialize(ostream& os) const override
	{
		size_t pins = _functions ? (size_t)_functions->cols() : (size_t)_components;
		os << *_pinnedSolver << " (singular matrix: " << (pins == 1 ? "one DoF" : to_string(pins) + " DoFs") << " pinned, projections onto the image)";
	}

	// The pinned DoFs: for each vector of the kernel in turn, minus its multiples of the previous ones that make it zero on
	// their pinned DoFs, the DoF of its largest |coefficient|. The kernel restricted to the pinned DoFs is then triangular
	// with a non-zero diagonal: regular. On vectors of disjoint supports (the constants of the components), the DoF of
	// largest |z_i| of each.
	void Setup(const SparseMatrix& A) override
	{
		Solver::Setup(A);
		_z = Kernel(_interpOfOne, _blockSize, A.rows(), _components, _functions, _functionDoFs);
		_pinnedA = A;
		_pins.clear();
		vector<Vector> eliminated = _z;
		for (size_t j = 0; j < _z.size(); ++j)
		{
			Vector& w = eliminated[j];
			for (size_t i = 0; i < j; ++i)
				w -= w[_pins[i]] / eliminated[i][_pins[i]] * eliminated[i];
			Eigen::Index pin;
			w.cwiseAbs().maxCoeff(&pin);
			_pins.push_back(pin);
			Pin(_pinnedA, pin, pin);
		}
		_pinnedSolver->Setup(_pinnedA);
		this->SetupComputationalWork = _pinnedSolver->SetupComputationalWork;
	}

	// On several processes: the pinned DoFs as on one process, the DoF of largest |coefficient| on the smallest rank that
	// has it, whose coefficients the eliminations take from its owner. Its owner pins its row and its column of Diag; the
	// processes that have it as a ghost zero its column of Offd. Collective.
	void Setup(const DistributedMatrix& A) override
	{
		this->Matrix = nullptr;
		this->Distributed = &A;
		const Communicator& comm = A.Comm();
		_z = Kernel(_interpOfOne, _blockSize, A.LocalRows(), _components, _functions, _functionDoFs);
		_pinnedDistributedA = A;
		_pins.clear();
		vector<Vector> eliminated = _z;
		vector<Eigen::Index> pins;  // of each vector: the local pinned DoF on its owner, -1 on the other processes
		vector<double> pivots;      // of each vector: its eliminated coefficient on its pinned DoF
		for (size_t j = 0; j < _z.size(); ++j)
		{
			Vector& w = eliminated[j];
			for (size_t i = 0; i < j; ++i)
				w -= comm.Sum(pins[i] >= 0 ? w[pins[i]] : 0.0) / pivots[i] * eliminated[i];
			Eigen::Index localPin = -1;
			double localMax = w.rows() > 0 ? w.cwiseAbs().maxCoeff(&localPin) : -1;
			int owner = comm.Min(localMax == comm.Max(localMax) ? comm.Rank() : comm.Size());
			Eigen::Index pin = comm.Rank() == owner ? localPin : -1;
			pins.push_back(pin);
			pivots.push_back(comm.Sum(pin >= 0 ? w[pin] : 0.0));
			if (pin >= 0)
				_pins.push_back(pin);
			GlobalNumber pinGlobal = comm.Sum<GlobalNumber>(pin >= 0 ? A.RowMap->Begin * A.RowBlockSize + pin : 0);

			BigNumber ghostColumn = -1;
			const vector<GlobalNumber>& ghosts = A.ColMap->Ghosts;
			auto ghost = find(ghosts.begin(), ghosts.end(), pinGlobal / A.ColBlockSize);
			if (ghost != ghosts.end())
				ghostColumn = (BigNumber)((ghost - ghosts.begin()) * A.ColBlockSize + pinGlobal % A.ColBlockSize);
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
