#pragma once
#include <numeric>
#include "../Solver.h"
using namespace std;

#ifdef MUMPS_ENABLED
#include <dmumps_c.h>

// MUMPS's sparse factorization (analysis and factorization in the setup, a solve per Solve()): LDL^T without pivoting
// of a symmetric positive definite matrix (SPD: MUMPS's SYM=1, from its lower triangle), LU otherwise (SYM=0). On
// several processes (the coarse solver of U-AMG), nothing is gathered: each process gives its rows of the matrix and of
// the right-hand side, MUMPS factorizes on all of them, and each process gets back its rows of the solution.
class MUMPSSolver : public Solver
{
public:
	const bool SPD;

private:
	DMUMPS_STRUC_C _id = {};
	bool _initialized = false;
	Communicator _self; // the processes of a matrix that is not distributed: this one
	const Communicator* _comm = nullptr;
	BigNumber _firstRow = 0;
	// The entries of the rows of the process (1-based global numbers; the lower triangle if SPD), read by the setup
	vector<MUMPS_INT> _irn, _jcn;
	vector<double> _a;
	// The rows of the process (1-based global numbers) and their right-hand side. MUMPS's arrays have at least one
	// element (a process may have no row).
	vector<MUMPS_INT> _rhsRows;
	vector<double> _rhs;
	// The solution as MUMPS distributes it: _nSolution rows, whose global numbers (1-based) MUMPS gives at the first
	// solve (fixed by the factorization). They go to their owners through _solutionMap: the rows of the process, then
	// as ghosts those of other processes, at _solutionPositions.
	MUMPS_INT _nSolution = 0;
	vector<double> _solution;
	vector<MUMPS_INT> _solutionRows;
	unique_ptr<IndexMap> _solutionMap;
	vector<BigNumber> _solutionPositions;

public:
	explicit MUMPSSolver(bool spd) : Solver(), SPD(spd) {}

	MUMPSSolver(const MUMPSSolver&) = delete; // MUMPS keeps pointers to the arrays
	MUMPSSolver& operator=(const MUMPSSolver&) = delete;

	void Serialize(ostream& os) const override
	{
		os << (SPD ? "LDL^T" : "LU") << " factorization (MUMPS)";
	}

	void Setup(const SparseMatrix& A) override
	{
		Solver::Setup(A);
		Factorize(A, 0, A.rows(), _self);
	}

	// Collective
	void Setup(const DistributedMatrix& A) override
	{
		Factorize(A.GlobalColumns(), A.RowMap->Begin * A.RowBlockSize, A.GlobalRows(), A.Comm());
		this->Matrix = &A.Diag;
		this->Distributed = &A;
	}

	// b: the rows of the process. Collective.
	Vector Solve(const Vector& b) override
	{
		BigNumber n = b.rows();
		copy(b.data(), b.data() + n, _rhs.begin());
		_id.job = 3;
		Run("the solve");
		if (!_solutionMap)
		{
			vector<BigNumber> ghosts;
			for (MUMPS_INT p = 0; p < _nSolution; p++)
			{
				BigNumber row = _solutionRows[p] - 1;
				if (row < _firstRow || row >= _firstRow + n)
					ghosts.push_back(row);
			}
			_solutionMap = make_unique<IndexMap>(IndexMap(*_comm, n), std::move(ghosts));
			_solutionPositions.resize(_nSolution);
			for (MUMPS_INT p = 0; p < _nSolution; p++)
				_solutionPositions[p] = _solutionMap->Local(_solutionRows[p] - 1);
		}
		Vector x = Vector::Zero(_solutionMap->NLocal());
		for (MUMPS_INT p = 0; p < _nSolution; p++)
			x[_solutionPositions[p]] = _solution[p];
		_solutionMap->AddGhostsToOwners(x.data() + n, x.data(), 1);
		this->SolvingComputationalWork = 0;
		return x.head(n);
	}

	~MUMPSSolver()
	{
		Terminate();
	}

private:
	// The matrix of size n, of which the process has the rows [firstRow, firstRow + rows.rows()) (global columns), on the
	// processes of comm. Collective.
	void Factorize(const SparseMatrix& rows, BigNumber firstRow, BigNumber n, const Communicator& comm)
	{
		Terminate();
		_comm = &comm;
		_firstRow = firstRow;
		_irn.clear();
		_jcn.clear();
		_a.clear();
		for (BigNumber i = 0; i < rows.rows(); i++)
		{
			for (SparseMatrix::InnerIterator it(rows, i); it; ++it)
			{
				if (!SPD || it.col() <= firstRow + i)
				{
					_irn.push_back((MUMPS_INT)(firstRow + i + 1));
					_jcn.push_back((MUMPS_INT)(it.col() + 1));
					_a.push_back(it.value());
				}
			}
		}
		BigNumber nRows = rows.rows();
		_rhsRows.resize(max<BigNumber>(1, nRows));
		iota(_rhsRows.begin(), _rhsRows.end(), (MUMPS_INT)(firstRow + 1));
		_rhs.resize(max<BigNumber>(1, nRows));

		_id = {};
		_id.comm_fortran = comm.FortranHandle();
		_id.par = 1; // the host process factorizes too
		_id.sym = SPD ? 1 : 0;
		_id.job = -1;
		Run("the initialization");
		_initialized = true;
		_id.icntl[0] = -1; _id.icntl[1] = -1; _id.icntl[2] = -1; _id.icntl[3] = 0; // silent
		_id.icntl[17] = 3;  // ICNTL(18): the entries given by each process
		_id.icntl[19] = 10; // ICNTL(20): the right-hand side given by each process
		_id.icntl[20] = 1;  // ICNTL(21): the solution left distributed
		_id.n = (MUMPS_INT)n;
		_id.nnz_loc = _a.size();
		_id.irn_loc = _irn.data();
		_id.jcn_loc = _jcn.data();
		_id.a_loc = _a.data();
		_id.job = 4; // analysis and factorization
		for (int attempt = 1; ; attempt++)
		{
			dmumps_c(&_id);
			// INFOG(1) = -9: MUMPS's estimate of its working space too small (e.g. MUMPS 5.7.3 on 3 processes): more space
			// (ICNTL(14), the percentage added to the estimate) and the factorization again, as MUMPS's guide says
			if (_id.infog[0] != -9 || attempt == 5)
				break;
			_id.icntl[13] = max(2 * _id.icntl[13], 40);
			_id.job = 2; // the factorization, on the same analysis
		}
		Check("the factorization");
		this->SetupComputationalWork = _id.rinfog[2] * 1e-6; // RINFOG(3): the flops of the elimination

		_nSolution = _id.info[22]; // INFO(23)
		_solution.resize(max<MUMPS_INT>(1, _nSolution));
		_solutionRows.resize(_solution.size());
		_id.sol_loc = _solution.data();
		_id.isol_loc = _solutionRows.data();
		_id.lsol_loc = _solution.size();
		_id.nrhs = 1;
		_id.nloc_rhs = (MUMPS_INT)nRows;
		_id.lrhs_loc = _rhs.size();
		_id.rhs_loc = _rhs.data();
		_id.irhs_loc = _rhsRows.data();
		_solutionMap.reset();
	}

	// One call of MUMPS (job set). Its errors (INFOG) are the same on all the processes.
	void Run(const string& what)
	{
		dmumps_c(&_id);
		Check(what);
	}

	void Check(const string& what) const
	{
		if (_id.infog[0] < 0)
			Utils::FatalError("MUMPS: " + what + " failed with INFOG(1) = " + to_string(_id.infog[0]) + ", INFOG(2) = " + to_string(_id.infog[1]) + " (see MUMPS's user guide).");
	}

	// Releases the instance of MUMPS (collective), unless MPI is finalized (e.g. at the exit of a Python program: the
	// memory goes with the process)
	void Terminate()
	{
		if (_initialized && !Communicator::Finalized())
		{
			_id.job = -2;
			dmumps_c(&_id);
		}
		_initialized = false;
	}
};

#else

// Compiled without MUMPS
class MUMPSSolver : public Solver
{
public:
	const bool SPD;

	explicit MUMPSSolver(bool spd) : Solver(), SPD(spd)
	{
		Utils::FatalError("fhhos4 is compiled without MUMPS: reconfigure it with -DENABLE_MUMPS=ON (or ENABLE_MPI=ON), or choose another coarse solver.");
	}

	void Serialize(ostream& os) const override {}
	Vector Solve(const Vector& b) override { return b; }
};

#endif
