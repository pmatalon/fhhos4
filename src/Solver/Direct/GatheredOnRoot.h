#pragma once
#include <limits>
#include "../Solver.h"
using namespace std;

// The coarse solve of the reproducible modes on several processes (ADR-0024), in place of MUMPS, whose solutions vary
// from run to run there: the matrix gathered on the process of rank 0, its rows in the global numbering (the rows of
// GlobalColumns(), MUMPS's input) in rank order, so that it doesn't depend on the order in which the messages arrive;
// factorized there by the solver's own code, the one of one process (Eigen's); each right-hand side gathered there,
// solved, and the solution scattered to the owners of its rows. Every step is collective. The process 0 holds and
// factorizes the whole matrix: a serial step, which doesn't scale to many processes (ADR-0021: MUMPS by default).
class GatheredOnRoot
{
private:
	const Communicator& _comm;
	vector<int> _rows; // of each process, on the process 0

public:
	// factorize(M): M, the whole matrix, on the process 0 only; its errors raised on all the processes
	template <class Factorize>
	GatheredOnRoot(const DistributedMatrix& A, Factorize factorize)
		: _comm(A.Comm())
	{
		GlobalColumnsMatrix rows = A.GlobalColumns();
		vector<GlobalNumber> rowNonZeros(rows.rows());
		for (BigNumber i = 0; i < rows.rows(); i++)
			rowNonZeros[i] = rows.outerIndexPtr()[i + 1] - rows.outerIndexPtr()[i];
		const GlobalNumber* columns = rows.innerIndexPtr();
		vector<GlobalNumber> allRowNonZeros = _comm.GatherV(rowNonZeros, 0, &_rows);
		vector<GlobalNumber> allColumns = _comm.GatherV(vector<GlobalNumber>(columns, columns + rows.nonZeros()));
		vector<double> allValues = _comm.GatherV(vector<double>(rows.valuePtr(), rows.valuePtr() + rows.nonZeros()));
		GlobalNumber n = A.GlobalRows();
		CollectiveCheck(_comm, [&]()
		{
			if (!_comm.IsRoot())
				return;
			if (n > numeric_limits<SparseMatrixIndex>::max() || (GlobalNumber)allValues.size() > numeric_limits<SparseMatrixIndex>::max())
				Utils::FatalError("The coarse matrix gathered on the process 0 (" + to_string(n) + " rows, " + to_string(allValues.size()) + " non-zeros) exceeds the index type: compile fhhos4 with -DENABLE_64BIT_INDEX=ON.");
			vector<SparseMatrixIndex> rowPointers(n + 1, 0);
			for (GlobalNumber i = 0; i < n; i++)
				rowPointers[i + 1] = rowPointers[i] + (SparseMatrixIndex)allRowNonZeros[i];
			vector<SparseMatrixIndex> columnIndices(allColumns.begin(), allColumns.end());
			SparseMatrix M = Eigen::Map<const SparseMatrix>(n, n, allValues.size(), rowPointers.data(), columnIndices.data(), allValues.data());
			factorize(M);
		});
	}

	// b: the rows of the process; solve(all): on the process 0, the solution of the whole matrix
	template <class Solve>
	Vector SolveWith(const Vector& b, Solve solve) const
	{
		vector<double> all = _comm.GatherV(vector<double>(b.data(), b.data() + b.rows()));
		if (_comm.IsRoot())
		{
			Vector x = solve(Eigen::Map<const Vector>(all.data(), all.size()));
			all.assign(x.data(), x.data() + x.rows());
		}
		vector<double> local = _comm.ScatterV(all, _rows, b.rows());
		return Eigen::Map<const Vector>(local.data(), local.size());
	}
};
