#pragma once
#include <algorithm>
#include "../Utils/SparseMatrixOps.h"
#include "../Utils/Utils.h"
using namespace std;

// Static condensation of the cell unknowns of a hybrid system, computed by cells, in parallel. The system
//     [ A_TT  A_TF ] [x_T]   [b_T]
//     [ A_FT  A_FF ] [x_F] = [b_F],    A_TT block diagonal (one block of cellBlockSize rows per cell),
// becomes A x_F = b, with A = A_FF - A_FT A_TT^{-1} A_TF and b = b_F - A_FT A_TT^{-1} b_T; then
// x_T = A_TT^{-1} b_T - A_TT^{-1} A_TF x_F.
// The constructor factorizes the blocks of A_TT (Cholesky for a symmetric system, as the program's condensation; LU with
// partial pivoting otherwise) and computes W = A_TT^{-1} A_TF by their solves, whose rows of a cell are dense on the
// union of the columns of its rows in A_TF (the faces of the cell); the solves of A_TT^{-1} b_T use the factorizations
// too. Explicit inverses are not used: their rounding errors grow with the condition number of the blocks and would
// limit the attainable accuracy of the discrete solution. The rows of A are then computed in the pattern of A_FF (no
// separate sum), completed in the rows where A_FT W has other columns (couplings of faces that are zero in A_FF, e.g.
// dropped by its assembly), sorted only there. Bit for bit the results of the global product
// FullFromLower(A_FF) - Multiply(A_FT, W) (FullFromLower for a symmetric system), which computes each coefficient in
// the same order.
// Symmetric system (A_FT not given): A_TT and A_FF are read from their lower triangular parts, and A_FT = A_TF^T.
// Non-symmetric system: A_FT is read in place, by the constructor and by CondensedRHS() (not copied): it must stay alive
// and unchanged until the destruction of the object.
class StaticCondensation
{
private:
	int _blockSize;                // DoFs per cell
	bool _symmetric;
	// The factorizations of the blocks of A_TT, column-major, one after the other: Cholesky (the factor L in the lower
	// triangular part) for a symmetric system, LU with partial pivoting otherwise (with the row permutation of each block)
	vector<double> _factors;
	vector<int> _permutations;
	SparseMatrix _W;               // A_TT^{-1} A_TF
	SparseMatrix _transposedA_TF;  // A_FT of a symmetric system; empty otherwise
	const SparseMatrix& _A_FT;     // the caller's, or _transposedA_TF
public:
	// A_FT: null if the system is symmetric
	StaticCondensation(const SparseMatrix& A_TT, const SparseMatrix& A_TF, const SparseMatrix* A_FT, int cellBlockSize)
		: _blockSize(cellBlockSize), _symmetric(!A_FT), _transposedA_TF(A_FT ? SparseMatrix() : SparseMatrixOps::Transpose(A_TF)),
		  _A_FT(A_FT ? *A_FT : _transposedA_TF)
	{
		assert(A_TT.rows() == A_TF.rows() && A_TT.rows() % cellBlockSize == 0);
		bool symmetric = !A_FT;
		BigNumber nCells = A_TT.rows() / cellBlockSize;
		int nT = cellBlockSize;

		// The columns of a cell: the union of the columns of its rows in A_TF, sorted
		auto cellColumns = [&](BigNumber T, vector<SparseMatrixIndex>& cols)
		{
			cols.clear();
			for (int r = 0; r < nT; ++r)
			{
				const SparseMatrixIndex* begin = A_TF.innerIndexPtr() + A_TF.outerIndexPtr()[T * nT + r];
				const SparseMatrixIndex* end = A_TF.innerIndexPtr() + A_TF.outerIndexPtr()[T * nT + r + 1];
				if (r == 0)
					cols.assign(begin, end);
				else if (!equal(begin, end, cols.begin(), cols.end())) // all the same in HHO
				{
					vector<SparseMatrixIndex> merged;
					set_union(cols.begin(), cols.end(), begin, end, back_inserter(merged));
					cols = std::move(merged);
				}
			}
		};

		// Pattern of W: the rows of a cell have its columns
		_W = SparseMatrix(A_TF.rows(), A_TF.cols());
		SparseMatrixIndex* outer = _W.outerIndexPtr();
		ThreadLocal<vector<SparseMatrixIndex>> columns;
		#pragma omp parallel for
		for (BigNumber T = 0; T < nCells; ++T)
		{
			vector<SparseMatrixIndex>& cols = columns.Local();
			cellColumns(T, cols);
			for (int r = 0; r < nT; ++r)
				outer[T * nT + r + 1] = cols.size();
		}
		outer[0] = 0;
		for (BigNumber i = 0; i < (BigNumber)A_TF.rows(); ++i)
			outer[i + 1] += outer[i];
		_W.resizeNonZeros(outer[A_TF.rows()]);

		// Factorizations of the blocks of A_TT and values of W
		_factors.resize(nCells * nT * nT);
		if (!symmetric)
			_permutations.resize(nCells * nT);
		ThreadLocal<DenseMatrix> rowsOfA_TF;
		BigNumber notPositiveDefinite = nCells; // the first cell whose block is not (symmetric system)
		#pragma omp parallel for reduction(min: notPositiveDefinite)
		for (BigNumber T = 0; T < nCells; ++T)
		{
			BigNumber first = T * nT;
			DenseMatrix block = DenseMatrix::Zero(nT, nT);
			for (int r = 0; r < nT; ++r)
			{
				for (SparseMatrix::InnerIterator it(A_TT, first + r); it; ++it)
				{
					BigNumber c = it.col() - first;
					if (it.col() < (SparseMatrixIndex)first || c >= (BigNumber)nT || (symmetric && (int)c > r))
						continue;
					block(r, c) = it.value();
					if (symmetric)
						block(c, r) = it.value();
				}
			}

			vector<SparseMatrixIndex>& cols = columns.Local();
			cellColumns(T, cols);
			int m = cols.size();
			DenseMatrix& R = rowsOfA_TF.Local(); // the rows of the cell in A_TF, on its columns
			R.setZero(nT, m);
			for (int r = 0; r < nT; ++r)
				for (SparseMatrix::InnerIterator it(A_TF, first + r); it; ++it)
					R(r, lower_bound(cols.begin(), cols.end(), it.col()) - cols.begin()) = it.value();

			Eigen::Map<DenseMatrix> factor(_factors.data() + T * nT * nT, nT, nT);
			DenseMatrix W_T;
			if (!CellSolve(block, R, factor, symmetric ? nullptr : _permutations.data() + T * nT, W_T))
				notPositiveDefinite = min(notPositiveDefinite, T);
			for (int r = 0; r < nT; ++r)
			{
				SparseMatrixIndex k = outer[first + r];
				for (int j = 0; j < m; ++j)
				{
					_W.innerIndexPtr()[k + j] = cols[j];
					_W.valuePtr()[k + j] = W_T(r, j);
				}
			}
		}
		if (notPositiveDefinite < nCells)
			Utils::FatalError("fhhos4::Solver: the block of A_TT of the cell " + to_string(notPositiveDefinite) + " is not positive definite (symmetric system: the Cholesky factorization of the blocks of A_TT).");
	}

	// W_T = A_TT_T^{-1} R by the factorization of the block A_TT_T: Cholesky (permutation null), its factor written in
	// factor, or LU with partial pivoting, its factors in factor and its row permutation in permutation. Returns false if
	// the Cholesky factorization fails (the block is not positive definite). Also used by the tests, as the reference.
	static bool CellSolve(const DenseMatrix& block, const DenseMatrix& R, Eigen::Map<DenseMatrix>& factor, int* permutation, DenseMatrix& W_T)
	{
		if (!permutation)
		{
			Eigen::LLT<DenseMatrix> llt(block);
			factor = llt.matrixLLT();
			W_T = llt.solve(R);
			return llt.info() == Eigen::Success;
		}
		Eigen::PartialPivLU<DenseMatrix> lu(block);
		factor = lu.matrixLU();
		for (int k = 0; k < block.rows(); ++k)
			permutation[k] = lu.permutationP().indices()(k);
		W_T = lu.solve(R);
		return true;
	}

	// Not copyable: _A_FT may refer to a member
	StaticCondensation(const StaticCondensation&) = delete;
	StaticCondensation& operator=(const StaticCondensation&) = delete;

	// A = A_FF - A_FT A_TT^{-1} A_TF, both triangles. lower: A_FF read from its lower triangular part (symmetric system).
	SparseMatrix CondensedMatrix(const SparseMatrix& A_FF, bool lower) const
	{
		return AddProduct(A_FF, lower, -1);
	}

	// A_FF = A + A_FT A_TT^{-1} A_TF, A read in full
	SparseMatrix FaceBlock(const SparseMatrix& A) const
	{
		return AddProduct(A, false, 1);
	}

	// A_TT^{-1} b_T, by the factorizations of the blocks (the steps of Eigen's LLT::solve() and PartialPivLU::solve())
	Vector SolveCells(const Vector& b_T) const
	{
		int nT = _blockSize;
		Vector y(b_T.rows());
		#pragma omp parallel for
		for (BigNumber T = 0; T < (BigNumber)b_T.rows() / nT; ++T)
		{
			Eigen::Map<const DenseMatrix> factor(_factors.data() + T * nT * nT, nT, nT);
			auto y_T = y.segment(T * nT, nT);
			auto b = b_T.segment(T * nT, nT);
			if (_symmetric)
			{
				y_T = b;
				factor.triangularView<Eigen::Lower>().solveInPlace(y_T);
				factor.triangularView<Eigen::Lower>().adjoint().solveInPlace(y_T);
			}
			else
			{
				const int* permutation = _permutations.data() + T * nT;
				for (int k = 0; k < nT; ++k)
					y_T(permutation[k]) = b(k);
				factor.triangularView<Eigen::UnitLower>().solveInPlace(y_T);
				factor.triangularView<Eigen::Upper>().solveInPlace(y_T);
			}
		}
		return y;
	}

	// b = b_F - A_FT A_TT^{-1} b_T, from invA_TT_b_T = A_TT^{-1} b_T
	Vector CondensedRHS(const Vector& b_F, const Vector& invA_TT_b_T) const
	{
		return b_F - _A_FT * invA_TT_b_T;
	}

	// x_T = A_TT^{-1} b_T - A_TT^{-1} A_TF x_F, from invA_TT_b_T = A_TT^{-1} b_T
	Vector CellUnknowns(const Vector& invA_TT_b_T, const Vector& x_F) const
	{
		return invA_TT_b_T - _W * x_F;
	}

private:
	// M + sign A_FT W, M read in full or from its lower triangular part: in the pattern of M (in full), completed in the
	// rows where A_FT W has other columns (couplings of the faces of a cell that are zero in M, e.g. dropped by its
	// assembly), whose coefficients are then 0 + sign (A_FT W)(i, j), as in the sum of the global products
	SparseMatrix AddProduct(const SparseMatrix& M, bool lower, int sign) const
	{
		// M may have fewer rows than A_FT (read in full): the others are empty (the faces of other processes, see
		// DistributedStaticCondensation)
		assert((M.rows() == _A_FT.rows() || (M.rows() < _A_FT.rows() && !lower)) && M.cols() == _W.cols());
		BigNumber n = _A_FT.rows();
		BigNumber nRowsOfM = M.rows();
		SparseMatrix Mt; // lower: row i of Mt, column i of M, whose part under the diagonal gives the part of row i above it
		if (lower)
			Mt = SparseMatrixOps::Transpose(M);
		// The row i of M in full: [begin1, end1) in M, [begin2, end2) in Mt (empty if !lower)
		auto rowOfM = [&](BigNumber i, SparseMatrixIndex& begin1, SparseMatrixIndex& end1, SparseMatrixIndex& begin2, SparseMatrixIndex& end2)
		{
			begin1 = i < nRowsOfM ? M.outerIndexPtr()[i] : 0;
			end1 = i < nRowsOfM ? M.outerIndexPtr()[i + 1] : 0;
			begin2 = end2 = 0;
			if (lower)
			{
				end1 = upper_bound(M.innerIndexPtr() + begin1, M.innerIndexPtr() + end1, (SparseMatrixIndex)i) - M.innerIndexPtr();
				begin2 = upper_bound(Mt.innerIndexPtr() + Mt.outerIndexPtr()[i], Mt.innerIndexPtr() + Mt.outerIndexPtr()[i + 1], (SparseMatrixIndex)i) - Mt.innerIndexPtr();
				end2 = Mt.outerIndexPtr()[i + 1];
			}
		};

		ThreadLocal<vector<SparseMatrixIndex>> positions((size_t)M.cols(), (SparseMatrixIndex)-1); // in the row of S, -1 if absent
		ThreadLocal<vector<SparseMatrixIndex>> extraColumns;
		// The columns of the row i of A_FT W that are not in the row of M, in extra (unsorted). pos: all -1 before and after.
		auto extraColumnsOf = [&](BigNumber i, SparseMatrixIndex begin1, SparseMatrixIndex end1, SparseMatrixIndex begin2, SparseMatrixIndex end2, vector<SparseMatrixIndex>& pos, vector<SparseMatrixIndex>& extra)
		{
			for (SparseMatrixIndex k = begin1; k < end1; ++k)
				pos[M.innerIndexPtr()[k]] = 0;
			for (SparseMatrixIndex k = begin2; k < end2; ++k)
				pos[Mt.innerIndexPtr()[k]] = 0;
			extra.clear();
			BigNumber currentCell = (BigNumber)-1;
			for (SparseMatrix::InnerIterator a(_A_FT, i); a; ++a)
			{
				if ((BigNumber)a.col() / _blockSize == currentCell) // the rows of W of a cell have the same columns
					continue;
				currentCell = a.col() / _blockSize;
				for (SparseMatrixIndex k = _W.outerIndexPtr()[a.col()]; k < _W.outerIndexPtr()[a.col() + 1]; ++k)
				{
					SparseMatrixIndex col = _W.innerIndexPtr()[k];
					if (pos[col] < 0)
					{
						pos[col] = 0;
						extra.push_back(col);
					}
				}
			}
			for (SparseMatrixIndex k = begin1; k < end1; ++k)
				pos[M.innerIndexPtr()[k]] = -1;
			for (SparseMatrixIndex k = begin2; k < end2; ++k)
				pos[Mt.innerIndexPtr()[k]] = -1;
			for (SparseMatrixIndex col : extra)
				pos[col] = -1;
		};

		// Sizes of the rows
		SparseMatrix S(n, M.cols());
		SparseMatrixIndex* outer = S.outerIndexPtr();
		#pragma omp parallel for schedule(dynamic, 1024)
		for (BigNumber i = 0; i < n; ++i)
		{
			SparseMatrixIndex begin1, end1, begin2, end2;
			rowOfM(i, begin1, end1, begin2, end2);
			vector<SparseMatrixIndex>& extra = extraColumns.Local();
			extraColumnsOf(i, begin1, end1, begin2, end2, positions.Local(), extra);
			outer[i + 1] = (end1 - begin1) + (end2 - begin2) + extra.size();
		}
		outer[0] = 0;
		for (BigNumber i = 0; i < n; ++i)
			outer[i + 1] += outer[i];
		S.resizeNonZeros(outer[n]);

		ThreadLocal<vector<double>> products;
		ThreadLocal<vector<SparseMatrixIndex>> cellPositions; // of the columns of the current cell in the row of S
		ThreadLocal<vector<pair<SparseMatrixIndex, double>>> completedRows;
		#pragma omp parallel for schedule(dynamic, 1024)
		for (BigNumber i = 0; i < n; ++i)
		{
			vector<SparseMatrixIndex>& pos = positions.Local();
			vector<double>& p = products.Local();
			vector<SparseMatrixIndex>& cellPos = cellPositions.Local();

			// The row of M, completed by the other columns of A_FT W (zeros of M), sorted
			SparseMatrixIndex begin1, end1, begin2, end2;
			rowOfM(i, begin1, end1, begin2, end2);
			SparseMatrixIndex first = outer[i];
			SparseMatrixIndex* cols = S.innerIndexPtr() + first;
			double* values = S.valuePtr() + first;
			copy(M.innerIndexPtr() + begin1, M.innerIndexPtr() + end1, cols);
			copy(M.valuePtr() + begin1, M.valuePtr() + end1, values);
			if (lower)
			{
				copy(Mt.innerIndexPtr() + begin2, Mt.innerIndexPtr() + end2, cols + (end1 - begin1));
				copy(Mt.valuePtr() + begin2, Mt.valuePtr() + end2, values + (end1 - begin1));
			}
			SparseMatrixIndex sizeOfM = (end1 - begin1) + (end2 - begin2);
			SparseMatrixIndex size = outer[i + 1] - first;
			if (size > sizeOfM)
			{
				vector<SparseMatrixIndex>& extra = extraColumns.Local();
				extraColumnsOf(i, begin1, end1, begin2, end2, pos, extra);
				vector<pair<SparseMatrixIndex, double>>& row = completedRows.Local();
				row.clear();
				for (SparseMatrixIndex k = 0; k < sizeOfM; ++k)
					row.push_back({ cols[k], values[k] });
				for (SparseMatrixIndex col : extra)
					row.push_back({ col, 0.0 });
				sort(row.begin(), row.end(), [](const pair<SparseMatrixIndex, double>& x, const pair<SparseMatrixIndex, double>& y) { return x.first < y.first; });
				for (SparseMatrixIndex k = 0; k < size; ++k)
				{
					cols[k] = row[k].first;
					values[k] = row[k].second;
				}
			}
			for (SparseMatrixIndex k = 0; k < size; ++k)
				pos[cols[k]] = k;

			// (A_FT W)(i, :), in the order of Multiply(A_FT, W): the columns c of the row of A_FT, then the columns of W
			p.assign(size, 0.0);
			BigNumber currentCell = (BigNumber)-1;
			for (SparseMatrix::InnerIterator a(_A_FT, i); a; ++a)
			{
				SparseMatrixIndex start = _W.outerIndexPtr()[a.col()];
				SparseMatrixIndex m = _W.outerIndexPtr()[a.col() + 1] - start;
				if ((BigNumber)a.col() / _blockSize != currentCell) // the rows of W of a cell have the same columns
				{
					currentCell = a.col() / _blockSize;
					cellPos.resize(m);
					for (SparseMatrixIndex j = 0; j < m; ++j)
						cellPos[j] = pos[_W.innerIndexPtr()[start + j]];
				}
				const double* w = _W.valuePtr() + start;
				for (SparseMatrixIndex j = 0; j < m; ++j)
					p[cellPos[j]] += w[j] * a.value();
			}
			for (SparseMatrixIndex k = 0; k < size; ++k)
			{
				values[k] = sign < 0 ? values[k] - p[k] : values[k] + p[k];
				pos[cols[k]] = -1;
			}
		}
		return S;
	}
};
