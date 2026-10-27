#pragma once
#include <memory>
#include "IndexMap.h"
#include "SparseMatrixOps.h"
using namespace std;

// A sparse matrix distributed by rows over the processes (on several of them; on one, the solvers use the SparseMatrix
// itself). Each process holds the rows of the entities it owns in RowMap (blocks of RowBlockSize rows per entity),
// split into the columns of the entities it owns in ColMap (Diag) and those of its ghosts in ColMap (Offd, the columns
// in the order of the ghosts, by blocks of ColBlockSize). For a square operator (RowMap and ColMap with the same
// ownership), Diag is the square block of the process: the smoothers work on it.
class DistributedMatrix
{
public:
	shared_ptr<const IndexMap> RowMap;
	shared_ptr<const IndexMap> ColMap;
	int RowBlockSize = 1;
	int ColBlockSize = 1;
	SparseMatrix Diag;
	SparseMatrix Offd;

private:
	// The transposes of Diag and Offd, for MultiplyTransposed() (SetupTransposed())
	SparseMatrix _diagT;
	SparseMatrix _offdT;

public:
	DistributedMatrix() {}

	const Communicator& Comm() const { return *RowMap->Comm; }
	BigNumber LocalRows() const { return Diag.rows(); }
	BigNumber GlobalRows() const { return RowMap->GlobalSize * RowBlockSize; }
	BigNumber GlobalCols() const { return ColMap->GlobalSize * ColBlockSize; }
	BigNumber LocalNonZeros() const { return Diag.nonZeros() + Offd.nonZeros(); }
	BigNumber GlobalNonZeros() const { return Comm().Sum(LocalNonZeros()); }

	//-------------------------------------------------------------------//
	//                          Construction                             //
	//-------------------------------------------------------------------//

	// From the rows of the owned entities of rowMap, whose columns are global (in [0, colOwners.GlobalSize * colBlockSize),
	// sorted in each row): the ghosts of the column map are the entities of the columns not owned. Collective.
	static DistributedMatrix FromGlobalColumns(shared_ptr<const IndexMap> rowMap, const IndexMap& colOwners, const SparseMatrix& rows, int rowBlockSize, int colBlockSize)
	{
		assert(rows.rows() == rowMap->NOwned * rowBlockSize);
		// The entities of the columns not owned, each one once per thread (sorted and unique), then merged
		BigNumber first = colOwners.Begin * colBlockSize, end = (colOwners.Begin + colOwners.NOwned) * colBlockSize;
		ThreadLocal<vector<BigNumber>> ghostsTL;
		#pragma omp parallel
		{
			vector<BigNumber>& ghosts = ghostsTL.Local();
			#pragma omp for schedule(static)
			for (BigNumber i = 0; i < (BigNumber)rows.rows(); i++)
			{
				BigNumber previous = -1;
				for (SparseMatrix::InnerIterator it(rows, i); it; ++it)
				{
					BigNumber col = it.col();
					if ((col < first || col >= end) && col / colBlockSize != previous)
					{
						previous = col / colBlockSize;
						ghosts.push_back(previous);
					}
				}
			}
			sort(ghosts.begin(), ghosts.end());
			ghosts.erase(unique(ghosts.begin(), ghosts.end()), ghosts.end());
		}
		vector<BigNumber> ghosts;
		for (vector<BigNumber>& g : ghostsTL)
			ghosts.insert(ghosts.end(), g.begin(), g.end());
		auto colMap = make_shared<const IndexMap>(colOwners, std::move(ghosts));
		return FromGlobalColumns(rowMap, colMap, rows, rowBlockSize, colBlockSize);
	}

	// Same, with the column map given (its ghosts must contain the entities of the columns not owned)
	static DistributedMatrix FromGlobalColumns(shared_ptr<const IndexMap> rowMap, shared_ptr<const IndexMap> colMap, const SparseMatrix& rows, int rowBlockSize, int colBlockSize)
	{
		DistributedMatrix M;
		M.RowMap = rowMap;
		M.ColMap = colMap;
		M.RowBlockSize = rowBlockSize;
		M.ColBlockSize = colBlockSize;
		int cbs = colBlockSize;
		const IndexMap& cm = *colMap;
		BigNumber first = cm.Begin * cbs, end = (cm.Begin + cm.NOwned) * cbs;
		auto local = [&](BigNumber col) // local column: owned, or NOwned * cbs + ghost column
		{
			if (col >= first && col < end)
				return col - first;
			BigNumber l = cm.Local(col / cbs);
			if (l < 0)
				Utils::FatalError("DistributedMatrix: a column refers to an entity that is neither owned nor a ghost of the column map.");
			return l * cbs + col % cbs;
		};
		BigNumber nOwnedCols = cm.NOwned * cbs;
		M.Diag = SparseMatrixOps::BuildByRows(rows.rows(), nOwnedCols, 1, [&](BigNumber firstRow, BigNumber lastRow, SparseMatrixOps::RowChunk& chunk)
		{
			for (BigNumber i = firstRow; i < lastRow; ++i)
			{
				for (SparseMatrix::InnerIterator it(rows, i); it; ++it)
					if (it.col() >= first && it.col() < end)
						chunk.Add(it.col() - first, it.value());
				chunk.EndRow();
			}
		});
		ThreadLocal<vector<pair<SparseMatrixIndex, double>>> rowsTL;
		M.Offd = SparseMatrixOps::BuildByRows(rows.rows(), cm.NGhosts() * cbs, 1, [&](BigNumber firstRow, BigNumber lastRow, SparseMatrixOps::RowChunk& chunk)
		{
			vector<pair<SparseMatrixIndex, double>>& row = rowsTL.Local();
			for (BigNumber i = firstRow; i < lastRow; ++i)
			{
				row.clear();
				for (SparseMatrix::InnerIterator it(rows, i); it; ++it)
					if (it.col() < first || it.col() >= end)
						row.push_back({ (SparseMatrixIndex)(local(it.col()) - nOwnedCols), it.value() });
				// Sorted if the ghosts follow the global numbers (the maps of fhhos4)
				auto byColumn = [](const pair<SparseMatrixIndex, double>& x, const pair<SparseMatrixIndex, double>& y) { return x.first < y.first; };
				if (!is_sorted(row.begin(), row.end(), byColumn))
					sort(row.begin(), row.end(), byColumn);
				for (const auto& coeff : row)
					chunk.Add(coeff.first, coeff.second);
				chunk.EndRow();
			}
		});
		return M;
	}

	// From the local matrix of the rows of the owned entities, whose columns are local in colMap (owned, then ghosts)
	static DistributedMatrix FromLocalColumns(shared_ptr<const IndexMap> rowMap, shared_ptr<const IndexMap> colMap, const SparseMatrix& rows, int rowBlockSize, int colBlockSize)
	{
		assert(rows.rows() == rowMap->NOwned * rowBlockSize && rows.cols() == colMap->NLocal() * colBlockSize);
		DistributedMatrix M;
		M.RowMap = rowMap;
		M.ColMap = colMap;
		M.RowBlockSize = rowBlockSize;
		M.ColBlockSize = colBlockSize;
		BigNumber nOwnedCols = colMap->NOwned * colBlockSize;
		M.Diag = SplitColumns(rows, 0, nOwnedCols);
		M.Offd = SplitColumns(rows, nOwnedCols, rows.cols());
		return M;
	}

	// The columns [begin, end) of M, renumbered from 0
	static SparseMatrix SplitColumns(const SparseMatrix& M, BigNumber begin, BigNumber end)
	{
		return SparseMatrixOps::BuildByRows(M.rows(), end - begin, 1, [&](BigNumber firstRow, BigNumber lastRow, SparseMatrixOps::RowChunk& chunk)
		{
			for (BigNumber i = firstRow; i < lastRow; ++i)
			{
				for (SparseMatrix::InnerIterator it(M, i); it; ++it)
					if (it.col() >= begin && it.col() < end)
						chunk.Add(it.col() - begin, it.value());
				chunk.EndRow();
			}
		});
	}

	// The local matrix [Diag Offd]: the columns local in ColMap (owned, then ghosts)
	SparseMatrix LocalColumns() const
	{
		BigNumber nDiag = Diag.cols();
		return SparseMatrixOps::BuildByRows(Diag.rows(), nDiag + Offd.cols(), 1, [&](BigNumber firstRow, BigNumber lastRow, SparseMatrixOps::RowChunk& chunk)
		{
			for (BigNumber i = firstRow; i < lastRow; ++i)
			{
				for (SparseMatrix::InnerIterator it(Diag, i); it; ++it)
					chunk.Add(it.col(), it.value());
				for (SparseMatrix::InnerIterator it(Offd, i); it; ++it)
					chunk.Add(nDiag + it.col(), it.value());
				chunk.EndRow();
			}
		});
	}

	// The rows of the process with global columns
	SparseMatrix GlobalColumns() const
	{
		int cbs = ColBlockSize;
		const IndexMap& cm = *ColMap;
		BigNumber first = cm.Begin * cbs;
		ThreadLocal<vector<pair<SparseMatrixIndex, double>>> rowsTL;
		return SparseMatrixOps::BuildByRows(Diag.rows(), GlobalCols(), 1, [&](BigNumber firstRow, BigNumber lastRow, SparseMatrixOps::RowChunk& chunk)
		{
			vector<pair<SparseMatrixIndex, double>>& row = rowsTL.Local();
			for (BigNumber i = firstRow; i < lastRow; ++i)
			{
				row.clear();
				for (SparseMatrix::InnerIterator it(Diag, i); it; ++it)
					row.push_back({ (SparseMatrixIndex)(first + it.col()), it.value() });
				for (SparseMatrix::InnerIterator it(Offd, i); it; ++it)
					row.push_back({ (SparseMatrixIndex)(cm.Ghosts[it.col() / cbs] * cbs + it.col() % cbs), it.value() });
				sort(row.begin(), row.end(), [](const auto& x, const auto& y) { return x.first < y.first; });
				for (const auto& coeff : row)
					chunk.Add(coeff.first, coeff.second);
				chunk.EndRow();
			}
		});
	}

	//-------------------------------------------------------------------//
	//              Rows exchanged between the processes                 //
	//-------------------------------------------------------------------//

	// The rows of the ghosts of map, from their owners: ownedRows has the rows of the owned entities (NOwned * blockSize
	// rows, columns in any numbering common to the processes, e.g. global); returns those of the ghosts (NGhosts *
	// blockSize rows, same columns). Collective on the processes of the halo.
	static SparseMatrix FetchGhostRows(const IndexMap& map, const SparseMatrix& ownedRows, int blockSize)
	{
		assert(ownedRows.rows() == map.NOwned * blockSize);
		vector<vector<SparseMatrixIndex>> sendIndices(map.SendRanks.size());
		vector<vector<double>> sendValues(map.SendRanks.size());
		for (size_t k = 0; k < map.SendRanks.size(); k++)
		{
			for (BigNumber s = map.SendStart[k]; s < map.SendStart[k + 1]; s++)
				for (int r = 0; r < blockSize; r++)
					PackRow(ownedRows, map.SendIndices[s] * blockSize + r, sendIndices[k], sendValues[k]);
		}
		vector<vector<SparseMatrixIndex>> recvIndices = map.ForwardExchange(sendIndices);
		vector<vector<double>> recvValues = map.ForwardExchange(sendValues);
		vector<Eigen::Triplet<double, SparseMatrixIndex>> triplets;
		for (size_t k = 0; k < map.RecvRanks.size(); k++)
		{
			size_t pi = 0, pv = 0;
			for (BigNumber q = map.RecvStart[k]; q < map.RecvStart[k + 1]; q++)
				for (int r = 0; r < blockSize; r++)
					UnpackRow(map.RecvIndices[q] * blockSize + r, recvIndices[k], pi, recvValues[k], pv, triplets);
		}
		SparseMatrix ghostRows(map.NGhosts() * blockSize, ownedRows.cols());
		ghostRows.setFromTriplets(triplets.begin(), triplets.end());
		return ghostRows;
	}

	// The rows of all the local entities of rowMap (the ownership of RowMap: its owned ones, then its ghosts, whose rows
	// are fetched from their owners), with their columns of the local entities of colMap (the others dropped), numbered
	// locally in colMap. Collective.
	SparseMatrix RowsWithGhosts(const IndexMap& rowMap, const IndexMap& colMap) const
	{
		SparseMatrix owned = GlobalColumns();
		SparseMatrix ghosts = FetchGhostRows(rowMap, owned, RowBlockSize);
		int cbs = ColBlockSize;
		BigNumber nOwnedRows = owned.rows();
		ThreadLocal<vector<pair<SparseMatrixIndex, double>>> rowsTL;
		return SparseMatrixOps::BuildByRows(rowMap.NLocal() * RowBlockSize, colMap.NLocal() * cbs, 1, [&](BigNumber firstRow, BigNumber lastRow, SparseMatrixOps::RowChunk& chunk)
		{
			vector<pair<SparseMatrixIndex, double>>& row = rowsTL.Local();
			for (BigNumber i = firstRow; i < lastRow; ++i)
			{
				row.clear();
				const SparseMatrix& M = i < nOwnedRows ? owned : ghosts;
				for (SparseMatrix::InnerIterator it(M, i < nOwnedRows ? i : i - nOwnedRows); it; ++it)
				{
					BigNumber c = colMap.Local(it.col() / cbs);
					if (c >= 0)
						row.push_back({ (SparseMatrixIndex)(c * cbs + it.col() % cbs), it.value() });
				}
				sort(row.begin(), row.end(), [](const auto& x, const auto& y) { return x.first < y.first; });
				for (const auto& coeff : row)
					chunk.Add(coeff.first, coeff.second);
				chunk.EndRow();
			}
		});
	}

	// The operator assembled by the processes from the rows of their local entities in map (owned, then ghosts), with
	// columns local in map: a square operator on the entities of map, distributed by rows. The rows of the ghosts are
	// sent to their owners, which add them to theirs; the blocks Diag and Offd built in one pass on the rows. Collective.
	static DistributedMatrix AssembleFromLocalRows(shared_ptr<const IndexMap> map, const SparseMatrix& localRows, int blockSize)
	{
		const IndexMap& m = *map;
		int bs = blockSize;
		BigNumber nOwnedRows = m.NOwned * bs, nOwnedCols = m.NOwned * bs;
		assert(localRows.rows() == m.NLocal() * bs && localRows.cols() == m.NLocal() * bs);
		auto globalColumn = [&m, bs](SparseMatrixIndex c) { return (SparseMatrixIndex)(m.Global(c / bs) * bs + c % bs); };

		// The rows of the ghosts to their owners, with global columns
		vector<vector<SparseMatrixIndex>> sendIndices(m.RecvRanks.size());
		vector<vector<double>> sendValues(m.RecvRanks.size());
		for (size_t k = 0; k < m.RecvRanks.size(); k++)
		{
			for (BigNumber q = m.RecvStart[k]; q < m.RecvStart[k + 1]; q++)
				for (int r = 0; r < bs; r++)
				{
					BigNumber i = (m.NOwned + m.RecvIndices[q]) * bs + r;
					SparseMatrixIndex begin = localRows.outerIndexPtr()[i], end = localRows.outerIndexPtr()[i + 1];
					sendIndices[k].push_back(end - begin);
					for (SparseMatrixIndex p = begin; p < end; p++)
						sendIndices[k].push_back(globalColumn(localRows.innerIndexPtr()[p]));
					sendValues[k].insert(sendValues[k].end(), localRows.valuePtr() + begin, localRows.valuePtr() + end);
				}
		}
		vector<vector<SparseMatrixIndex>> recvIndices = m.ReverseExchange(sendIndices);
		vector<vector<double>> recvValues = m.ReverseExchange(sendValues);
		struct ReceivedRow { size_t Buffer; size_t Indices; size_t Values; SparseMatrixIndex Length; };
		vector<vector<ReceivedRow>> received(nOwnedRows);
		for (size_t k = 0; k < m.SendRanks.size(); k++)
		{
			size_t pi = 0, pv = 0;
			for (BigNumber s = m.SendStart[k]; s < m.SendStart[k + 1]; s++)
				for (int r = 0; r < bs; r++)
				{
					SparseMatrixIndex n = recvIndices[k][pi];
					received[m.SendIndices[s] * bs + r].push_back({ k, pi + 1, pv, n });
					pi += 1 + n;
					pv += n;
				}
		}

		// The column map: the ghosts of map in the owned rows, and the entities of other processes in the received rows
		vector<char> ghostUsed(m.NGhosts(), 0);
		#pragma omp parallel for
		for (BigNumber p = 0; p < (BigNumber)localRows.outerIndexPtr()[nOwnedRows]; p++)
		{
			SparseMatrixIndex c = localRows.innerIndexPtr()[p];
			if (c >= nOwnedCols)
				ghostUsed[c / bs - m.NOwned] = 1;
		}
		vector<BigNumber> ghosts;
		for (BigNumber g = 0; g < m.NGhosts(); g++)
			if (ghostUsed[g])
				ghosts.push_back(m.Ghosts[g]);
		BigNumber firstOwnedCol = m.Begin * bs, endOwnedCol = (m.Begin + m.NOwned) * bs;
		for (size_t k = 0; k < recvIndices.size(); k++)
		{
			size_t pi = 0;
			while (pi < recvIndices[k].size())
			{
				SparseMatrixIndex n = recvIndices[k][pi++];
				for (SparseMatrixIndex q = 0; q < n; q++, pi++)
				{
					SparseMatrixIndex c = recvIndices[k][pi];
					if (c < firstOwnedCol || c >= endOwnedCol)
						ghosts.push_back(c / bs);
				}
			}
		}
		auto owners = make_shared<const IndexMap>(m.WithoutGhosts());
		auto colMap = make_shared<const IndexMap>(m, std::move(ghosts)); // sorted (collective)
		const IndexMap& cm = *colMap;
		// The local columns of map and the global ones in the numbering of the blocks: Diag (< nOwnedCols) or Offd
		// (nOwnedCols + its column)
		vector<BigNumber> ghostPosition(m.NGhosts(), -1);
		for (BigNumber g = 0; g < m.NGhosts(); g++)
			if (ghostUsed[g])
				ghostPosition[g] = cm.Local(m.Ghosts[g]) - cm.NOwned;
		auto fromLocal = [&](SparseMatrixIndex c) -> SparseMatrixIndex
		{
			return c < nOwnedCols ? c : (SparseMatrixIndex)(nOwnedCols + ghostPosition[c / bs - m.NOwned] * bs + c % bs);
		};
		auto fromGlobal = [&](SparseMatrixIndex c) -> SparseMatrixIndex
		{
			if (c >= firstOwnedCol && c < endOwnedCol)
				return c - firstOwnedCol;
			return (SparseMatrixIndex)(nOwnedCols + (cm.Local(c / bs) - cm.NOwned) * bs + c % bs);
		};

		// The rows: the local ones plus the received ones, sorted (the duplicates summed), split
		ThreadLocal<vector<pair<SparseMatrixIndex, double>>> rowsTL;
		ThreadLocal<vector<pair<SparseMatrixIndex, double>>> offdTL;
		DistributedMatrix M;
		M.RowMap = owners;
		M.ColMap = colMap;
		M.RowBlockSize = bs;
		M.ColBlockSize = bs;
		auto rowEntries = [&](BigNumber i, vector<pair<SparseMatrixIndex, double>>& row)
		{
			row.clear();
			for (SparseMatrix::InnerIterator it(localRows, i); it; ++it)
				row.push_back({ fromLocal(it.col()), it.value() });
			for (const ReceivedRow& rr : received[i])
				for (SparseMatrixIndex q = 0; q < rr.Length; q++)
					row.push_back({ fromGlobal(recvIndices[rr.Buffer][rr.Indices + q]), recvValues[rr.Buffer][rr.Values + q] });
			auto byColumn = [](const pair<SparseMatrixIndex, double>& x, const pair<SparseMatrixIndex, double>& y) { return x.first < y.first; };
			if (!is_sorted(row.begin(), row.end(), byColumn))
				stable_sort(row.begin(), row.end(), byColumn);
		};
		auto addRange = [](const vector<pair<SparseMatrixIndex, double>>& row, SparseMatrixIndex from, SparseMatrixIndex to, SparseMatrixIndex shift, SparseMatrixOps::RowChunk& chunk)
		{
			for (size_t p = 0; p < row.size(); )
			{
				SparseMatrixIndex col = row[p].first;
				double v = row[p].second;
				for (p++; p < row.size() && row[p].first == col; p++)
					v += row[p].second;
				if (col >= from && col < to)
					chunk.Add(col - shift, v);
			}
		};
		SparseMatrixIndex nAllCols = nOwnedCols + cm.NGhosts() * bs;
		M.Diag = SparseMatrixOps::BuildByRows(nOwnedRows, nOwnedCols, 1, [&](BigNumber firstRow, BigNumber lastRow, SparseMatrixOps::RowChunk& chunk)
		{
			vector<pair<SparseMatrixIndex, double>>& row = rowsTL.Local();
			for (BigNumber i = firstRow; i < lastRow; ++i)
			{
				rowEntries(i, row);
				addRange(row, 0, nOwnedCols, 0, chunk);
				chunk.EndRow();
			}
		});
		M.Offd = SparseMatrixOps::BuildByRows(nOwnedRows, cm.NGhosts() * bs, 1, [&](BigNumber firstRow, BigNumber lastRow, SparseMatrixOps::RowChunk& chunk)
		{
			vector<pair<SparseMatrixIndex, double>>& row = offdTL.Local();
			for (BigNumber i = firstRow; i < lastRow; ++i)
			{
				rowEntries(i, row);
				addRange(row, nOwnedCols, nAllCols, nOwnedCols, chunk);
				chunk.EndRow();
			}
		});
		return M;
	}

	// M with its columns renumbered by blocks of blockSize: the column c * blockSize + k becomes newBlock(c) * blockSize + k,
	// nCols columns, and is dropped if newBlock(c) < 0; the rows sorted
	template <class NewBlock>
	static SparseMatrix RenumberColumns(const SparseMatrix& M, NewBlock newBlock, BigNumber nCols, int blockSize)
	{
		ThreadLocal<vector<pair<SparseMatrixIndex, double>>> rowsTL;
		return SparseMatrixOps::BuildByRows(M.rows(), nCols, 1, [&](BigNumber firstRow, BigNumber lastRow, SparseMatrixOps::RowChunk& chunk)
		{
			vector<pair<SparseMatrixIndex, double>>& row = rowsTL.Local();
			for (BigNumber i = firstRow; i < lastRow; ++i)
			{
				row.clear();
				for (SparseMatrix::InnerIterator it(M, i); it; ++it)
				{
					BigNumber block = newBlock(it.col() / blockSize);
					if (block >= 0)
						row.push_back({ (SparseMatrixIndex)(block * blockSize + it.col() % blockSize), it.value() });
				}
				auto byColumn = [](const pair<SparseMatrixIndex, double>& x, const pair<SparseMatrixIndex, double>& y) { return x.first < y.first; };
				if (!is_sorted(row.begin(), row.end(), byColumn))
					sort(row.begin(), row.end(), byColumn);
				for (const auto& coeff : row)
					chunk.Add(coeff.first, coeff.second);
				chunk.EndRow();
			}
		});
	}

private:
	// The row i of M appended to indices (its length, then its columns) and values
	static void PackRow(const SparseMatrix& M, BigNumber i, vector<SparseMatrixIndex>& indices, vector<double>& values)
	{
		SparseMatrixIndex begin = M.outerIndexPtr()[i], end = M.outerIndexPtr()[i + 1];
		indices.push_back(end - begin);
		indices.insert(indices.end(), M.innerIndexPtr() + begin, M.innerIndexPtr() + end);
		values.insert(values.end(), M.valuePtr() + begin, M.valuePtr() + end);
	}

	// The next row of indices and values (from positions pi and pv), as triplets of the row i
	static void UnpackRow(BigNumber i, const vector<SparseMatrixIndex>& indices, size_t& pi, const vector<double>& values, size_t& pv, vector<Eigen::Triplet<double, SparseMatrixIndex>>& triplets)
	{
		SparseMatrixIndex n = indices[pi++];
		for (SparseMatrixIndex q = 0; q < n; q++)
			triplets.emplace_back(i, indices[pi++], values[pv++]);
	}

public:
	//-------------------------------------------------------------------//
	//                            Products                             //
	//-------------------------------------------------------------------//

	// The values of the ghosts of ColMap for the owned values x (halo exchange). Collective on the halo.
	Vector GhostValues(const Vector& x) const
	{
		assert(x.rows() == ColMap->NOwned * ColBlockSize);
		return GhostValues(x.data());
	}

	// Same, x of size ColMap->NOwned * ColBlockSize
	Vector GhostValues(const double* x) const
	{
		Vector ghosts(ColMap->NGhosts() * ColBlockSize);
		ColMap->UpdateGhosts(x, ghosts.data(), ColBlockSize);
		return ghosts;
	}

	// y = A x (x: the owned values of the column entities). Collective on the halo.
	void Multiply(const Vector& x, Vector& y) const
	{
		Vector ghosts = GhostValues(x);
		y.noalias() = Diag * x;
		if (Offd.cols() > 0)
			y.noalias() += Offd * ghosts;
	}

	Vector operator*(const Vector& x) const
	{
		Vector y(Diag.rows());
		Multiply(x, y);
		return y;
	}

	// For MultiplyTransposed()
	void SetupTransposed()
	{
		_diagT = SparseMatrixOps::Transpose(Diag);
		_offdT = SparseMatrixOps::Transpose(Offd);
	}

	// y = A^T x (x: the owned rows; y: the owned values of the column entities): the contributions to the ghost columns
	// are added to their owners. SetupTransposed() first. Collective on the halo.
	Vector MultiplyTransposed(const Vector& x) const
	{
		assert(_diagT.rows() == Diag.cols() && _offdT.rows() == Offd.cols());
		Vector y = _diagT * x;
		Vector ghosts = _offdT * x;
		ColMap->AddGhostsToOwners(ghosts.data(), y.data(), ColBlockSize);
		return y;
	}
};
