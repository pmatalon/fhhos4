#pragma once
#include <memory>
#include "../Utils/DistributedMatrix.h"
using namespace std;

// On several processes (MPI): the hybrid system assembled by every process (the whole system), distributed by rows for
// the solver. The cells of each process are given by a partition (GMSH's); each face belongs to the process of the
// smaller rank among those of its cells. The cells and faces are renumbered so that each process owns contiguous ones,
// in rank order (in their order on each process): the global numbers of the distributed blocks.
class DistributedHybridSystem
{
public:
	const Communicator& Comm;
	int CellBlockSize;
	int FaceBlockSize;
	vector<BigNumber> NewCell, NewFace; // the new numbers of the cells and faces
	shared_ptr<const IndexMap> Cells;   // their ownership, without ghosts
	shared_ptr<const IndexMap> Faces;
	shared_ptr<const IndexMap> FaceMap; // the faces of the rows of the process: its own, then the others (by increasing new number)

	// cellRank, faceRank: the process of each cell and face, by number; A_T_F, A_F_F: the whole blocks, to find the faces
	// of the process's rows (the faces of its cells, and those coupled to its faces through the cells of other processes)
	DistributedHybridSystem(const Communicator& comm, const vector<int>& cellRank, const vector<int>& faceRank, const SparseMatrix& A_T_F, const SparseMatrix& A_F_F, int cellBlockSize, int faceBlockSize)
		: Comm(comm), CellBlockSize(cellBlockSize), FaceBlockSize(faceBlockSize)
	{
		NewCell = Renumbering(cellRank, comm.Size());
		NewFace = Renumbering(faceRank, comm.Size());
		Cells = make_shared<const IndexMap>(comm, (BigNumber)count(cellRank.begin(), cellRank.end(), comm.Rank()));
		Faces = make_shared<const IndexMap>(comm, (BigNumber)count(faceRank.begin(), faceRank.end(), comm.Rank()));
		vector<BigNumber> ghosts;
		for (BigNumber T = 0; T < (BigNumber)cellRank.size(); T++)
		{
			if (cellRank[T] != comm.Rank())
				continue;
			for (int k = 0; k < cellBlockSize; k++)
				for (SparseMatrix::InnerIterator it(A_T_F, T * cellBlockSize + k); it; ++it)
				{
					BigNumber f = NewFace[it.col() / faceBlockSize];
					if (!Faces->IsOwned(f))
						ghosts.push_back(f);
				}
		}
		for (BigNumber F = 0; F < (BigNumber)faceRank.size(); F++)
		{
			if (faceRank[F] != comm.Rank())
				continue;
			for (int k = 0; k < faceBlockSize; k++)
				for (SparseMatrix::InnerIterator it(A_F_F, F * faceBlockSize + k); it; ++it)
				{
					BigNumber f = NewFace[it.col() / faceBlockSize];
					if (!Faces->IsOwned(f))
						ghosts.push_back(f);
				}
		}
		FaceMap = make_shared<const IndexMap>(*Faces, std::move(ghosts));
	}

	// The rows of M of the cells (rowBlockSize = CellBlockSize) or faces of the process, in their new order, the columns
	// of the faces (colsAreFaces) or of the cells renumbered: global (local false), or in FaceMap (local true)
	SparseMatrix Rows(const SparseMatrix& M, bool rowsAreCells, bool colsAreFaces, bool local) const
	{
		const vector<BigNumber>& newRow = rowsAreCells ? NewCell : NewFace;
		const vector<BigNumber>& newCol = colsAreFaces ? NewFace : NewCell;
		const IndexMap& rowMap = rowsAreCells ? *Cells : *Faces;
		int rbs = rowsAreCells ? CellBlockSize : FaceBlockSize;
		int cbs = colsAreFaces ? FaceBlockSize : CellBlockSize;
		vector<BigNumber> oldRow(rowMap.NOwned);
		for (BigNumber i = 0; i < (BigNumber)newRow.size(); i++)
			if (rowMap.IsOwned(newRow[i]))
				oldRow[newRow[i] - rowMap.Begin] = i;
		BigNumber nCols = local ? FaceMap->NLocal() * cbs : (BigNumber)newCol.size() * cbs;
		ThreadLocal<vector<pair<SparseMatrixIndex, double>>> rowsTL;
		return SparseMatrixOps::BuildByRows(rowMap.NOwned * rbs, nCols, 1, [&](BigNumber firstRow, BigNumber lastRow, SparseMatrixOps::RowChunk& chunk)
		{
			vector<pair<SparseMatrixIndex, double>>& row = rowsTL.Local();
			for (BigNumber i = firstRow; i < lastRow; ++i)
			{
				row.clear();
				for (SparseMatrix::InnerIterator it(M, oldRow[i / rbs] * rbs + i % rbs); it; ++it)
				{
					BigNumber block = newCol[it.col() / cbs];
					if (local)
						block = FaceMap->Local(block);
					row.push_back({ (SparseMatrixIndex)(block * cbs + it.col() % cbs), it.value() });
				}
				sort(row.begin(), row.end(), [](const auto& x, const auto& y) { return x.first < y.first; });
				for (const auto& coeff : row)
					chunk.Add(coeff.first, coeff.second);
				chunk.EndRow();
			}
		});
	}

	// The diagonal blocks of the cells of the process of the block diagonal matrix M (cells x cells), numbered locally
	SparseMatrix CellBlocks(const SparseMatrix& M) const
	{
		SparseMatrix rows = Rows(M, true, false, false);
		BigNumber first = Cells->Begin * CellBlockSize;
		return DistributedMatrix::RenumberColumns(rows, [first, this](BigNumber c) { return c - first / CellBlockSize; }, Cells->NOwned * CellBlockSize, CellBlockSize);
	}

	// The values of v (by blocks of blockSize, numbered as the cells or the faces) of the cells or faces of the process
	Vector Own(const Vector& v, bool cells, int blockSize) const
	{
		const vector<BigNumber>& newNumber = cells ? NewCell : NewFace;
		const IndexMap& map = cells ? *Cells : *Faces;
		Vector own(map.NOwned * blockSize);
		for (BigNumber i = 0; i < (BigNumber)newNumber.size(); i++)
			if (map.IsOwned(newNumber[i]))
				own.segment((newNumber[i] - map.Begin) * blockSize, blockSize) = v.segment(i * blockSize, blockSize);
		return own;
	}

	// The whole vector, in the original numbering, on every process, from the values of the cells or faces of each
	// process. Collective.
	Vector Gather(const Vector& own, bool cells, int blockSize) const
	{
		const vector<BigNumber>& newNumber = cells ? NewCell : NewFace;
		vector<double> values(own.data(), own.data() + own.rows());
		vector<double> all = Comm.AllGatherV(values);
		Vector v(newNumber.size() * blockSize);
		for (BigNumber i = 0; i < (BigNumber)newNumber.size(); i++)
			v.segment(i * blockSize, blockSize) = Eigen::Map<const Vector>(all.data() + newNumber[i] * blockSize, blockSize);
		return v;
	}

	// The global numbers of the ghosts of FaceMap
	vector<int64_t> GhostFaces() const
	{
		return vector<int64_t>(FaceMap->Ghosts.begin(), FaceMap->Ghosts.end());
	}

private:
	// The new numbers: by rank, then in the old order
	static vector<BigNumber> Renumbering(const vector<int>& rank, int nRanks)
	{
		vector<BigNumber> next(nRanks + 1, 0);
		for (int r : rank)
			next[r + 1]++;
		for (int r = 0; r < nRanks; r++)
			next[r + 1] += next[r];
		vector<BigNumber> renumbering(rank.size());
		for (size_t i = 0; i < rank.size(); i++)
			renumbering[i] = next[rank[i]]++;
		return renumbering;
	}
};
