#pragma once
#include <memory>
#include "../Utils/DistributedMatrix.h"
using namespace std;

// On several processes (MPI): the hybrid system of the part of the mesh of each process (its cells and the other cells
// of their faces: GMSHMesh), in its local numbering, whose rows of its cells and faces are complete (the others may be
// missing), distributed by rows for the solver. Each cell and face has a global number, those of each process
// contiguous, in rank order: the global numbers of the distributed blocks.
class DistributedHybridSystem
{
public:
	const Communicator& Comm;
	int CellBlockSize;
	int FaceBlockSize;
	vector<BigNumber> NewCell, NewFace; // the global numbers of the cells and faces of the system of the process
	shared_ptr<const IndexMap> Cells;   // their ownership, without ghosts
	shared_ptr<const IndexMap> Faces;
	shared_ptr<const IndexMap> FaceMap; // the faces of the rows of the process: its own, then the others (by increasing global number)

	// cellNumbers, faceNumbers: the global numbers of the cells and faces of the process, by local number, nOwnCells and
	// nOwnFaces of which are its own; A_T_F, A_F_F: the blocks (their rows of the process), to find the faces of the
	// process's rows (the faces of its cells, and those coupled to its faces through the cells of other processes)
	DistributedHybridSystem(const Communicator& comm, vector<BigNumber> cellNumbers, BigNumber nOwnCells, vector<BigNumber> faceNumbers, BigNumber nOwnFaces,
	                        const SparseMatrix& A_T_F, const SparseMatrix& A_F_F, int cellBlockSize, int faceBlockSize)
		: Comm(comm), CellBlockSize(cellBlockSize), FaceBlockSize(faceBlockSize), NewCell(std::move(cellNumbers)), NewFace(std::move(faceNumbers))
	{
		Cells = make_shared<const IndexMap>(comm, nOwnCells);
		Faces = make_shared<const IndexMap>(comm, nOwnFaces);
		vector<BigNumber> ghosts;
		for (BigNumber T = 0; T < (BigNumber)NewCell.size(); T++)
		{
			if (!Cells->IsOwned(NewCell[T]))
				continue;
			for (int k = 0; k < cellBlockSize; k++)
				for (SparseMatrix::InnerIterator it(A_T_F, T * cellBlockSize + k); it; ++it)
				{
					BigNumber f = NewFace[it.col() / faceBlockSize];
					if (!Faces->IsOwned(f))
						ghosts.push_back(f);
				}
		}
		for (BigNumber F = 0; F < (BigNumber)NewFace.size(); F++)
		{
			if (!Faces->IsOwned(NewFace[F]))
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

	// The rows of M of the cells (rowBlockSize = CellBlockSize) or faces of the process, in their global order, the columns
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
		BigNumber nCols = (local ? FaceMap->NLocal() : (colsAreFaces ? Faces : Cells)->GlobalSize) * cbs;
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

	// The vector of the cells or faces of the process (by local number), from the values of its own ones: those of the
	// faces of its rows (FaceMap: among which the faces of its cells) from their processes, the others 0. Collective.
	Vector Local(const Vector& own, bool cells, int blockSize) const
	{
		const vector<BigNumber>& globalNumber = cells ? NewCell : NewFace;
		const IndexMap& map = cells ? *Cells : *FaceMap;
		Vector withGhosts = cells ? own : FaceMap->WithGhosts(own, blockSize);
		Vector v = Vector::Zero(globalNumber.size() * blockSize);
		for (BigNumber i = 0; i < (BigNumber)globalNumber.size(); i++)
		{
			BigNumber j = map.Local(globalNumber[i]);
			if (j >= 0)
				v.segment(i * blockSize, blockSize) = withGhosts.segment(j * blockSize, blockSize);
		}
		return v;
	}

	// The global numbers of the ghosts of FaceMap
	vector<int64_t> GhostFaces() const
	{
		return vector<int64_t>(FaceMap->Ghosts.begin(), FaceMap->Ghosts.end());
	}
};
