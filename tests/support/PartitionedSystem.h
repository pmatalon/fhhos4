#pragma once
#include <memory>
#include <vector>
#include "Utils/DistributedMatrix.h"
#include "CartesianHybridSystem.h"

namespace fhhos4_tests
{
	// The cells of each process: horizontal strips, or blocks (px x py processes, px the largest divisor of the number
	// of processes not greater than its square root)
	enum class PartitionShape { Strips, Blocks };

	// The Cartesian system of n x n cells (CartesianHybridSystem) partitioned among the processes of comm,
	// renumbered so that each process owns contiguous cells and faces: the cells of the processes in rank order (in the
	// order of the system on each one), each face owned by the process of its first cell (the smaller rank of its cells).
	// Global: the renumbered system, on every process; the distributed blocks: the rows of the process.
	struct PartitionedSystem
	{
		Communicator Self; // a single process
		struct
		{
			SparseMatrix A_TT, A_TF, A_FT, A_FF, S;
		} Global;
		shared_ptr<const IndexMap> Cells, Faces; // the ownership, without ghosts
		SparseMatrix A_TT;                       // the cells of the process
		DistributedMatrix A_TF, S;

		PartitionedSystem(const CartesianHybridSystem& system, int n, PartitionShape shape, const Communicator& comm)
		{
			int P = comm.Size();
			int px = 1;
			for (int d = 1; d * d <= P; d++)
				if (P % d == 0)
					px = d;
			int py = P / px;
			BigNumber nCells = system.A_TT.rows(), nFaces = system.S.rows();
			std::vector<int> cellRank(nCells);
			for (int j = 0; j < n; j++)
				for (int i = 0; i < n; i++)
					cellRank[j * n + i] = shape == PartitionShape::Strips ? j * P / n : (j * py / n) * px + i * px / n;
			// The faces: their smaller rank among their cells (the columns of A_TF)
			std::vector<int> faceRank(nFaces, P);
			for (BigNumber T = 0; T < nCells; T++)
				for (SparseMatrix::InnerIterator it(system.A_TF, T); it; ++it)
					faceRank[it.col()] = std::min(faceRank[it.col()], cellRank[T]);

			std::vector<BigNumber> newCell = Renumbering(cellRank, P), newFace = Renumbering(faceRank, P);
			Global.A_TT = Permute(system.A_TT, newCell, newCell);
			Global.A_TF = Permute(system.A_TF, newCell, newFace);
			Global.A_FT = Permute(system.A_FT, newFace, newCell);
			Global.A_FF = Permute(system.A_FF, newFace, newFace);
			Global.S = Permute(system.S, newFace, newFace);

			BigNumber nOwnedCells = std::count(cellRank.begin(), cellRank.end(), comm.Rank());
			BigNumber nOwnedFaces = std::count(faceRank.begin(), faceRank.end(), comm.Rank());
			Cells = std::make_shared<const IndexMap>(comm, nOwnedCells);
			Faces = std::make_shared<const IndexMap>(comm, nOwnedFaces);
			A_TT = Global.A_TT.block(Cells->Begin, Cells->Begin, nOwnedCells, nOwnedCells);
			A_TF = DistributedMatrix::FromGlobalColumns(Cells, *Faces, SparseMatrix(Global.A_TF.middleRows(Cells->Begin, nOwnedCells)), 1, 1);
			S = DistributedMatrix::FromGlobalColumns(Faces, *Faces, SparseMatrix(Global.S.middleRows(Faces->Begin, nOwnedFaces)), 1, 1);
		}

	private:
		// The new numbers: by rank, then in the old order
		static std::vector<BigNumber> Renumbering(const std::vector<int>& rank, int P)
		{
			std::vector<BigNumber> next(P + 1, 0);
			for (int r : rank)
				next[r + 1]++;
			for (int r = 0; r < P; r++)
				next[r + 1] += next[r];
			std::vector<BigNumber> renumbering(rank.size());
			for (size_t i = 0; i < rank.size(); i++)
				renumbering[i] = next[rank[i]]++;
			return renumbering;
		}

		static SparseMatrix Permute(const SparseMatrix& M, const std::vector<BigNumber>& newRow, const std::vector<BigNumber>& newCol)
		{
			std::vector<Eigen::Triplet<double, SparseMatrixIndex>> t;
			for (BigNumber i = 0; i < M.rows(); i++)
				for (SparseMatrix::InnerIterator it(M, i); it; ++it)
					t.emplace_back(newRow[i], newCol[it.col()], it.value());
			SparseMatrix permuted(M.rows(), M.cols());
			permuted.setFromTriplets(t.begin(), t.end());
			return permuted;
		}
	};
}
