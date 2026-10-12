#pragma once
#include <algorithm>
#include "IndexLists.h"
#include "HybridAlgebraicMesh.h"
#include "LocalMatrices.h"
#include "../../FixedPoint/BlockDiagonalSolver.h"
#include "../../../Utils/SparseMatrixOps.h"
using namespace std;

// A coarsening pass of U-AMG with its default options, on flat arrays: pairwise aggregation of the cells, interface
// collapsing, coarsening prolongation ReconstructSmoothedTraceOrInject (-coarsening-prolong 6) computed per aggregate
// from the local matrices of the operator. Same aggregates, coarse faces and matrices as HybridAlgebraicMesh::Build(),
// Coarsen() and UncondensedLevel::Coarsen() with these options (same numbering, sums in the same order).
class CoarseningPass
{
private:
	int _cbs;
	int _fbs;
	bool _symmetric;
	double _strongCouplingThreshold;

public:
	BigNumber NCells = 0;
	BigNumber NFaces = 0;
	IndexLists CellFaces;            // faces of each cell, in the order of their first non-zero in the rows of A_T_F
	vector<double> FaceTraces;       // for each item of CellFaces, the trace of the block A_T_F(cell, face)
	vector<double> FaceToCellTraces; // non-symmetric matrix: for each item of CellFaces, the trace of A_F_T(face, cell)
	IndexLists FaceCells;            // cells of each face, increasing

	// Neighbours of each cell through its faces, by increasing coupling (the strongest first), and their couplings.
	// Only the first NNeighbours[e] items of the list of e are used. The strong neighbours of e are its first
	// NStrongNeighbours[e] neighbours.
	IndexLists Neighbours;
	vector<double> Couplings;
	vector<int> NNeighbours;
	vector<int> NStrongNeighbours;
	vector<int> NStrongNeighbourOf;  // number of cells of which the cell is a strong neighbour

	// Aggregation
	BigNumber NAggregates = 0;
	IndexLists AggregateCells;        // in the order of their aggregation
	vector<BigNumber> CellAggregate;
	vector<char> FaceRemoved;         // face interior to an aggregate
	BigNumber NCoarseFaces = 0;
	vector<BigNumber> FaceCoarseFace; // coarse face of each kept face (-1 for the removed faces)
	IndexLists AggregateCoarseFaces;  // coarse faces of each aggregate, increasing

	CoarseningPass(int cellBlockSize, int faceBlockSize, bool symmetric, double strongCouplingThreshold)
		: _cbs(cellBlockSize), _fbs(faceBlockSize), _symmetric(symmetric), _strongCouplingThreshold(strongCouplingThreshold)
	{}

	//-------------------------------------------------------------//
	//      Cells, faces and couplings (HybridAlgebraicMesh::Build) //
	//-------------------------------------------------------------//

	// A_F_T: non-symmetric matrix only (null otherwise)
	void Build(const SparseMatrix& A_T_F, const SparseMatrix* A_F_T)
	{
		int cbs = _cbs, fbs = _fbs;
		NCells = A_T_F.rows() / cbs;
		NFaces = A_T_F.cols() / fbs;

		// Faces of the cells (HybridAlgebraicMesh::BuildElementFaces())
		vector<int> nCellFaces(NCells);
		#pragma omp parallel
		{
			vector<BigNumber> faces;
			#pragma omp for
			for (BigNumber e = 0; e < NCells; ++e)
			{
				faces.clear();
				for (int k = 0; k < cbs; ++k)
				{
					for (SparseMatrix::InnerIterator it(A_T_F, e * cbs + k); it; ++it)
					{
						BigNumber f = it.col() / fbs;
						if (find(faces.rbegin(), faces.rend(), f) == faces.rend())
							faces.push_back(f);
					}
				}
				nCellFaces[e] = faces.size();
			}
		}
		CellFaces.Allocate(NCells, [&](BigNumber e) { return nCellFaces[e]; });
		FaceTraces.assign(CellFaces.Items.size(), 0);
		if (A_F_T)
			FaceToCellTraces.assign(CellFaces.Items.size(), 0);
		#pragma omp parallel for
		for (BigNumber e = 0; e < NCells; ++e)
		{
			BigNumber* faces = CellFaces[e];
			double* traces = FaceTraces.data() + CellFaces.Start[e];
			int n = 0;
			for (int k = 0; k < cbs; ++k)
			{
				for (SparseMatrix::InnerIterator it(A_T_F, e * cbs + k); it; ++it)
				{
					BigNumber f = it.col() / fbs;
					int a = n - 1;
					while (a >= 0 && faces[a] != f)
						a--;
					if (a < 0)
					{
						a = n;
						faces[n++] = f;
					}
					// Diagonal coefficient (k, k) of the block: summed in the order of the rows
					if (it.col() - f * fbs == k)
						traces[a] += it.value();
				}
			}
			// Non-symmetric matrix: the traces of the blocks A_F_T(face, e), summed in the order of the columns of e
			if (A_F_T)
			{
				double* traces2 = FaceToCellTraces.data() + CellFaces.Start[e];
				for (int k = 0; k < min(cbs, fbs); ++k)
				{
					for (int a = 0; a < n; ++a)
					{
						double v;
						if (Coefficient(*A_F_T, faces[a] * fbs + k, e * cbs + k, v))
							traces2[a] += v;
					}
				}
			}
		}

		// Cells of the faces, from the first column of each face (HybridAlgebraicMesh::Build())
		vector<int> nFaceCells(NFaces, 0);
		vector<BigNumber> last(NFaces, -1);
		for (BigNumber e = 0; e < NCells; ++e)
		{
			for (int k = 0; k < cbs; ++k)
			{
				for (SparseMatrix::InnerIterator it(A_T_F, e * cbs + k); it; ++it)
				{
					if (it.col() % fbs != 0)
						continue;
					BigNumber f = it.col() / fbs;
					if (last[f] != e)
					{
						last[f] = e;
						nFaceCells[f]++;
					}
				}
			}
		}
		FaceCells.Allocate(NFaces, [&](BigNumber f) { return nFaceCells[f]; });
		vector<BigNumber> pos(FaceCells.Start.begin(), FaceCells.Start.end() - 1);
		fill(last.begin(), last.end(), -1);
		for (BigNumber e = 0; e < NCells; ++e)
		{
			for (int k = 0; k < cbs; ++k)
			{
				for (SparseMatrix::InnerIterator it(A_T_F, e * cbs + k); it; ++it)
				{
					if (it.col() % fbs != 0)
						continue;
					BigNumber f = it.col() / fbs;
					if (last[f] != e)
					{
						last[f] = e;
						FaceCells.Items[pos[f]++] = e;
					}
				}
			}
		}

		// Neighbours and couplings
		Neighbours.Allocate(NCells, [&](BigNumber e)
		{
			int bound = 0;
			for (int a = 0; a < CellFaces.Size(e); ++a)
				bound += FaceCells.Size(CellFaces[e][a]);
			return bound;
		});
		Couplings.resize(Neighbours.Items.size());
		NNeighbours.assign(NCells, 0);
		NStrongNeighbours.assign(NCells, 0);
		NStrongNeighbourOf.assign(NCells, 0);
		#pragma omp parallel
		{
			vector<pair<BigNumber, double>> nb;
			#pragma omp for
			for (BigNumber e = 0; e < NCells; ++e)
			{
				nb.clear();
				const BigNumber* faces = CellFaces[e];
				for (int i = 0; i < CellFaces.Size(e); ++i)
				{
					BigNumber f = faces[i];
					for (int q = 0; q < FaceCells.Size(f); ++q)
					{
						BigNumber n = FaceCells[f][q];
						if (n == e)
							continue;
						int j = max(CellFaces.Find(n, f), 0);
						BigNumber ie = CellFaces.Start[e] + i, in = CellFaces.Start[n] + j;
						// Non-symmetric matrix: from the averages of the traces of the blocks A_T_F(cell, face) and
						// A_F_T(face, cell) of each cell (see HybridAlgebraicMesh::Build())
						double coupling = A_F_T ? HybridAlgebraicMesh::CouplingValue(0.5 * (FaceTraces[ie] + FaceToCellTraces[ie]), 0.5 * (FaceTraces[in] + FaceToCellTraces[in]))
						                        : HybridAlgebraicMesh::CouplingValue(FaceTraces[ie], FaceTraces[in]);
						// A neighbour sharing several faces is listed once, with the sum of the couplings
						auto it = find_if(nb.begin(), nb.end(), [n](const pair<BigNumber, double>& p) { return p.first == n; });
						if (it == nb.end())
							nb.push_back({ n, coupling });
						else
							it->second += coupling;
					}
				}
				// Same sort as HybridAlgebraicMesh::Build(): same order of the equal couplings
				sort(nb.begin(), nb.end(), [](const pair<BigNumber, double>& n1, const pair<BigNumber, double>& n2) { return n1.second < n2.second; });
				NNeighbours[e] = nb.size();
				for (int i = 0; i < (int)nb.size(); ++i)
				{
					Neighbours[e][i] = nb[i].first;
					Couplings[Neighbours.Start[e] + i] = nb[i].second;
				}
				for (int i = 0; i < (int)nb.size(); ++i)
				{
					if (!HybridAlgebraicMesh::IsStronglyCoupled(nb[i].second, nb[0].second, _strongCouplingThreshold))
						break;
					NStrongNeighbours[e]++;
					#pragma omp atomic
					NStrongNeighbourOf[nb[i].first]++;
				}
			}
		}
	}

	//-------------------------------------------------------------//
	//         Pairwise aggregation (PairwiseAggregation)           //
	//-------------------------------------------------------------//

	// Returns true if the coarsest possible mesh is reached
	bool Aggregate()
	{
		vector<int> count = NStrongNeighbourOf;
		CellAggregate.assign(NCells, -1);

		// Aggregation order: from a cell that is the strong neighbour of the fewest cells, the neighbours of the cells
		// already ordered, each time by increasing number of cells of which they are a strong neighbour
		vector<BigNumber> order(NCells);
		vector<char> ordered(NCells, 0);
		BigNumber start = 0;
		for (BigNumber e = 0; e < NCells; ++e)
		{
			if (count[e] < count[start])
				start = e;
		}
		BigNumber nOrdered = 0;
		order[nOrdered++] = start;
		ordered[start] = 1;
		vector<BigNumber> unordered;
		BigNumber nextUnordered = 0;
		for (BigNumber i = 0; i < NCells; ++i)
		{
			if (i == nOrdered) // disconnected cells: continue from the first cell not ordered yet
			{
				while (ordered[nextUnordered])
					nextUnordered++;
				order[nOrdered++] = nextUnordered;
				ordered[nextUnordered] = 1;
			}
			BigNumber e = order[i];
			unordered.clear();
			for (int k = 0; k < NNeighbours[e]; ++k)
			{
				BigNumber n = Neighbours[e][k];
				if (!ordered[n])
					unordered.push_back(n);
			}
			if (unordered.empty())
				continue;
			// Same sort as PairwiseAggregation: same order of the equal counts
			sort(unordered.begin(), unordered.end(), [&count](BigNumber n1, BigNumber n2) { return count[n1] < count[n2]; });
			for (BigNumber n : unordered)
			{
				order[nOrdered++] = n;
				ordered[n] = 1;
			}
		}

		// Aggregation: each cell not aggregated yet with its strongest available neighbour
		NAggregates = 0;
		AggregateCells.Start.assign(1, 0);
		AggregateCells.Start.reserve(NCells + 1);
		AggregateCells.Items.clear();
		AggregateCells.Items.reserve(NCells);
		auto addToAggregate = [&](BigNumber e)
		{
			AggregateCells.Items.push_back(e);
			CellAggregate[e] = NAggregates;
			for (int k = 0; k < NStrongNeighbours[e]; ++k)
			{
				BigNumber n = Neighbours[e][k];
				if (CellAggregate[n] < 0)
					count[n]--;
			}
		};
		BigNumber remaining = NCells;
		for (BigNumber e : order)
		{
			if (CellAggregate[e] >= 0)
				continue;
			addToAggregate(e);
			remaining--;

			double min = 100000000;
			BigNumber chosen = -1;
			for (int k = 0; k < NStrongNeighbours[e]; ++k)
			{
				BigNumber n = Neighbours[e][k];
				if (CellAggregate[n] < 0 && count[n] < min)
				{
					chosen = n;
					min = count[n];
				}
			}
			if (chosen >= 0)
			{
				addToAggregate(chosen);
				remaining--;
			}
			NAggregates++;
			AggregateCells.Start.push_back(AggregateCells.Items.size());
			if (remaining == 0)
				break;
		}

		double coarseningRatio = (double)NCells / (double)NAggregates;
		return coarseningRatio <= 1.01;
	}

	//-------------------------------------------------------------//
	//    Removed faces, interface collapsing (Coarsen())          //
	//-------------------------------------------------------------//

	// The coarse faces are the interfaces between two aggregates, numbered by increasing first aggregate, then second
	// aggregate (the order of the interface collapsing of HybridAlgebraicMesh::Coarsen()). Returns true if there are none.
	bool CollapseInterfaces()
	{
		// A face is removed if two cells of an aggregate share it
		FaceRemoved.assign(NFaces, 0);
		#pragma omp parallel for
		for (BigNumber f = 0; f < NFaces; ++f)
		{
			const BigNumber* cells = FaceCells[f];
			int n = FaceCells.Size(f);
			for (int i = 0; i < n && !FaceRemoved[f]; ++i)
				for (int j = i + 1; j < n; ++j)
					if (CellAggregate[cells[i]] == CellAggregate[cells[j]])
						FaceRemoved[f] = 1;
		}

		// Neighbouring aggregates of each aggregate (through its kept faces), increasing
		IndexLists neighbours;
		neighbours.Allocate(NAggregates, [&](BigNumber k)
		{
			int bound = 0;
			for (int c = 0; c < AggregateCells.Size(k); ++c)
			{
				BigNumber e = AggregateCells[k][c];
				for (int a = 0; a < CellFaces.Size(e); ++a)
					bound += FaceCells.Size(CellFaces[e][a]);
			}
			return bound;
		});
		vector<int> nNeighbours(NAggregates), nUpper(NAggregates);
		#pragma omp parallel for
		for (BigNumber k = 0; k < NAggregates; ++k)
		{
			BigNumber* list = neighbours[k];
			int n = 0;
			for (int c = 0; c < AggregateCells.Size(k); ++c)
			{
				BigNumber e = AggregateCells[k][c];
				for (int a = 0; a < CellFaces.Size(e); ++a)
				{
					BigNumber f = CellFaces[e][a];
					if (FaceRemoved[f])
						continue;
					for (int q = 0; q < FaceCells.Size(f); ++q)
					{
						BigNumber b = CellAggregate[FaceCells[f][q]];
						if (b != k)
							list[n++] = b;
					}
				}
			}
			sort(list, list + n);
			n = unique(list, list + n) - list;
			nNeighbours[k] = n;
			nUpper[k] = list + n - upper_bound(list, list + n, k);
		}
		vector<BigNumber> firstCoarseFace(NAggregates + 1, 0);
		for (BigNumber k = 0; k < NAggregates; ++k)
			firstCoarseFace[k + 1] = firstCoarseFace[k] + nUpper[k];
		NCoarseFaces = firstCoarseFace[NAggregates];

		// Coarse face between the aggregates a and b
		auto coarseFace = [&](BigNumber a, BigNumber b)
		{
			BigNumber lo = min(a, b), hi = max(a, b);
			const BigNumber* upper = neighbours[lo] + nNeighbours[lo] - nUpper[lo];
			return firstCoarseFace[lo] + (lower_bound(upper, upper + nUpper[lo], hi) - upper);
		};
		AggregateCoarseFaces.Allocate(NAggregates, [&](BigNumber k) { return nNeighbours[k]; });
		#pragma omp parallel for
		for (BigNumber k = 0; k < NAggregates; ++k)
		{
			for (int i = 0; i < nNeighbours[k]; ++i)
				AggregateCoarseFaces[k][i] = coarseFace(k, neighbours[k][i]);
		}
		FaceCoarseFace.assign(NFaces, -1);
		bool faceOfOneCell = false;
		#pragma omp parallel for reduction(||: faceOfOneCell)
		for (BigNumber f = 0; f < NFaces; ++f)
		{
			if (FaceRemoved[f])
				continue;
			if (FaceCells.Size(f) >= 2)
				FaceCoarseFace[f] = coarseFace(CellAggregate[FaceCells[f][0]], CellAggregate[FaceCells[f][1]]);
			else
				faceOfOneCell = true;
		}
		// As in HybridAlgebraicMesh::Coarsen(), a face of a single cell would get no coarse face
		if (faceOfOneCell)
			Utils::FatalError("UncondensedAMG: a face belongs to a single cell: it has no coarse face.");
		return NCoarseFaces == 0;
	}

	//-------------------------------------------------------------//
	//    Prolongation and blocks of the coarse cells              //
	//-------------------------------------------------------------//

	// The prolongation P of the pass (ReconstructSmoothedTraceOrInject: injection on the kept faces; on the removed
	// faces, the trace of the reconstruction -A_TcTc^-1 A_TcFc in their aggregate, smoothed by one block Jacobi
	// iteration (omega = 2/3) on the operator S, given by its local matrices on the cells), and the blocks
	// A_T_Tc = Q_T^T A_T_T Q_T, A_T_Fc = Q_T^T A_T_F P. All of them are computed per aggregate.
	void ComputeProlongation(const SparseMatrix& A_T_T, const SparseMatrix& A_T_F, const LocalMatrices& S, SparseMatrix& P, SparseMatrix& A_T_Tc, SparseMatrix& A_T_Fc)
	{
		int cbs = _cbs, fbs = _fbs;
		const double omega = 2.0 / 3.0;
		const double zero = NonZeroCoefficients::ZeroThreshold;
		DenseMatrix oneMinusOmegaIdentity = (1 - omega) * DenseMatrix::Identity(fbs, fbs);

		vector<double> blocksTT((size_t)NAggregates * cbs * cbs); // the blocks of A_T_Tc, row-major
		// Rows of P on the removed faces, stored by the thread that computes them
		struct RemovedRows
		{
			vector<SparseMatrixIndex> Cols;
			vector<double> Values;
		};
		ThreadLocal<RemovedRows> removedRowsTL;
		vector<int> removedRowsThread(NFaces, -1);
		vector<size_t> removedRowStart((size_t)NFaces * fbs + 1, 0); // and its end at removedRowEnd
		vector<size_t> removedRowEnd((size_t)NFaces * fbs, 0);

		struct Work
		{
			vector<pair<SparseMatrixIndex, double>> qta;  // rows of Q_T^T A_T_F of the aggregate (cbs rows, from qtaStart)
			vector<size_t> qtaStart;
			vector<SparseMatrixIndex> qtaPos;            // position of a column in the current row of Q_T^T A_T_F
			vector<int> coarseFacePos;                    // local number of a coarse face in the aggregate
			vector<BigNumber> cells, removed;
			vector<SparseMatrixIndex> thetaCols;          // row 0 of Theta of the aggregate
			vector<double> theta;
			vector<double> acc;
			vector<char> touched;
			vector<double> A_TcFc;                        // cbs x m*fbs, column-major
			vector<double> blockRow;                      // fbs x nBlocks*fbs, column-major
			vector<pair<BigNumber, int>> occurrences;     // (cell, local number of the removed face in it)
			vector<int> sorted, cursor;                   // local numbers of the faces of each occurrence, by global number
			vector<size_t> sortedStart;
			vector<BigNumber> jBlocks;
			Eigen::LLT<DenseMatrix> llt;
			Eigen::PartialPivLU<DenseMatrix> plu;
			Eigen::FullPivLU<DenseMatrix> lu;
			DenseMatrix A_TcTc, rhs, x, Df, A_ij, jacobiBlock;
		};
		ThreadLocal<Work> workTL;

		A_T_Fc = SparseMatrixOps::BuildByRows(NAggregates * cbs, NCoarseFaces * fbs, cbs, [&](BigNumber firstRow, BigNumber lastRow, SparseMatrixOps::RowChunk& chunk)
		{
			Work& w = workTL.Local();
			RemovedRows& rr = removedRowsTL.Local();
			int thread = Parallelism::ThreadNumber();
			if (w.coarseFacePos.empty())
			{
				w.qtaPos.assign(A_T_F.cols(), -1);
				w.coarseFacePos.assign(NCoarseFaces, -1);
				w.A_TcTc.resize(cbs, cbs);
				w.rhs.resize(cbs, fbs);
				w.x.resize(cbs, fbs);
				w.Df.resize(fbs, fbs);
				w.A_ij.resize(fbs, fbs);
				w.jacobiBlock.resize(fbs, fbs);
			}
			auto add = [&](size_t L, double v)
			{
				if (w.touched[L])
					w.acc[L] += v;
				else
				{
					w.acc[L] = v;
					w.touched[L] = 1;
				}
			};
			for (BigNumber K = firstRow / cbs; K < lastRow / cbs; ++K)
			{
				w.cells.assign(AggregateCells[K], AggregateCells[K] + AggregateCells.Size(K));
				sort(w.cells.begin(), w.cells.end());
				const BigNumber* cfs = AggregateCoarseFaces[K];
				int m = AggregateCoarseFaces.Size(K);
				size_t M = (size_t)m * fbs;
				for (int B = 0; B < m; ++B)
					w.coarseFacePos[cfs[B]] = B;
				auto local = [&](SparseMatrixIndex col) { return w.coarseFacePos[col / fbs] * fbs + col % fbs; };
				w.removed.clear();
				for (BigNumber e : w.cells)
					for (int a = 0; a < CellFaces.Size(e); ++a)
					{
						BigNumber f = CellFaces[e][a];
						if (FaceRemoved[f] && FaceCells[f][0] == e)
							w.removed.push_back(f);
					}

				// Q_T^T A_T_F: sum of the rows of the cells, in increasing order of the cells
				w.qta.clear();
				w.qtaStart.assign(1, 0);
				for (int i = 0; i < cbs; ++i)
				{
					if (w.cells.size() == 1)
					{
						for (SparseMatrix::InnerIterator it(A_T_F, w.cells[0] * cbs + i); it; ++it)
							w.qta.push_back({ it.col(), it.value() });
					}
					else if (w.cells.size() == 2) // merge of the two sorted rows
					{
						SparseMatrix::InnerIterator it1(A_T_F, w.cells[0] * cbs + i), it2(A_T_F, w.cells[1] * cbs + i);
						while (it1 || it2)
						{
							if (it1 && (!it2 || it1.col() < it2.col()))
							{
								w.qta.push_back({ it1.col(), it1.value() });
								++it1;
							}
							else if (!it1 || it2.col() < it1.col())
							{
								w.qta.push_back({ it2.col(), it2.value() });
								++it2;
							}
							else
							{
								w.qta.push_back({ it1.col(), it1.value() + it2.value() });
								++it1;
								++it2;
							}
						}
					}
					else
					{
						size_t rowStart = w.qta.size();
						for (BigNumber e : w.cells)
						{
							for (SparseMatrix::InnerIterator it(A_T_F, e * cbs + i); it; ++it)
							{
								SparseMatrixIndex& p = w.qtaPos[it.col()];
								if (p < 0)
								{
									p = w.qta.size() - rowStart;
									w.qta.push_back({ it.col(), it.value() });
								}
								else
									w.qta[rowStart + p].second += it.value();
							}
						}
						sort(w.qta.begin() + rowStart, w.qta.end(), [](const pair<SparseMatrixIndex, double>& x, const pair<SparseMatrixIndex, double>& y) { return x.first < y.first; });
						for (size_t p = rowStart; p < w.qta.size(); ++p)
							w.qtaPos[w.qta[p].first] = -1;
					}
					w.qtaStart.push_back(w.qta.size());
				}

				// A_T_Tc: sum of the blocks of the cells (symmetric: of their lower triangular parts)
				double* blockTT = blocksTT.data() + (size_t)K * cbs * cbs;
				fill(blockTT, blockTT + cbs * cbs, 0.0);
				for (BigNumber e : w.cells)
				{
					for (int i = 0; i < cbs; ++i)
					{
						for (SparseMatrix::InnerIterator it(A_T_T, e * cbs + i); it; ++it)
						{
							int j = it.col() - e * cbs;
							assert(j >= 0 && j < cbs);
							if (_symmetric && j > i)
								continue;
							blockTT[i * cbs + j] += it.value();
						}
					}
				}
				if (_symmetric)
					for (int i = 0; i < cbs; ++i)
						for (int j = i + 1; j < cbs; ++j)
							blockTT[i * cbs + j] = blockTT[j * cbs + i];

				// Theta: row 0 of the reconstruction -A_TcTc^-1 A_TcFc in the aggregate, where A_TcFc = Q_T^T A_T_F Q_F
				// (Q_F: the coarse face of a kept face, the average of the coarse faces of the aggregate on a removed face)
				w.thetaCols.clear();
				w.theta.clear();
				if (!w.removed.empty())
				{
					w.A_TcFc.assign(cbs * M, 0.0);
					double average = 1.0 / m;
					for (int i = 0; i < cbs; ++i)
					{
						for (size_t p = w.qtaStart[i]; p < w.qtaStart[i + 1]; ++p)
						{
							BigNumber f = w.qta[p].first / fbs;
							int r = w.qta[p].first % fbs;
							double v = w.qta[p].second;
							if (!FaceRemoved[f])
								w.A_TcFc[(w.coarseFacePos[FaceCoarseFace[f]] * fbs + r) * cbs + i] += 1.0 * v;
							else
								for (int B = 0; B < m; ++B)
									w.A_TcFc[(B * fbs + r) * cbs + i] += average * v;
						}
					}
					if (cbs == 1 && fbs == 1)
					{
						// The arithmetic of Eigen's solves with 1x1 matrices (multiplications by the inverse of the diagonal):
						// LLT, -(b / sqrt(a) / sqrt(a)); PartialPivLU, -(b / a)
						double a = blockTT[0];
						double inverse = 1.0 / a;
						if (_symmetric)
						{
							double sqrtA = a <= 0 ? a : sqrt(a); // LLT leaves a non-positive pivot as it is
							inverse = 1.0 / sqrtA;
						}
						for (int B = 0; B < m; ++B)
						{
							double x = w.A_TcFc[B] * inverse;
							if (_symmetric)
								x = (x - 0.0) * inverse;
							double t = -x;
							if (abs(t) > zero)
							{
								w.thetaCols.push_back(cfs[B]);
								w.theta.push_back(t);
							}
						}
					}
					else
					{
						for (int i = 0; i < cbs; ++i)
							for (int j = 0; j < cbs; ++j)
								w.A_TcTc(i, j) = blockTT[i * cbs + j];
						auto reconstruct = [&](const auto& solverA_TcTc)
						{
							for (int B = 0; B < m; ++B)
							{
								for (int j = 0; j < fbs; ++j)
									for (int i = 0; i < cbs; ++i)
										w.rhs(i, j) = w.A_TcFc[(B * fbs + j) * cbs + i];
								w.x = solverA_TcTc.solve(w.rhs);
								for (int j = 0; j < fbs; ++j)
								{
									double t = -w.x(0, j);
									if (abs(t) > zero)
									{
										w.thetaCols.push_back(cfs[B] * fbs + j);
										w.theta.push_back(t);
									}
								}
							}
						};
						if (_symmetric)
						{
							w.llt.compute(w.A_TcTc);
							reconstruct(w.llt);
						}
						else
						{
							w.plu.compute(w.A_TcTc);
							reconstruct(w.plu);
						}
					}
				}

				// Rows of P on the removed faces: J Y, where J is the block Jacobi iteration matrix of the operator on the
				// row of the face, and Y the injection on the kept faces, Theta on the removed ones
				for (BigNumber f : w.removed)
				{
					// Row of the operator: sum of the rows of the local matrices of the cells of f, in increasing order of the
					// cells (LocalMatrices::Assemble()): their faces merged by increasing global number
					w.occurrences.clear();
					w.sorted.clear();
					w.sortedStart.assign(1, 0);
					for (int q = 0; q < FaceCells.Size(f); ++q)
					{
						BigNumber t = FaceCells[f][q];
						const BigNumber* tFaces = S.CellFaces(t);
						int n = S.NFaces(t);
						int a = find(tFaces, tFaces + n, f) - tFaces;
						assert(a < n);
						w.occurrences.push_back({ t, a });
						size_t s0 = w.sorted.size();
						for (int b = 0; b < n; ++b) // insertion sort (the faces are usually sorted already)
						{
							size_t p = w.sorted.size();
							w.sorted.push_back(b);
							while (p > s0 && tFaces[w.sorted[p - 1]] > tFaces[b])
							{
								w.sorted[p] = w.sorted[p - 1];
								p--;
							}
							w.sorted[p] = b;
						}
						w.sortedStart.push_back(w.sorted.size());
					}
					int nOcc = w.occurrences.size();
					w.cursor.assign(nOcc, 0);
					w.jBlocks.clear();
					w.blockRow.resize(fbs * fbs * w.sorted.size());
					double* br = w.blockRow.data();
					while (true)
					{
						BigNumber g = 0;
						bool found = false;
						for (int q = 0; q < nOcc; ++q)
						{
							size_t pq = w.sortedStart[q] + w.cursor[q];
							if (pq < w.sortedStart[q + 1])
							{
								BigNumber gq = S.CellFaces(w.occurrences[q].first)[w.sorted[pq]];
								if (!found || gq < g)
									g = gq;
								found = true;
							}
						}
						if (!found)
							break;
						int B = w.jBlocks.size();
						w.jBlocks.push_back(g);
						bool first = true;
						for (int q = 0; q < nOcc; ++q)
						{
							size_t pq = w.sortedStart[q] + w.cursor[q];
							BigNumber t = w.occurrences[q].first;
							if (pq < w.sortedStart[q + 1] && S.CellFaces(t)[w.sorted[pq]] == g)
							{
								int a = w.occurrences[q].second;
								int b = w.sorted[pq];
								size_t n = (size_t)S.NFaces(t) * fbs;
								const double* St = S.Values(t);
								for (int c = 0; c < fbs; ++c)
									for (int r = 0; r < fbs; ++r)
									{
										double v = St[(b * fbs + c) * n + a * fbs + r];
										if (first)
											br[(B * fbs + c) * fbs + r] = v;
										else
											br[(B * fbs + c) * fbs + r] += v;
									}
								first = false;
								w.cursor[q]++;
							}
						}
					}
					int nBlocks = w.jBlocks.size();

					// Block Jacobi: J_f = [ -omega D_f^-1 S_fg ... (1-omega) I ... ] (BlockJacobi::IterationMatrix())
					int Bf = lower_bound(w.jBlocks.begin(), w.jBlocks.end(), f) - w.jBlocks.begin();
					if (fbs == 1)
					{
						// The arithmetic of BlockDiagonalSolver::SolveLU() with a 1x1 block d: s * (1/d), 0 if d = 0 (rank 0)
						double d = br[Bf];
						double inverse = 1.0 / d;
						for (int Bg = 0; Bg < nBlocks; ++Bg)
						{
							if (Bg == Bf)
								br[Bg] = oneMinusOmegaIdentity(0, 0);
							else
								br[Bg] = (d == 0 ? 0.0 : br[Bg] * inverse) * -omega;
						}
					}
					else
					{
						for (int c = 0; c < fbs; ++c)
							for (int r = 0; r < fbs; ++r)
								w.Df(r, c) = br[(Bf * fbs + c) * fbs + r];
						w.lu.compute(w.Df);
						const int* permP = w.lu.permutationP().indices().data();
						const int* permQ = w.lu.permutationQ().indices().data();
						int rank = w.lu.nonzeroPivots();
						for (int Bg = 0; Bg < nBlocks; ++Bg)
						{
							double* block = br + Bg * fbs * fbs;
							if (Bg == Bf)
							{
								for (int c = 0; c < fbs; ++c)
									for (int r = 0; r < fbs; ++r)
										block[c * fbs + r] = oneMinusOmegaIdentity(r, c);
							}
							else
							{
								for (int c = 0; c < fbs; ++c)
									for (int r = 0; r < fbs; ++r)
										w.A_ij(r, c) = block[c * fbs + r];
								BlockDiagonalSolver::SolveLU(w.lu.matrixLU(), permP, permQ, rank, w.A_ij, w.jacobiBlock);
								w.jacobiBlock *= -omega;
								for (int c = 0; c < fbs; ++c)
									for (int r = 0; r < fbs; ++r)
										block[c * fbs + r] = w.jacobiBlock(r, c);
							}
						}
					}

					// J_f Y, summed in the order of the columns of J_f
					removedRowsThread[f] = thread;
					for (int r = 0; r < fbs; ++r)
					{
						w.acc.assign(M, 0);
						w.touched.assign(M, 0);
						for (int c = 0; c < nBlocks * fbs; ++c)
						{
							double j = br[c * fbs + r];
							if (!(abs(j) > zero))
								continue;
							BigNumber g = w.jBlocks[c / fbs];
							int cg = c % fbs;
							if (!FaceRemoved[g])
								add(w.coarseFacePos[FaceCoarseFace[g]] * fbs + cg, 1.0 * j);
							else if (cg == 0)
								for (size_t x = 0; x < w.theta.size(); ++x)
									add(local(w.thetaCols[x]), w.theta[x] * j);
						}
						removedRowStart[f * fbs + r] = rr.Cols.size();
						for (size_t L = 0; L < M; ++L)
						{
							if (w.touched[L] && abs(w.acc[L]) > zero)
							{
								rr.Cols.push_back(cfs[L / fbs] * fbs + L % fbs);
								rr.Values.push_back(w.acc[L]);
							}
						}
						removedRowEnd[f * fbs + r] = rr.Cols.size();
					}
				}

				// A_T_Fc = (Q_T^T A_T_F) P, summed in the order of the columns of Q_T^T A_T_F
				for (int i = 0; i < cbs; ++i)
				{
					w.acc.assign(M, 0);
					w.touched.assign(M, 0);
					for (size_t p = w.qtaStart[i]; p < w.qtaStart[i + 1]; ++p)
					{
						BigNumber f = w.qta[p].first / fbs;
						int r = w.qta[p].first % fbs;
						double v = w.qta[p].second;
						if (!FaceRemoved[f])
							add(w.coarseFacePos[FaceCoarseFace[f]] * fbs + r, 1.0 * v);
						else
							for (size_t q = removedRowStart[f * fbs + r]; q < removedRowEnd[f * fbs + r]; ++q)
								add(local(rr.Cols[q]), rr.Values[q] * v);
					}
					for (size_t L = 0; L < M; ++L)
						if (w.touched[L])
							chunk.Add(cfs[L / fbs] * fbs + L % fbs, w.acc[L]);
					chunk.EndRow();
				}

				for (int B = 0; B < m; ++B)
					w.coarseFacePos[cfs[B]] = -1;
			}
		});

		// A_T_Tc
		A_T_Tc = SparseMatrix(NAggregates * cbs, NAggregates * cbs);
		A_T_Tc.resizeNonZeros((size_t)NAggregates * cbs * cbs);
		#pragma omp parallel for
		for (BigNumber i = 0; i < NAggregates * cbs; ++i)
		{
			BigNumber K = i / cbs;
			A_T_Tc.outerIndexPtr()[i + 1] = (i + 1) * cbs;
			for (int j = 0; j < cbs; ++j)
			{
				A_T_Tc.innerIndexPtr()[i * cbs + j] = K * cbs + j;
				A_T_Tc.valuePtr()[i * cbs + j] = blocksTT[(size_t)i * cbs + j];
			}
		}
		A_T_Tc.outerIndexPtr()[0] = 0;

		// P: injection on the kept faces
		vector<RemovedRows*> stores;
		for (RemovedRows& rr : removedRowsTL)
			stores.push_back(&rr);
		P = SparseMatrixOps::BuildByRows(NFaces * fbs, NCoarseFaces * fbs, fbs, [&](BigNumber firstRow, BigNumber lastRow, SparseMatrixOps::RowChunk& chunk)
		{
			for (BigNumber i = firstRow; i < lastRow; ++i)
			{
				BigNumber f = i / fbs;
				if (FaceRemoved[f])
				{
					const RemovedRows& rr = *stores[removedRowsThread[f]];
					for (size_t p = removedRowStart[i]; p < removedRowEnd[i]; ++p)
						chunk.Add(rr.Cols[p], rr.Values[p]);
				}
				else
					chunk.Add(FaceCoarseFace[f] * fbs + i % fbs, 1.0);
				chunk.EndRow();
			}
		});
	}

	// Q_T: the aggregate of each cell
	SparseMatrix Q_T() const
	{
		int cbs = _cbs;
		SparseMatrix Q(NCells * cbs, NAggregates * cbs);
		Q.resizeNonZeros(NCells * cbs);
		#pragma omp parallel for
		for (BigNumber i = 0; i < NCells * cbs; ++i)
		{
			Q.outerIndexPtr()[i + 1] = i + 1;
			Q.innerIndexPtr()[i] = CellAggregate[i / cbs] * cbs + i % cbs;
			Q.valuePtr()[i] = 1;
		}
		Q.outerIndexPtr()[0] = 0;
		return Q;
	}

private:
	// The coefficient (i, j) of M if it is stored
	static bool Coefficient(const SparseMatrix& M, BigNumber i, BigNumber j, double& value)
	{
		for (SparseMatrix::InnerIterator it(M, i); it; ++it)
		{
			if (it.col() == j)
			{
				value = it.value();
				return true;
			}
			if (it.col() > j)
				break;
		}
		return false;
	}
};
