#pragma once
#include <algorithm>
#include "IndexLists.h"
#include "HybridAlgebraicMesh.h"
#include "LocalMatrices.h"
#include "../../FixedPoint/BlockDiagonalSolver.h"
#include "../../../Utils/Communicator.h"
#include "../../../Utils/SparseMatrixOps.h"
using namespace std;

// A coarsening pass of U-AMG on flat arrays: aggregation of the cells (pairwise, or with all their neighbours),
// interface collapsing, coarsening prolongation computed per aggregate from the local matrices of the operator
// (-coarsening-prolong 3 to 6, 10, -face-prolong 1, 2). Same aggregates, coarse faces and matrices as
// HybridAlgebraicMesh::Build(), Coarsen() and UncondensedLevel::Coarsen() with these options (same numbering, sums in
// the same order), except 10 (HarmonicExtensionOrInject), only implemented here. See UncondensedLevel::FlatPassesApply().
// Hint to load the memory at address into the cache
inline void Prefetch(const void* address)
{
#if defined(__GNUC__) || defined(__clang__)
	__builtin_prefetch(address);
#endif
}

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
	// On several processes (DistributedPass): the faces from NOwnedFaces are those of other processes, which may belong
	// to none of the cells of the pass (couplings of the operator); -1: all the faces are owned (one process)
	BigNumber NOwnedFaces = -1;
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
	// The two aggregates of each coarse face (CollapseInterfaces()): the first one, then the second one, greater (an
	// aggregate of another process on several processes, see DistributedPass), or -1 for the boundary of the first one
	vector<pair<BigNumber, BigNumber>> CoarseFaceAggregates;

	// The faces of a single cell (the boundary of the domain in the system, e.g. Neumann faces): false (default), removed
	// as the faces interior to an aggregate (their coarse values from the prolongation); true, kept and those of an
	// aggregate collapsed into one coarse face (see CollapseInterfaces())
	bool CollapseBoundaryFaces = false;

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

		// Cells of the faces, from the first column of each face (HybridAlgebraicMesh::Build()): the cells with a non-zero
		// in it, increasing
		auto firstColumnFaces = [&](BigNumber e, vector<BigNumber>& faces)
		{
			faces.clear();
			for (int k = 0; k < cbs; ++k)
				for (SparseMatrix::InnerIterator it(A_T_F, e * cbs + k); it; ++it)
					if (it.col() % fbs == 0 && (cbs == 1 || find(faces.begin(), faces.end(), it.col() / fbs) == faces.end()))
						faces.push_back(it.col() / fbs);
		};
		vector<int> nFaceCells(NFaces, 0);
		#pragma omp parallel
		{
			vector<BigNumber> faces;
			#pragma omp for
			for (BigNumber e = 0; e < NCells; ++e)
			{
				firstColumnFaces(e, faces);
				for (BigNumber f : faces)
				{
					#pragma omp atomic
					nFaceCells[f]++;
				}
			}
		}
		FaceCells.Allocate(NFaces, [&](BigNumber f) { return nFaceCells[f]; });
		fill(nFaceCells.begin(), nFaceCells.end(), 0);
		#pragma omp parallel
		{
			vector<BigNumber> faces;
			#pragma omp for
			for (BigNumber e = 0; e < NCells; ++e)
			{
				firstColumnFaces(e, faces);
				for (BigNumber f : faces)
				{
					int q;
					#pragma omp atomic capture
					q = nFaceCells[f]++;
					FaceCells[f][q] = e;
				}
			}
		}
		#pragma omp parallel for
		for (BigNumber f = 0; f < NFaces; ++f)
		{
			BigNumber* cells = FaceCells[f];
			if (FaceCells.Size(f) == 2)
			{
				if (cells[1] < cells[0])
					swap(cells[0], cells[1]);
			}
			else if (FaceCells.Size(f) > 2)
				sort(cells, cells + FaceCells.Size(f));
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
	//  Aggregation (PairwiseAggregation, AllNeighbourAggregation)  //
	//-------------------------------------------------------------//

	// Each cell not aggregated yet with its strongest available neighbour (pairwise aggregation), or with all its
	// neighbours not aggregated yet (allNeighbours). Returns true if the coarsest possible mesh is reached.
	bool Aggregate(bool allNeighbours = false)
	{
		// The data of each cell, together: the cells are visited in the order of the aggregation, which has no locality
		// in their numbering (each access to a cell is a cache miss), so they are also prefetched a few cells ahead
		struct CellData
		{
			BigNumber NeighboursStart;
			int NNeighbours;
			int NStrong;
			int Count;           // number of cells not aggregated yet of which the cell is a strong neighbour
			BigNumber Aggregate; // -1 if not aggregated yet
			bool Ordered;
		};
		const Communicator& comm = Communicator::Current(); // on several processes, each one aggregates its own cells
		if (NCells == 0)
		{
			NAggregates = 0;
			AggregateCells.Start.assign(1, 0);
			AggregateCells.Items.clear();
			CellAggregate.clear();
			return CoarsestMeshReached(comm);
		}
		vector<CellData> cells(NCells);
		#pragma omp parallel for
		for (BigNumber e = 0; e < NCells; ++e)
			cells[e] = { Neighbours.Start[e], NNeighbours[e], NStrongNeighbours[e], NStrongNeighbourOf[e], -1, false };
		const BigNumber* neighbours = Neighbours.Items.data();
		const int ahead = 8;
		// The cell order[i + 3*ahead], the neighbours of order[i + 2*ahead], their data for order[i + ahead]
		auto prefetch = [&](const vector<BigNumber>& order, BigNumber i, BigNumber n, bool strongOnly)
		{
			if (i + 3 * ahead < n)
				Prefetch(&cells[order[i + 3 * ahead]]);
			if (i + 2 * ahead < n)
				Prefetch(neighbours + cells[order[i + 2 * ahead]].NeighboursStart);
			if (i + ahead < n)
			{
				const CellData& c = cells[order[i + ahead]];
				for (int k = 0; k < (strongOnly ? c.NStrong : c.NNeighbours); ++k)
					Prefetch(&cells[neighbours[c.NeighboursStart + k]]);
			}
		};

		// Aggregation order: from a cell that is the strong neighbour of the fewest cells, the neighbours of the cells
		// already ordered, each time by increasing number of cells of which they are a strong neighbour
		vector<BigNumber> order(NCells);
		BigNumber start = 0;
		for (BigNumber e = 0; e < NCells; ++e)
		{
			if (cells[e].Count < cells[start].Count)
				start = e;
		}
		BigNumber nOrdered = 0;
		order[nOrdered++] = start;
		cells[start].Ordered = true;
		vector<BigNumber> unordered;
		BigNumber nextUnordered = 0;
		for (BigNumber i = 0; i < NCells; ++i)
		{
			if (i == nOrdered) // disconnected cells: continue from the first cell not ordered yet
			{
				while (cells[nextUnordered].Ordered)
					nextUnordered++;
				order[nOrdered++] = nextUnordered;
				cells[nextUnordered].Ordered = true;
			}
			prefetch(order, i, nOrdered, false);
			const CellData& c = cells[order[i]];
			unordered.clear();
			for (int k = 0; k < c.NNeighbours; ++k)
			{
				BigNumber n = neighbours[c.NeighboursStart + k];
				if (!cells[n].Ordered)
					unordered.push_back(n);
			}
			if (unordered.empty())
				continue;
			// Same sort as PairwiseAggregation: same order of the equal counts
			sort(unordered.begin(), unordered.end(), [&cells](BigNumber n1, BigNumber n2) { return cells[n1].Count < cells[n2].Count; });
			for (BigNumber n : unordered)
			{
				order[nOrdered++] = n;
				cells[n].Ordered = true;
			}
		}

		// Aggregation
		NAggregates = 0;
		AggregateCells.Start.assign(1, 0);
		AggregateCells.Start.reserve(NCells + 1);
		AggregateCells.Items.clear();
		AggregateCells.Items.reserve(NCells);
		auto addToAggregate = [&](BigNumber e)
		{
			AggregateCells.Items.push_back(e);
			CellData& c = cells[e];
			c.Aggregate = NAggregates;
			for (int k = 0; k < c.NStrong; ++k)
			{
				CellData& n = cells[neighbours[c.NeighboursStart + k]];
				if (n.Aggregate < 0)
					n.Count--;
			}
		};
		BigNumber remaining = NCells;
		for (BigNumber i = 0; i < NCells; ++i)
		{
			prefetch(order, i, NCells, !allNeighbours);
			BigNumber e = order[i];
			if (cells[e].Aggregate >= 0)
				continue;
			addToAggregate(e);
			remaining--;

			const CellData& c = cells[e];
			if (allNeighbours)
			{
				for (int k = 0; k < c.NNeighbours; ++k)
				{
					BigNumber n = neighbours[c.NeighboursStart + k];
					if (cells[n].Aggregate < 0)
					{
						addToAggregate(n);
						remaining--;
					}
				}
			}
			else
			{
				double min = 100000000;
				BigNumber chosen = -1;
				for (int k = 0; k < c.NStrong; ++k)
				{
					BigNumber n = neighbours[c.NeighboursStart + k];
					if (cells[n].Aggregate < 0 && cells[n].Count < min)
					{
						chosen = n;
						min = cells[n].Count;
					}
				}
				if (chosen >= 0)
				{
					addToAggregate(chosen);
					remaining--;
				}
			}
			NAggregates++;
			AggregateCells.Start.push_back(AggregateCells.Items.size());
			if (remaining == 0)
				break;
		}
		CellAggregate.resize(NCells);
		#pragma omp parallel for
		for (BigNumber e = 0; e < NCells; ++e)
			CellAggregate[e] = cells[e].Aggregate;

		return CoarsestMeshReached(comm);
	}

	// The coarsening ratio of the aggregation (of all the processes) is too small
	bool CoarsestMeshReached(const Communicator& comm) const
	{
		double coarseningRatio = (double)comm.Sum(NCells) / (double)comm.Sum(NAggregates);
		return coarseningRatio <= 1.01;
	}

	//-------------------------------------------------------------//
	//    Removed faces, interface collapsing (Coarsen())          //
	//-------------------------------------------------------------//

	// The coarse faces are the interfaces between two aggregates, and with CollapseBoundaryFaces, the boundary of an
	// aggregate in the system: its faces of a single cell (e.g. Neumann faces), collapsed into one coarse face of the
	// aggregate alone, as the interfaces (otherwise, they are removed). They are numbered by increasing first aggregate,
	// then second aggregate, the boundary after the others (the order of the interface collapsing of
	// HybridAlgebraicMesh::Coarsen()). Returns true if there are none.
	bool CollapseInterfaces()
	{
		// A face is removed if two cells of an aggregate share it, or if it belongs to a single cell (without
		// CollapseBoundaryFaces)
		FaceRemoved.assign(NFaces, 0);
		#pragma omp parallel for
		for (BigNumber f = 0; f < NFaces; ++f)
		{
			const BigNumber* cells = FaceCells[f];
			int n = FaceCells.Size(f);
			if (n == 1 && !CollapseBoundaryFaces)
				FaceRemoved[f] = 1;
			for (int i = 0; i < n && !FaceRemoved[f]; ++i)
				for (int j = i + 1; j < n; ++j)
					if (CellAggregate[cells[i]] == CellAggregate[cells[j]])
						FaceRemoved[f] = 1;
		}

		// Neighbouring aggregates of each aggregate (through its kept faces), increasing, and whether it has kept faces of a
		// single cell (its boundary, with CollapseBoundaryFaces)
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
		vector<int> nNeighbours(NAggregates), nUpper(NAggregates), nBoundary(NAggregates, 0);
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
					if (FaceCells.Size(f) == 1)
						nBoundary[k] = 1;
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
			firstCoarseFace[k + 1] = firstCoarseFace[k] + nUpper[k] + nBoundary[k];
		NCoarseFaces = firstCoarseFace[NAggregates];

		// Coarse face between the aggregates a and b, and of the boundary of the aggregate k
		auto coarseFace = [&](BigNumber a, BigNumber b)
		{
			BigNumber lo = min(a, b), hi = max(a, b);
			const BigNumber* upper = neighbours[lo] + nNeighbours[lo] - nUpper[lo];
			return firstCoarseFace[lo] + (lower_bound(upper, upper + nUpper[lo], hi) - upper);
		};
		auto boundaryCoarseFace = [&](BigNumber k) { return firstCoarseFace[k] + nUpper[k]; };
		AggregateCoarseFaces.Allocate(NAggregates, [&](BigNumber k) { return nNeighbours[k] + nBoundary[k]; });
		CoarseFaceAggregates.resize(NCoarseFaces);
		#pragma omp parallel for
		for (BigNumber k = 0; k < NAggregates; ++k)
		{
			for (int i = 0; i < nNeighbours[k]; ++i)
				AggregateCoarseFaces[k][i] = coarseFace(k, neighbours[k][i]);
			const BigNumber* upper = neighbours[k] + nNeighbours[k] - nUpper[k];
			for (int i = 0; i < nUpper[k]; ++i)
				CoarseFaceAggregates[firstCoarseFace[k] + i] = { k, upper[i] };
			if (nBoundary[k])
			{
				AggregateCoarseFaces[k][nNeighbours[k]] = boundaryCoarseFace(k); // the greatest of k
				CoarseFaceAggregates[boundaryCoarseFace(k)] = { k, -1 };
			}
		}
		FaceCoarseFace.assign(NFaces, -1);
		bool faceOfNoCell = false, faceOfMoreCells = false;
		#pragma omp parallel for reduction(||: faceOfNoCell, faceOfMoreCells)
		for (BigNumber f = 0; f < NFaces; ++f)
		{
			faceOfMoreCells = faceOfMoreCells || FaceCells.Size(f) > 2;
			if (FaceRemoved[f])
				continue;
			if (FaceCells.Size(f) >= 2)
				FaceCoarseFace[f] = coarseFace(CellAggregate[FaceCells[f][0]], CellAggregate[FaceCells[f][1]]);
			else if (FaceCells.Size(f) == 1)
				FaceCoarseFace[f] = boundaryCoarseFace(CellAggregate[FaceCells[f][0]]);
			else if (NOwnedFaces < 0 || f < NOwnedFaces) // a face of another process may belong to none of the cells
				faceOfNoCell = true;
		}
		// As in HybridAlgebraicMesh::Coarsen(), a face of no cell (no non-zero in its columns of A_T_F) would get no coarse
		// face (on several processes: the errors of all of them)
		const Communicator& comm = Communicator::Current();
		if (comm.Any(faceOfNoCell))
			Utils::FatalError("UncondensedAMG: a face belongs to no cell (its columns of A_T_F are zero): it has no coarse face.");
		// The blocks of a coarse level where the cells are coupled to the faces of non-neighbouring cells (e.g. with the
		// multigrid prolongation ReconstructionTrace, -prolong 1): not a mesh, which the interface collapsing assumes
		// (HybridAlgebraicMesh::Coarsen() left faces without coarse face, and crashed)
		if (comm.Any(faceOfMoreCells))
			Utils::FatalError("UncondensedAMG: a face belongs to more than two cells: the coarsening by interface collapsing is not possible on this level.");
		return comm.Sum(NCoarseFaces) == 0;
	}

	//-------------------------------------------------------------//
	//    Prolongation and blocks of the coarse cells              //
	//-------------------------------------------------------------//

	// The coarsening prolongation P of the pass, and the blocks A_T_Tc = Q_T^T A_T_T Q_T, A_T_Fc = Q_T^T A_T_F P, all
	// computed per aggregate (as UncondensedLevel::BuildProlongation() and Coarsen()). P is the injection on the kept
	// faces; on the removed faces, the face prolongation Q_F (FaceProlongation), Q_F smoothed by one block Jacobi
	// iteration (omega = 2/3) on the operator S, given by its local matrices on the cells
	// (FaceProlongationAndInteriorSmoothing), the trace of the reconstruction -A_TcTc^-1 A_TcFc in their aggregate
	// (ReconstructTraceOrInject), that trace smoothed by block Jacobi (ReconstructSmoothedTraceOrInject), or the local
	// solve of their rows of S, S_RR^-1 (-S_RK Y), where R are the removed faces of the aggregate, K the kept faces of its
	// cells and Y their injection (HarmonicExtensionOrInject): the discrete harmonic extension of the values of the kept
	// faces, i.e. FaceProlongationAndInteriorSmoothing with omega = 1 when the aggregate has one removed face (pairwise
	// aggregation). It converts no cell value into face values: it uses neither A_T_T (null) nor the blocks of the
	// coarse cells (A_T_Tc is left empty), and it doesn't depend on the bases of the cells (S is the condensed operator).
	// Q_F: the coarse face of a kept face; on a removed face, the average of the coarse faces of its aggregate
	// (BoundaryAggregatesInteriorAverage) or 0 (BoundaryAggregatesInteriorZero).
	void ComputeProlongation(UAMGProlongation prolongation, UAMGFaceProlongation faceProlongation, const SparseMatrix* A_T_T, const SparseMatrix& A_T_F,
	                         const LocalMatrices& S, SparseMatrix& P, SparseMatrix& A_T_Tc, SparseMatrix& A_T_Fc)
	{
		assert(prolongation == UAMGProlongation::FaceProlongation || prolongation == UAMGProlongation::FaceProlongationAndInteriorSmoothing
			|| prolongation == UAMGProlongation::ReconstructTraceOrInject || prolongation == UAMGProlongation::ReconstructSmoothedTraceOrInject
			|| prolongation == UAMGProlongation::HarmonicExtensionOrInject);
		assert(faceProlongation == UAMGFaceProlongation::BoundaryAggregatesInteriorAverage || faceProlongation == UAMGFaceProlongation::BoundaryAggregatesInteriorZero);
		bool harmonic = prolongation == UAMGProlongation::HarmonicExtensionOrInject;
		assert(harmonic || A_T_T);
		bool reconstruction = prolongation == UAMGProlongation::ReconstructTraceOrInject || prolongation == UAMGProlongation::ReconstructSmoothedTraceOrInject;
		bool smoothing = prolongation == UAMGProlongation::FaceProlongationAndInteriorSmoothing || prolongation == UAMGProlongation::ReconstructSmoothedTraceOrInject;
		bool averageOnRemovedFaces = faceProlongation == UAMGFaceProlongation::BoundaryAggregatesInteriorAverage;
		bool inverseDiagonal = BlockDiagonalSolver::DefaultMethod == DiagBlockSolveMethod::Inverse;
		int cbs = _cbs, fbs = _fbs;
		const double omega = 2.0 / 3.0;
		const double zero = NonZeroCoefficients::ZeroThreshold;
		DenseMatrix oneMinusOmegaIdentity = (1 - omega) * DenseMatrix::Identity(fbs, fbs);

		vector<double> blocksTT(harmonic ? 0 : (size_t)NAggregates * cbs * cbs); // the blocks of A_T_Tc, row-major
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
			vector<double> S_RR, B_RK;                    // HarmonicExtensionOrInject: S_RR and -S_RK Y, column-major
			Eigen::LLT<DenseMatrix> llt;
			Eigen::PartialPivLU<DenseMatrix> plu;
			Eigen::FullPivLU<DenseMatrix> lu;
			DenseMatrix A_TcTc, rhs, x, Df, inverseDf, A_ij, jacobiBlock;
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
			// The row of S on the face f: the sum of the rows of the local matrices of the cells of f, in increasing order
			// of the cells (LocalMatrices::Assemble()), their faces merged by increasing global number. Its faces in
			// w.jBlocks, the blocks S(f, g) in w.blockRow (fbs x fbs, column-major, one after the other).
			auto mergeRow = [&](BigNumber f)
			{
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

				// A_T_Tc: sum of the blocks of the cells (symmetric: of their lower triangular parts). Not with
				// HarmonicExtensionOrInject.
				double* blockTT = harmonic ? nullptr : blocksTT.data() + (size_t)K * cbs * cbs;
				if (!harmonic)
				{
					fill(blockTT, blockTT + cbs * cbs, 0.0);
					for (BigNumber e : w.cells)
					{
						for (int i = 0; i < cbs; ++i)
						{
							for (SparseMatrix::InnerIterator it(*A_T_T, e * cbs + i); it; ++it)
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
				}

				// Theta: row 0 of the reconstruction -A_TcTc^-1 A_TcFc in the aggregate, where A_TcFc = Q_T^T A_T_F Q_F
				double average = 1.0 / m; // Q_F on the removed faces (BoundaryAggregatesInteriorAverage)
				w.thetaCols.clear();
				w.theta.clear();
				if (!w.removed.empty() && reconstruction)
				{
					w.A_TcFc.assign(cbs * M, 0.0);
					for (int i = 0; i < cbs; ++i)
					{
						for (size_t p = w.qtaStart[i]; p < w.qtaStart[i + 1]; ++p)
						{
							BigNumber f = w.qta[p].first / fbs;
							int r = w.qta[p].first % fbs;
							double v = w.qta[p].second;
							if (!FaceRemoved[f])
								w.A_TcFc[(w.coarseFacePos[FaceCoarseFace[f]] * fbs + r) * cbs + i] += 1.0 * v;
							else if (averageOnRemovedFaces)
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

				// Rows of P on the removed faces
				if (harmonic && !w.removed.empty())
				{
					// S_RR X = -S_RK Y: the rows of S on the removed faces R of the aggregate, whose other faces are the kept
					// faces of its cells, injected (Y) into the coarse faces of the aggregate. The rows of the local
					// matrices of the cells of each removed face added directly into S_RR and -S_RK Y (no merge of the
					// rows of S: any order of the sums).
					int nR = w.removed.size();
					size_t n = (size_t)nR * fbs;
					w.S_RR.assign(n * n, 0.0);
					w.B_RK.assign(n * M, 0.0);
					for (int a = 0; a < nR; ++a)
					{
						BigNumber f = w.removed[a];
						for (int q = 0; q < FaceCells.Size(f); ++q)
						{
							BigNumber t = FaceCells[f][q];
							const BigNumber* tFaces = S.CellFaces(t);
							int nt = S.NFaces(t);
							int af = find(tFaces, tFaces + nt, f) - tFaces;
							assert(af < nt);
							size_t ld = (size_t)nt * fbs;
							// The block S_t(f, g) of the column-major local matrix: its coefficient (r, c) at block[c * ld + r]
							const double* rowF = S.Values(t) + (size_t)af * fbs;
							for (int b = 0; b < nt; ++b)
							{
								BigNumber g = tFaces[b];
								const double* block = rowF + (size_t)b * fbs * ld;
								if (FaceRemoved[g]) // a face of the cells of a removed face, removed: in the same aggregate
								{
									int bR = g == f ? a : find(w.removed.begin(), w.removed.end(), g) - w.removed.begin();
									assert(bR < nR);
									for (int c = 0; c < fbs; ++c)
										for (int r = 0; r < fbs; ++r)
											w.S_RR[(bR * fbs + c) * n + a * fbs + r] += block[c * ld + r];
								}
								else
								{
									size_t L = (size_t)w.coarseFacePos[FaceCoarseFace[g]] * fbs;
									for (int c = 0; c < fbs; ++c)
										for (int r = 0; r < fbs; ++r)
											w.B_RK[(L + c) * n + a * fbs + r] -= block[c * ld + r];
								}
							}
						}
					}
					// X, in B_RK. One removed face of 1 DoF (pairwise aggregation at k = 0): a division, 0 if the diagonal
					// coefficient d is 0, as the block Jacobi row of FaceProlongationAndInteriorSmoothing. Otherwise, a
					// rank-revealing LU (S_RR is small: the removed faces of an aggregate).
					if (n == 1)
					{
						double d = w.S_RR[0];
						double inverse = 1.0 / d;
						for (size_t L = 0; L < M; ++L)
							w.B_RK[L] = d == 0 ? 0.0 : w.B_RK[L] * inverse;
					}
					else
					{
						Eigen::Map<DenseMatrix> S_RR(w.S_RR.data(), n, n), B_RK(w.B_RK.data(), n, M);
						w.lu.compute(S_RR);
						B_RK = w.lu.solve(B_RK).eval();
					}
					for (int a = 0; a < nR; ++a)
					{
						BigNumber f = w.removed[a];
						removedRowsThread[f] = thread;
						for (int r = 0; r < fbs; ++r)
						{
							removedRowStart[f * fbs + r] = rr.Cols.size();
							for (size_t L = 0; L < M; ++L)
							{
								double v = w.B_RK[L * n + a * fbs + r];
								if (abs(v) > zero)
								{
									rr.Cols.push_back(cfs[L / fbs] * fbs + L % fbs);
									rr.Values.push_back(v);
								}
							}
							removedRowEnd[f * fbs + r] = rr.Cols.size();
						}
					}
				}
				else
				{
					for (BigNumber f : w.removed)
					{
						removedRowsThread[f] = thread;
						if (!smoothing)
						{
							// Q_F (FaceProlongation), or the trace of Theta on the constant mode (ReconstructTraceOrInject)
							for (int r = 0; r < fbs; ++r)
							{
								removedRowStart[f * fbs + r] = rr.Cols.size();
								if (prolongation == UAMGProlongation::FaceProlongation)
								{
									if (averageOnRemovedFaces)
										for (int B = 0; B < m; ++B)
										{
											rr.Cols.push_back(cfs[B] * fbs + r);
											rr.Values.push_back(average);
										}
								}
								else if (r == 0)
								{
									rr.Cols.insert(rr.Cols.end(), w.thetaCols.begin(), w.thetaCols.end());
									rr.Values.insert(rr.Values.end(), w.theta.begin(), w.theta.end());
								}
								removedRowEnd[f * fbs + r] = rr.Cols.size();
							}
							continue;
						}

						// J Y, where J is the block Jacobi iteration matrix of the operator on the row of the face, and Y the
						// injection on the kept faces; on the removed ones, Q_F (FaceProlongationAndInteriorSmoothing) or the
						// trace of Theta (ReconstructSmoothedTraceOrInject)
						mergeRow(f);
						double* br = w.blockRow.data();
						int nBlocks = w.jBlocks.size();

						// Block Jacobi: J_f = [ -omega D_f^-1 S_fg ... (1-omega) I ... ] (BlockJacobi::IterationMatrix())
						int Bf = lower_bound(w.jBlocks.begin(), w.jBlocks.end(), f) - w.jBlocks.begin();
						if (fbs == 1)
						{
							// The arithmetic of BlockDiagonalSolver::SolveLU() with a 1x1 block d: s * (1/d), 0 if d = 0 (rank 0).
							// With the inverse of the block (DiagBlockSolveMethod::Inverse), (1/d) * s: the same.
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
							if (inverseDiagonal)
								w.inverseDf = w.lu.inverse();
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
									if (inverseDiagonal)
										w.jacobiBlock.noalias() = w.inverseDf * w.A_ij;
									else
										BlockDiagonalSolver::SolveLU(w.lu.matrixLU(), permP, permQ, rank, w.A_ij, w.jacobiBlock);
									w.jacobiBlock *= -omega;
									for (int c = 0; c < fbs; ++c)
										for (int r = 0; r < fbs; ++r)
											block[c * fbs + r] = w.jacobiBlock(r, c);
								}
							}
						}

						// J_f Y, summed in the order of the columns of J_f
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
								else if (reconstruction)
								{
									if (cg == 0)
										for (size_t x = 0; x < w.theta.size(); ++x)
											add(local(w.thetaCols[x]), w.theta[x] * j);
								}
								else if (averageOnRemovedFaces)
									for (int B = 0; B < m; ++B)
										add(B * fbs + cg, average * j);
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

		// A_T_Tc (empty with HarmonicExtensionOrInject)
		A_T_Tc = SparseMatrix();
		if (!harmonic)
		{
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
		}

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
				else if (FaceCoarseFace[f] >= 0) // not a face of another process without cell in the pass (empty row)
					chunk.Add(FaceCoarseFace[f] * fbs + i % fbs, 1.0);
				chunk.EndRow();
			}
		});
	}

	// The face prolongation Q_F: the coarse face of each kept face; on a removed face, the average of the coarse faces of
	// its aggregate (BoundaryAggregatesInteriorAverage) or 0 (BoundaryAggregatesInteriorZero)
	SparseMatrix Q_F(UAMGFaceProlongation faceProlongation) const
	{
		int fbs = _fbs;
		return SparseMatrixOps::BuildByRows(NFaces * fbs, NCoarseFaces * fbs, fbs, [&](BigNumber firstRow, BigNumber lastRow, SparseMatrixOps::RowChunk& chunk)
		{
			for (BigNumber i = firstRow; i < lastRow; ++i)
			{
				BigNumber f = i / fbs;
				if (!FaceRemoved[f])
					chunk.Add(FaceCoarseFace[f] * fbs + i % fbs, 1.0);
				else if (faceProlongation == UAMGFaceProlongation::BoundaryAggregatesInteriorAverage)
				{
					BigNumber K = CellAggregate[FaceCells[f][0]];
					double average = 1.0 / AggregateCoarseFaces.Size(K);
					for (int B = 0; B < AggregateCoarseFaces.Size(K); ++B)
						chunk.Add(AggregateCoarseFaces[K][B] * fbs + i % fbs, average);
				}
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
