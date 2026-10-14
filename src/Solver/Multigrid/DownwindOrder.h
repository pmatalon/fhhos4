#pragma once
#include <algorithm>
#include <limits>
#include <queue>
#include "../../Utils/Utils.h"
using namespace std;

// Downwind order of the block rows of a non-symmetric matrix (convection), for the Gauss-Seidel smoothers: a sweep in
// this order follows the flow, and nearly solves the convection-dominated part of the system, whereas in an order
// unrelated to the flow, a sweep carries the information only a few cells downstream (the iterations then grow with
// the number of cells along the flow). J. Bey, G. Wittum, "Downwind numbering: robust multigrid for convection-diffusion
// problems", Appl. Numer. Math. 23 (1997); W. Hackbusch, T. Probst, "Downwind Gauss-Seidel smoothing for convection
// dominated problems", Numer. Linear Algebra Appl. 4 (1997).
// The coupling of the block rows I and J is that of their first DoFs (the constant modes of the HHO face bases),
// rescaled by the interpolation of 1 (c_I c_J a_IJ, one coefficient per block row; c = 1 if interpOfOne is empty). J is upstream of I if their
// coupling is, in the row of I:
// - strong, as Notay's strong couplings (ETNA 37, 2010, Algorithm 2.1): -a_IJ >= threshold * max_K -a_IK;
// - one-sided: (-a_IJ) - (-a_JI) >= threshold * max_K -a_IK. Without this condition, a symmetric matrix whose rows have
//   different maxima (heterogeneous, anisotropic) would have directed couplings without any flow. With the threshold
//   0.25, it means a cell Péclet number >= 1/3 for the upwind finite differences in 1D.
// The order is topological (each block row after its upstream ones): the next block row is the available one of fewest
// unprocessed upstream block rows (none, except to break a cycle: recirculating flow), then of smallest number. It
// stays close to the natural order, which the hybrid smoothers need: each thread sweeps a contiguous part of the order
// (a breadth-first order, which puts all the block rows without upstream ones first, takes 30-40 iterations instead of
// 4 on HArDCore3D's convection-dominated case on 8 threads). The couplings are computed once (stored as the rows of
// the matrix), the upstream lists filled by contiguous ranges of rows per thread.
// Empty if no coupling is one-sided (symmetric or diffusion-dominated matrix): the natural order is then kept.
// nOneSided, nCycleBreaks: statistics, for the log.
inline vector<BigNumber> DownwindOrder(const SparseMatrix& A, int blockSize, const Vector& interpOfOne, double threshold, BigNumber& nOneSided, BigNumber& nCycleBreaks)
{
	BigNumber nb = A.rows() / blockSize;
	auto c = [&](BigNumber I) { return interpOfOne.rows() > 0 ? interpOfOne[I] : 1.0; };
	const SparseMatrixIndex* outer = A.outerIndexPtr();
	const SparseMatrixIndex* col = A.innerIndexPtr();
	const double* val = A.valuePtr();

	// Couplings -c_I c_J a_IJ of the first DoFs, I != J, stored as rows (by increasing J: the rows of A are sorted by
	// column), and the maximum of each row
	vector<size_t> start(nb + 1, 0);
	#pragma omp parallel for
	for (BigNumber I = 0; I < nb; ++I)
	{
		size_t n = 0;
		for (SparseMatrixIndex p = outer[I * blockSize]; p < outer[I * blockSize + 1]; ++p)
			n += col[p] % blockSize == 0 && col[p] != I * blockSize;
		start[I + 1] = n;
	}
	for (BigNumber I = 0; I < nb; ++I)
		start[I + 1] += start[I];
	vector<BigNumber> couplingCol(start[nb]);
	vector<double> coupling(start[nb]);
	vector<double> maxCoupling(nb, 0);
	#pragma omp parallel for
	for (BigNumber I = 0; I < nb; ++I)
	{
		size_t q = start[I];
		for (SparseMatrixIndex p = outer[I * blockSize]; p < outer[I * blockSize + 1]; ++p)
		{
			if (col[p] % blockSize != 0 || col[p] == I * blockSize)
				continue;
			BigNumber J = col[p] / blockSize;
			double v = -c(I) * c(J) * val[p];
			couplingCol[q] = J;
			coupling[q++] = v;
			maxCoupling[I] = max(maxCoupling[I], v);
		}
	}

	// The upstream block rows of each block row: each thread fills its contiguous range of rows, in order
	vector<int> nUpstream(nb, 0);
	vector<vector<BigNumber>> threadUpstream(Parallelism::NThreads());
	#pragma omp parallel
	{
		vector<BigNumber>& upstream = threadUpstream[Parallelism::ThreadNumber()];
		auto [begin, end] = Parallelism::ThreadChunk(nb);
		for (BigNumber I = begin; I < end; ++I)
		{
			double strong = threshold * maxCoupling[I];
			for (size_t q = start[I]; q < start[I + 1]; ++q)
			{
				double v = coupling[q];
				if (v <= 0 || v < strong)
					continue;
				BigNumber J = couplingCol[q];
				auto it = lower_bound(couplingCol.begin() + start[J], couplingCol.begin() + start[J + 1], I);
				double v_JI = it != couplingCol.begin() + start[J + 1] && *it == I ? coupling[it - couplingCol.begin()] : 0;
				if (v - v_JI >= strong)
				{
					upstream.push_back(J);
					nUpstream[I]++;
				}
			}
		}
	}
	nOneSided = 0;
	for (const vector<BigNumber>& upstream : threadUpstream)
		nOneSided += upstream.size();
	nCycleBreaks = 0;
	if (nOneSided == 0)
		return {};

	// The downstream block rows of each block row (CSR), from the upstream lists in the order of the rows
	vector<size_t> downStart(nb + 1, 0);
	for (const vector<BigNumber>& upstream : threadUpstream)
		for (BigNumber J : upstream)
			downStart[J + 1]++;
	for (BigNumber J = 0; J < nb; ++J)
		downStart[J + 1] += downStart[J];
	vector<BigNumber> downstream(downStart[nb]);
	{
		vector<size_t> pos(downStart.begin(), downStart.end() - 1);
		BigNumber I = 0;
		int remaining = 0; // upstream block rows of I not read yet
		for (const vector<BigNumber>& upstream : threadUpstream)
		{
			for (BigNumber J : upstream)
			{
				while (remaining == 0)
					remaining = nUpstream[I++];
				downstream[pos[J]++] = I - 1;
				remaining--;
			}
		}
	}

	// Topological order: the available block row (fewest unprocessed upstream ones, then smallest number) first. The
	// queue may hold outdated entries (the number of unprocessed upstream ones decreases): skipped.
	vector<BigNumber> order;
	order.reserve(nb);
	vector<bool> done(nb, false);
	vector<pair<int, BigNumber>> heap(nb);
	for (BigNumber I = 0; I < nb; ++I)
		heap[I] = { nUpstream[I], I };
	priority_queue<pair<int, BigNumber>, vector<pair<int, BigNumber>>, greater<pair<int, BigNumber>>> available(greater<pair<int, BigNumber>>(), std::move(heap));
	while (!available.empty())
	{
		auto [n, I] = available.top();
		available.pop();
		if (done[I] || n != nUpstream[I])
			continue;
		if (n > 0)
			nCycleBreaks++;
		done[I] = true;
		order.push_back(I);
		for (size_t q = downStart[I]; q < downStart[I + 1]; ++q)
		{
			BigNumber K = downstream[q];
			if (!done[K])
				available.push({ --nUpstream[K], K });
		}
	}
	assert(order.size() == nb);
	return order;
}
