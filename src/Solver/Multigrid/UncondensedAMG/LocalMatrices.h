#pragma once
#include <algorithm>
#include <cstdlib>
#include <numeric>
#ifdef __linux__
#include <sys/mman.h>
#endif
#include "HybridAlgebraicMesh.h"
#include "../../../Utils/SparseMatrixOps.h"
using namespace std;

// Dense matrices on the faces of the cells of a mesh, one per cell, whose sum is a global sparse matrix:
// S = sum_T E_T^T S_T E_T, where E_T restricts a vector to the faces of the cell T (blocks of BlockSize unknowns).
//
// U-AMG computes its Galerkin products P^T S P with them: A_TT is block-diagonal, so the condensed operator S is a
// sum of local matrices, and P sends the faces of a cell only to the coarse faces of its aggregate. Hence
// P^T S P = sum_T P_T^T S_T P_T (P_T: rows of P on the faces of T), which is a sum of local matrices on the
// aggregates: the coarsening passes go on with them, and the operator of a level is assembled only once.
class LocalMatrices
{
private:
	// Values of all the matrices. The buffer is kept when the matrices are reallocated (it only grows): the first
	// access to new memory is slow (page faults, ~1 GB/s in WSL2), whereas the matrices of a coarsening pass take
	// about as much memory as those of the previous one.
	double* _values = nullptr;
	size_t _capacity = 0;
public:
	int BlockSize = 1;
	vector<size_t> FacesStart;  // the faces of cell t are Faces[FacesStart[t]], ..., Faces[FacesStart[t+1] - 1]
	vector<BigNumber> Faces;    // global numbers
	vector<size_t> ValuesStart; // the matrix of cell t is stored column-major from ValuesStart[t]

	LocalMatrices() {}
	LocalMatrices(const LocalMatrices&) = delete;
	LocalMatrices& operator=(const LocalMatrices&) = delete;
	~LocalMatrices() { std::free(_values); }

	void Swap(LocalMatrices& other)
	{
		swap(_values, other._values);
		swap(_capacity, other._capacity);
		swap(BlockSize, other.BlockSize);
		FacesStart.swap(other.FacesStart);
		Faces.swap(other.Faces);
		ValuesStart.swap(other.ValuesStart);
	}

	void Free()
	{
		std::free(_values);
		_values = nullptr;
		_capacity = 0;
		FacesStart = {};
		Faces = {};
		ValuesStart = {};
	}

	BigNumber NCells() const { return FacesStart.empty() ? 0 : FacesStart.size() - 1; }
	int NFaces(BigNumber t) const { return FacesStart[t + 1] - FacesStart[t]; }
	const BigNumber* CellFaces(BigNumber t) const { return Faces.data() + FacesStart[t]; }
	double* Values(BigNumber t) { return _values + ValuesStart[t]; }
	const double* Values(BigNumber t) const { return _values + ValuesStart[t]; }

	//-------------------------------------------//
	//               Construction                //
	//-------------------------------------------//

	// S = sum_T S_T on the cells of the mesh (whose faces' elements must be built): S_T(f, g) = S(f, g) / c, where c
	// is the number of cells containing both f and g (1 or 2: two cells share at most one face on simplicial and
	// Cartesian meshes, several on polygonal meshes built with -polymesh-fcs n). Exact for c <= 2. Only the lower
	// triangular part of S is read (S is symmetric).
	void Decompose(const SparseMatrix& S, const HybridAlgebraicMesh& mesh, int blockSize)
	{
		int bs = blockSize;
		BigNumber nCells = mesh.Elements.size();
		vector<int> nFaces(nCells);
		for (BigNumber t = 0; t < nCells; ++t)
			nFaces[t] = mesh.Elements[t].Faces.size();
		Allocate(bs, nCells, nFaces, [&](BigNumber t, BigNumber* faces)
		{
			for (int a = 0; a < nFaces[t]; ++a)
				faces[a] = mesh.Elements[t].Faces[a]->Number;
		});

		ThreadLocal<vector<int>> facePosTL((size_t)mesh.Faces.size(), -1); // local number of a face in the current cell
		#pragma omp parallel for schedule(dynamic, 256)
		for (BigNumber t = 0; t < nCells; ++t)
		{
			const HybridAlgebraicElement& T = mesh.Elements[t];
			vector<int>& facePos = facePosTL.Local();
			int n = T.Faces.size();
			size_t N = n * bs;
			for (int a = 0; a < n; ++a)
				facePos[T.Faces[a]->Number] = a;
			double* ST = Values(t);
			fill(ST, ST + N * N, 0.0);
			for (int a = 0; a < n; ++a)
			{
				const HybridAlgebraicFace& f = *T.Faces[a];
				for (int r = 0; r < bs; ++r)
				{
					BigNumber i = f.Number * bs + r;
					for (SparseMatrix::InnerIterator it(S, i); it; ++it)
					{
						BigNumber j = it.col();
						if (j > i)
							continue;
						int b = facePos[j / bs];
						if (b < 0) // face of another cell
							continue;
						int c = NCommonCells(f, mesh.Faces[j / bs], T);
						double v = c == 1 ? it.value() : it.value() / c;
						ST[(b * bs + j % bs) * N + a * bs + r] = v;
						ST[(a * bs + r) * N + b * bs + j % bs] = v;
					}
				}
			}
			for (int a = 0; a < n; ++a)
				facePos[T.Faces[a]->Number] = -1;
		}
	}

	// Local matrices of the aggregates of the mesh: S_K = sum_{T in K} P_T^T S_T P_T, where the S_T are cellMatrices
	// (the local matrices of the cells of the mesh) and P_T the rows of P on the faces of T. P_T is an injection on
	// the kept faces (rows of the identity), and its rows on the removed faces are read in P: the face prolongations
	// with interface collapsing (-face-prolong 1, 2) and the coarsening prolongations -coarsening-prolong 3 to 6
	// send the faces of a cell only to the coarse faces of its aggregate. Exactly symmetric (lower triangular part
	// mirrored).
	void GalerkinProduct(const LocalMatrices& cellMatrices, const HybridAlgebraicMesh& mesh, const SparseMatrix& P)
	{
		int bs = cellMatrices.BlockSize;
		BigNumber nAggregs = mesh.CoarseElements.size();
		vector<int> nCoarseFaces(nAggregs);
		int maxM = 0;
		for (BigNumber k = 0; k < nAggregs; ++k)
		{
			nCoarseFaces[k] = mesh.CoarseElements[k].CoarseFaces.size();
			maxM = max(maxM, nCoarseFaces[k]);
		}
		int maxN = 0;
		for (BigNumber t = 0; t < cellMatrices.NCells(); ++t)
			maxN = max(maxN, cellMatrices.NFaces(t));
		Allocate(bs, nAggregs, nCoarseFaces, [&](BigNumber k, BigNumber* faces)
		{
			for (int B = 0; B < nCoarseFaces[k]; ++B)
				faces[B] = mesh.CoarseElements[k].CoarseFaces[B]->Number;
		});

		struct Work
		{
			vector<int> coarseFacePos;     // local number of a coarse face in the current aggregate
			vector<int> removed;           // local numbers of the removed faces of the current cell
			vector<pair<int, int>> kept;   // local numbers of the kept faces of the current cell and of their coarse faces
			vector<double> Rt;             // transpose of the rows of P on the removed faces (m*bs x nRemoved*bs)
			vector<double> W;              // S_T P_T (n*bs x m*bs)
		};
		ThreadLocal<Work> workTL;
		#pragma omp parallel for schedule(dynamic, 64)
		for (BigNumber k = 0; k < nAggregs; ++k)
		{
			Work& w = workTL.Local();
			if (w.coarseFacePos.empty())
			{
				w.coarseFacePos.assign(mesh.CoarseFaces.size(), -1);
				w.Rt.resize((size_t)maxM * bs * maxN * bs);
				w.W.resize((size_t)maxN * bs * maxM * bs);
			}
			const HybridElementAggregate& K = mesh.CoarseElements[k];
			int m = nCoarseFaces[k];
			size_t M = m * bs;
			for (int B = 0; B < m; ++B)
				w.coarseFacePos[K.CoarseFaces[B]->Number] = B;
			double* SK = Values(k);
			fill(SK, SK + M * M, 0.0);
			double* Rt = w.Rt.data();
			double* W = w.W.data();

			for (const HybridAlgebraicElement* T : K.FineElements)
			{
				BigNumber t = T->Number;
				int n = cellMatrices.NFaces(t);
				size_t N = n * bs;
				const BigNumber* faces = cellMatrices.CellFaces(t);
				const double* ST = cellMatrices.Values(t);
				w.removed.clear();
				w.kept.clear();
				for (int a = 0; a < n; ++a)
				{
					const HybridAlgebraicFace& f = mesh.Faces[faces[a]];
					if (f.IsRemovedOnCoarseMesh)
						w.removed.push_back(a);
					else
					{
						assert(w.coarseFacePos[f.CoarseFace->Number] >= 0);
						w.kept.push_back({ a, w.coarseFacePos[f.CoarseFace->Number] });
					}
				}
				int nr = w.removed.size();

				// Rt
				fill(Rt, Rt + M * nr * bs, 0.0);
				for (int q = 0; q < nr; ++q)
				{
					for (int r = 0; r < bs; ++r)
					{
						for (SparseMatrix::InnerIterator it(P, faces[w.removed[q]] * bs + r); it; ++it)
						{
							int B = w.coarseFacePos[it.col() / bs];
							assert(B >= 0);
							Rt[(q * bs + r) * M + B * bs + it.col() % bs] = it.value();
						}
					}
				}

				// W = S_T P_T
				fill(W, W + N * M, 0.0);
				for (int q = 0; q < nr; ++q)
				{
					for (int kk = 0; kk < bs; ++kk)
					{
						const double* STcol = ST + (w.removed[q] * bs + kk) * N;
						const double* Rtcol = Rt + (q * bs + kk) * M;
						for (size_t j = 0; j < M; ++j)
						{
							double rj = Rtcol[j];
							if (rj == 0)
								continue;
							double* Wcol = W + j * N;
							for (size_t i = 0; i < N; ++i)
								Wcol[i] += STcol[i] * rj;
						}
					}
				}
				for (const pair<int, int>& ab : w.kept)
				{
					for (int c = 0; c < bs; ++c)
					{
						double* Wcol = W + (ab.second * bs + c) * N;
						const double* STcol = ST + (ab.first * bs + c) * N;
						for (size_t i = 0; i < N; ++i)
							Wcol[i] += STcol[i];
					}
				}

				// S_K += P_T^T W (lower triangular part)
				for (int q = 0; q < nr; ++q)
				{
					for (int kk = 0; kk < bs; ++kk)
					{
						const double* Rtcol = Rt + (q * bs + kk) * M;
						size_t row = w.removed[q] * bs + kk;
						for (size_t j = 0; j < M; ++j)
						{
							double wj = W[j * N + row];
							if (wj == 0)
								continue;
							double* SKcol = SK + j * M;
							for (size_t i = j; i < M; ++i)
								SKcol[i] += Rtcol[i] * wj;
						}
					}
				}
				for (const pair<int, int>& ab : w.kept)
				{
					for (int c = 0; c < bs; ++c)
					{
						size_t i = ab.second * bs + c;
						size_t row = ab.first * bs + c;
						for (size_t j = 0; j <= i; ++j)
							SK[j * M + i] += W[j * N + row];
					}
				}
			}

			for (size_t j = 0; j < M; ++j)
				for (size_t i = j + 1; i < M; ++i)
					SK[i * M + j] = SK[j * M + i];
			for (int B = 0; B < m; ++B)
				w.coarseFacePos[K.CoarseFaces[B]->Number] = -1;
		}
	}

	//-------------------------------------------//
	//                 Assembly                  //
	//-------------------------------------------//

	// Assembled matrix (nFaces faces), or only its block rows of the faces f such that (*rows)[f] (the others are
	// empty). Each row is the sum of the rows of the local matrices that contain it, in the order of the cells.
	// The rows are sorted by column; the explicit zeros are kept.
	SparseMatrix Assemble(BigNumber nFaces, const vector<bool>* rows = nullptr) const
	{
		int bs = BlockSize;
		BigNumber nCells = NCells();

		// Cells of each face, in increasing order, with the local number of the face in the cell
		vector<size_t> occStart(nFaces + 1, 0);
		for (BigNumber f : Faces)
			occStart[f + 1]++;
		for (BigNumber f = 0; f < nFaces; ++f)
			occStart[f + 1] += occStart[f];
		vector<pair<BigNumber, int>> occ(occStart[nFaces]);
		{
			vector<size_t> pos(occStart.begin(), occStart.end() - 1);
			for (BigNumber t = 0; t < nCells; ++t)
				for (int a = 0; a < NFaces(t); ++a)
					occ[pos[CellFaces(t)[a]]++] = { t, a };
		}

		// Local numbers of the faces of each cell, sorted by global number
		vector<int> sorted(Faces.size());
		#pragma omp parallel for
		for (BigNumber t = 0; t < nCells; ++t)
		{
			int* s = sorted.data() + FacesStart[t];
			const BigNumber* faces = CellFaces(t);
			iota(s, s + NFaces(t), 0);
			sort(s, s + NFaces(t), [faces](int a, int b) { return faces[a] < faces[b]; });
		}

		// Merge of the faces of the cells containing face f: calls block(g, q, b) for each face g, in increasing order,
		// and each cell q (number in the occurrences of f) where g has the local number b
		auto mergeFaces = [&](BigNumber f, vector<int>& cursor, auto block)
		{
			size_t nOcc = occStart[f + 1] - occStart[f];
			const pair<BigNumber, int>* o = occ.data() + occStart[f];
			cursor.assign(nOcc, 0);
			while (true)
			{
				BigNumber g = 0;
				bool found = false;
				for (size_t q = 0; q < nOcc; ++q)
				{
					BigNumber t = o[q].first;
					if (cursor[q] < NFaces(t))
					{
						BigNumber gq = CellFaces(t)[sorted[FacesStart[t] + cursor[q]]];
						if (!found || gq < g)
							g = gq;
						found = true;
					}
				}
				if (!found)
					break;
				for (size_t q = 0; q < nOcc; ++q)
				{
					BigNumber t = o[q].first;
					if (cursor[q] < NFaces(t) && CellFaces(t)[sorted[FacesStart[t] + cursor[q]]] == g)
					{
						block(g, q, sorted[FacesStart[t] + cursor[q]]);
						cursor[q]++;
					}
				}
			}
		};

		// Number of non-zeros of each row
		vector<SparseMatrixIndex> rowStart(nFaces * bs + 1, 0);
		#pragma omp parallel
		{
			vector<int> cursor;
			#pragma omp for schedule(dynamic, 1024)
			for (BigNumber f = 0; f < nFaces; ++f)
			{
				if (rows && !(*rows)[f])
					continue;
				SparseMatrixIndex nBlocks = 0;
				BigNumber lastG = 0;
				mergeFaces(f, cursor, [&](BigNumber g, size_t q, int b)
				{
					if (nBlocks == 0 || g != lastG)
						nBlocks++;
					lastG = g;
				});
				for (int r = 0; r < bs; ++r)
					rowStart[f * bs + r + 1] = nBlocks * bs;
			}
		}
		for (BigNumber i = 0; i < nFaces * bs; ++i)
			rowStart[i + 1] += rowStart[i];

		SparseMatrix M(nFaces * bs, nFaces * bs);
		M.resizeNonZeros(rowStart[nFaces * bs]);
		copy(rowStart.begin(), rowStart.end(), M.outerIndexPtr());
		SparseMatrixIndex* cols = M.innerIndexPtr();
		double* values = M.valuePtr();
		#pragma omp parallel
		{
			vector<int> cursor;
			#pragma omp for schedule(dynamic, 1024)
			for (BigNumber f = 0; f < nFaces; ++f)
			{
				if (rows && !(*rows)[f])
					continue;
				const pair<BigNumber, int>* o = occ.data() + occStart[f];
				SparseMatrixIndex rowLength = rowStart[f * bs + 1] - rowStart[f * bs];
				SparseMatrixIndex p = 0; // position of the current block in the row
				int nBlocks = 0;
				BigNumber lastG = 0;
				mergeFaces(f, cursor, [&](BigNumber g, size_t q, int b)
				{
					BigNumber t = o[q].first;
					int a = o[q].second;
					size_t n = NFaces(t) * bs;
					const double* St = Values(t);
					bool first = nBlocks == 0 || g != lastG;
					if (first)
					{
						p = nBlocks * bs;
						nBlocks++;
					}
					lastG = g;
					for (int r = 0; r < bs; ++r)
					{
						SparseMatrixIndex pos = rowStart[f * bs] + r * rowLength + p;
						for (int c = 0; c < bs; ++c)
						{
							double v = St[(b * bs + c) * n + a * bs + r];
							if (first)
							{
								cols[pos + c] = g * bs + c;
								values[pos + c] = v;
							}
							else
								values[pos + c] += v;
						}
					}
				});
			}
		}
		return M;
	}

private:
	// Allocates the matrices (not initialized); the faces of cell t are set by setFaces(t, faces) (called in parallel)
	template <class SetFaces>
	void Allocate(int blockSize, BigNumber nCells, const vector<int>& nFaces, SetFaces setFaces)
	{
		BlockSize = blockSize;
		FacesStart.assign(nCells + 1, 0);
		ValuesStart.assign(nCells + 1, 0);
		for (BigNumber t = 0; t < nCells; ++t)
		{
			FacesStart[t + 1] = FacesStart[t] + nFaces[t];
			ValuesStart[t + 1] = ValuesStart[t] + (size_t)nFaces[t] * nFaces[t] * blockSize * blockSize;
		}
		Faces.resize(FacesStart[nCells]);
		size_t size = ValuesStart[nCells];
		if (size > _capacity)
		{
			std::free(_values);
			_capacity = size + size / 4;
#ifdef __linux__
			size_t hugePage = 2 << 20;
			size_t bytes = (_capacity * sizeof(double) + hugePage - 1) / hugePage * hugePage;
			_values = (double*)aligned_alloc(hugePage, bytes);
			madvise(_values, bytes, MADV_HUGEPAGE); // the first access to the memory is then twice faster
#else
			_values = (double*)std::malloc(_capacity * sizeof(double));
#endif
		}
		#pragma omp parallel for
		for (BigNumber t = 0; t < nCells; ++t)
			setFaces(t, Faces.data() + FacesStart[t]);
	}

	// Number of cells containing the faces f and g of the cell T
	static int NCommonCells(const HybridAlgebraicFace& f, const HybridAlgebraicFace& g, const HybridAlgebraicElement& T)
	{
		if (&f == &g)
			return f.Elements.size();
		if (f.Elements.size() <= 2 && g.Elements.size() <= 2)
		{
			if (f.Elements.size() < 2 || g.Elements.size() < 2)
				return 1;
			const HybridAlgebraicElement* otherF = f.Elements[0] == &T ? f.Elements[1] : f.Elements[0];
			const HybridAlgebraicElement* otherG = g.Elements[0] == &T ? g.Elements[1] : g.Elements[0];
			return otherF == otherG ? 2 : 1;
		}
		int c = 0;
		for (const HybridAlgebraicElement* e : f.Elements)
			c += find(g.Elements.begin(), g.Elements.end(), e) != g.Elements.end();
		return c;
	}
};

// Local matrices of the operator of the current coarsening pass of U-AMG, shared by all the levels: the last
// coarsening pass of a level gives those of the next level's operator.
struct LocalOperator
{
	LocalMatrices Current;
	LocalMatrices Next;                     // buffer for those of the next coarsening pass
	const SparseMatrix* Assembled = nullptr; // the matrix equal to the sum of Current, if it has been assembled
};
