#include <gtest/gtest.h>
#include "Utils/DistributedMatrix.h"

namespace
{
	const Communicator& World()
	{
		return Communicator::Current();
	}

	// The 1D Laplacian of size n, tridiagonal, scaled so that the values differ
	SparseMatrix Laplacian(BigNumber n)
	{
		vector<Eigen::Triplet<double, SparseMatrixIndex>> t;
		for (BigNumber i = 0; i < n; i++)
		{
			t.emplace_back(i, i, 2.0 + 0.01 * i);
			if (i > 0)
				t.emplace_back(i, i - 1, -1.0 - 0.001 * i);
			if (i < n - 1)
				t.emplace_back(i, i + 1, -1.0 + 0.002 * i);
		}
		SparseMatrix A(n, n);
		A.setFromTriplets(t.begin(), t.end());
		return A;
	}

	// The aggregation of the n entities by pairs, with blocks of bs rows and columns: n*bs x ((n+1)/2)*bs
	SparseMatrix PairwiseProlongation(BigNumber n, int bs)
	{
		vector<Eigen::Triplet<double, SparseMatrixIndex>> t;
		for (BigNumber i = 0; i < n; i++)
			for (int k = 0; k < bs; k++)
				t.emplace_back(i * bs + k, (i / 2) * bs + k, 1.0 + 0.1 * k + 0.01 * i);
		SparseMatrix P(n * bs, ((n + 1) / 2) * bs);
		P.setFromTriplets(t.begin(), t.end());
		return P;
	}

	// The rows [first, first + n) of M
	SparseMatrix Rows(const SparseMatrix& M, BigNumber first, BigNumber n)
	{
		return M.middleRows(first, n);
	}

	Vector TestVector(BigNumber n)
	{
		Vector x(n);
		for (BigNumber i = 0; i < n; i++)
			x[i] = 1.0 + 0.5 * sin((double)i);
		return x;
	}

	// The entities of each process: 3 + 2*rank, so that the sizes differ
	BigNumber NOwned()
	{
		return 3 + 2 * World().Rank();
	}
}

TEST(IndexMap, Ownership)
{
	IndexMap map(World(), NOwned());
	BigNumber expectedSize = 0;
	for (int r = 0; r < World().Size(); r++)
		expectedSize += 3 + 2 * r;
	EXPECT_EQ(map.GlobalSize, expectedSize);
	EXPECT_EQ(map.Begin, World().ExclusiveScan(NOwned()));
	for (BigNumber g = 0; g < map.GlobalSize; g++)
	{
		int owner = map.Owner(g);
		EXPECT_GE(g, map.Offsets[owner]);
		EXPECT_LT(g, map.Offsets[owner + 1]);
	}
}

TEST(IndexMap, HaloExchanges)
{
	IndexMap owners(World(), NOwned());
	// Ghosts: the first entity of every other process, and the last one of the next process
	vector<BigNumber> ghosts;
	for (int r = 0; r < World().Size(); r++)
	{
		if (r == World().Rank())
			continue;
		ghosts.push_back(owners.Offsets[r]);
		ghosts.push_back(owners.Offsets[r]); // repeated: removed
	}
	if (World().Rank() + 1 < World().Size())
		ghosts.push_back(owners.Offsets[World().Rank() + 2] - 1);
	IndexMap map(owners, ghosts);
	for (BigNumber g = 0; g < map.NGhosts(); g++)
	{
		EXPECT_FALSE(map.IsOwned(map.Ghosts[g]));
		EXPECT_EQ(map.Local(map.Ghosts[g]), map.NOwned + g);
		EXPECT_EQ(map.Global(map.NOwned + g), map.Ghosts[g]);
	}

	// Forward: the ghosts get the values of their owners (2 per entity)
	int bs = 2;
	vector<double> owned(map.NOwned * bs), ghostValues(map.NGhosts() * bs, -1);
	for (BigNumber i = 0; i < map.NOwned; i++)
		for (int k = 0; k < bs; k++)
			owned[i * bs + k] = 10.0 * (map.Begin + i) + k;
	map.UpdateGhosts(owned.data(), ghostValues.data(), bs);
	for (BigNumber g = 0; g < map.NGhosts(); g++)
		for (int k = 0; k < bs; k++)
			EXPECT_EQ(ghostValues[g * bs + k], 10.0 * map.Ghosts[g] + k);

	// Reverse: each owner gets the number of processes where the entity is a ghost
	vector<BigNumber> counts(map.NOwned, 0), ones(map.NGhosts(), 1);
	map.AddGhostsToOwners(ones.data(), counts.data(), 1);
	for (BigNumber i = 0; i < map.NOwned; i++)
	{
		BigNumber global = map.Begin + i;
		BigNumber expected = 0;
		if (i == 0)
			expected += World().Size() - 1;
		if (global == map.Begin + map.NOwned - 1 && World().Rank() > 0)
			expected += 1; // the last one, ghost on the previous process
		EXPECT_EQ(counts[i], expected) << "entity " << global;
	}
}

TEST(DistributedMatrix, Products)
{
	IndexMap owners(World(), NOwned());
	BigNumber n = owners.GlobalSize;
	auto rowMap = make_shared<const IndexMap>(owners);
	SparseMatrix A = Laplacian(n);
	DistributedMatrix dA = DistributedMatrix::FromGlobalColumns(rowMap, owners, Rows(A, owners.Begin, owners.NOwned), 1, 1);
	EXPECT_EQ(dA.GlobalNonZeros(), A.nonZeros());
	EXPECT_EQ(dA.GlobalRows(), n);

	Vector x = TestVector(n);
	Vector y = dA * x.segment(owners.Begin, owners.NOwned);
	Vector expected = (A * x).segment(owners.Begin, owners.NOwned);
	EXPECT_LE((y - expected).norm(), 1e-14 * expected.norm());

	SparseMatrix local = dA.LocalColumns();
	EXPECT_EQ(local.rows(), owners.NOwned);
	EXPECT_EQ(local.cols(), dA.ColMap->NLocal());
	EXPECT_EQ((SparseMatrix(dA.GlobalColumns() - Rows(A, owners.Begin, owners.NOwned))).norm(), 0);
}

TEST(DistributedMatrix, TransposedProductWithBlocks)
{
	// P: fine entities (rows) x coarse entities (columns), blocks of 3: the coarse entities split differently
	int bs = 3;
	IndexMap fine(World(), NOwned());
	BigNumber nCoarse = (fine.GlobalSize + 1) / 2;
	BigNumber nCoarseOwned = nCoarse / World().Size() + (World().Rank() < nCoarse % World().Size() ? 1 : 0);
	IndexMap coarse(World(), nCoarseOwned);
	ASSERT_EQ(coarse.GlobalSize, nCoarse);
	SparseMatrix P = PairwiseProlongation(fine.GlobalSize, bs);
	DistributedMatrix dP = DistributedMatrix::FromGlobalColumns(make_shared<const IndexMap>(fine), coarse, Rows(P, fine.Begin * bs, fine.NOwned * bs), bs, bs);
	dP.SetupTransposed();

	Vector xc = TestVector(nCoarse * bs);
	Vector yf = dP * xc.segment(coarse.Begin * bs, coarse.NOwned * bs);
	Vector expectedF = (P * xc).segment(fine.Begin * bs, fine.NOwned * bs);
	EXPECT_LE((yf - expectedF).norm(), 1e-14 * expectedF.norm());

	Vector xf = TestVector(fine.GlobalSize * bs);
	Vector yc = dP.MultiplyTransposed(xf.segment(fine.Begin * bs, fine.NOwned * bs));
	Vector expectedC = (P.transpose() * xf).segment(coarse.Begin * bs, coarse.NOwned * bs);
	EXPECT_LE((yc - expectedC).norm(), 1e-14 * expectedC.norm());
}

TEST(Communicator, Reductions)
{
	const Communicator& comm = World();
	int P = comm.Size();
	EXPECT_EQ(comm.Sum(1), P);
	EXPECT_EQ(comm.Max(comm.Rank()), P - 1);
	EXPECT_EQ(comm.Min(comm.Rank()), 0);
	EXPECT_EQ(comm.ExclusiveScan((BigNumber)1), comm.Rank());
	vector<int> local(comm.Rank() + 1, comm.Rank());
	vector<int> counts;
	vector<int> all = comm.AllGatherV(local, &counts);
	EXPECT_EQ((int)all.size(), P * (P + 1) / 2);
	for (int r = 0; r < P; r++)
		EXPECT_EQ(counts[r], r + 1);
	vector<vector<int>> send(P);
	for (int r = 0; r < P; r++)
		send[r] = vector<int>(r + 1, comm.Rank());
	vector<vector<int>> received = comm.AllToAll(send);
	for (int r = 0; r < P; r++)
	{
		EXPECT_EQ((int)received[r].size(), comm.Rank() + 1);
		for (int v : received[r])
			EXPECT_EQ(v, r);
	}
	Vector x = Vector::Constant(4, comm.Rank() + 1.0);
	double expected = 0;
	for (int r = 0; r < P; r++)
		expected += 4 * (r + 1.0) * (r + 1.0);
	EXPECT_DOUBLE_EQ(Parallel::SquaredNorm(x), expected);
	EXPECT_DOUBLE_EQ(Parallel::Norm(x), sqrt(expected));
}
