#pragma once
#include <vector>
#include "../../../Utils/Parallelism.h"
#include "../../../Utils/Types.h"
using namespace std;

// Lists of indices stored contiguously (as the rows of a CSR matrix): the faces of each cell, the cells of each
// face, etc. List i is Items[Start[i]], ..., Items[Start[i + 1] - 1].
struct IndexLists
{
	UninitializedVector<BigNumber> Start = { 0 };
	UninitializedVector<BigNumber> Items;

	BigNumber Size() const { return Start.size() - 1; }
	int Size(BigNumber i) const { return Start[i + 1] - Start[i]; }
	const BigNumber* operator[](BigNumber i) const { return Items.data() + Start[i]; }
	BigNumber* operator[](BigNumber i) { return Items.data() + Start[i]; }

	// n lists of the given sizes (sizes(i), called in parallel), items 0
	template <class Sizes>
	void Allocate(BigNumber n, Sizes sizes)
	{
		Start.resize(n + 1);
		#pragma omp parallel for if (n > 20000)
		for (BigNumber i = 0; i < n; ++i)
			Start[i + 1] = sizes(i);
		Start[0] = 0;
		for (BigNumber i = 0; i < n; ++i)
			Start[i + 1] += Start[i];
		ParallelAssign(Items, Start[n], (BigNumber)0);
	}

	// Position of item in list i, -1 if absent
	int Find(BigNumber i, BigNumber item) const
	{
		const BigNumber* list = (*this)[i];
		for (int k = 0; k < Size(i); ++k)
		{
			if (list[k] == item)
				return k;
		}
		return -1;
	}
};
