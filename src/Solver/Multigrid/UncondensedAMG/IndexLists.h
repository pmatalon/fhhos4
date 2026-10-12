#pragma once
#include <vector>
#include "../../../Utils/Types.h"
using namespace std;

// Lists of indices stored contiguously (as the rows of a CSR matrix): the faces of each cell, the cells of each
// face, etc. List i is Items[Start[i]], ..., Items[Start[i + 1] - 1].
struct IndexLists
{
	vector<BigNumber> Start = { 0 };
	vector<BigNumber> Items;

	BigNumber Size() const { return Start.size() - 1; }
	int Size(BigNumber i) const { return Start[i + 1] - Start[i]; }
	const BigNumber* operator[](BigNumber i) const { return Items.data() + Start[i]; }
	BigNumber* operator[](BigNumber i) { return Items.data() + Start[i]; }

	// n lists of the given sizes (sizes(i)), items not initialized
	template <class Sizes>
	void Allocate(BigNumber n, Sizes sizes)
	{
		Start.assign(n + 1, 0);
		for (BigNumber i = 0; i < n; ++i)
			Start[i + 1] = Start[i] + sizes(i);
		Items.resize(Start[n]);
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
