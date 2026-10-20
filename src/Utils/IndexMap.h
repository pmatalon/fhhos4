#pragma once
#include <algorithm>
#include <memory>
#include "Communicator.h"
#include "Utils.h"
using namespace std;

// The distribution of a set of entities (cells, faces...) among the processes of a communicator: each process owns a
// contiguous range of their global numbers, [Begin, Begin + NOwned), the ranges following the ranks, and may refer to
// entities owned by others, its ghosts. Local numbers: the owned entities 0, ..., NOwned - 1 (global Begin + i), then
// the ghosts NOwned + g, in the order of Ghosts: increasing global numbers (hence grouped by owner, in rank order) for
// the maps built by fhhos4, any order for those of the caller of the library (InOrder()).
// The values attached to the entities (e.g. the DoFs of the faces) are stored by blocks of blockSize, in this order.
//
// The halo exchanges between the owners of the entities and the processes where they are ghosts: UpdateGhosts()
// (owners to ghosts), AddGhostsToOwners() (ghosts to owners, summed). Their communication pattern (the halo) is built
// by the constructor, collectively.
class IndexMap
{
public:
	const Communicator* Comm = nullptr;
	BigNumber NOwned = 0;
	BigNumber Begin = 0;
	BigNumber GlobalSize = 0;
	vector<BigNumber> Offsets; // the first global number of each process, then GlobalSize (Size() + 1 values)
	vector<BigNumber> Ghosts;  // global numbers, in the order of the local numbers

	// Halo. Receive: the ghosts RecvIndices[RecvStart[k]], ..., RecvIndices[RecvStart[k+1] - 1] (positions in Ghosts, by
	// increasing global number) are owned by RecvRanks[k]. Send: the owned entities SendIndices[SendStart[k]], ...,
	// SendIndices[SendStart[k+1] - 1] (local numbers) are ghosts on SendRanks[k], in the order of its RecvIndices.
	vector<int> RecvRanks;
	vector<BigNumber> RecvStart;
	vector<BigNumber> RecvIndices;
	vector<int> SendRanks;
	vector<BigNumber> SendStart;
	vector<BigNumber> SendIndices;

private:
	// The positions in Ghosts by increasing global number (empty if Ghosts is increasing), for Local()
	vector<BigNumber> _sortedGhostPositions;

public:
	IndexMap() {}

	// nOwned entities on this process, no ghosts. Collective.
	IndexMap(const Communicator& comm, BigNumber nOwned)
		: Comm(&comm), NOwned(nOwned)
	{
		vector<BigNumber> counts = comm.AllGather(nOwned);
		Offsets.assign(comm.Size() + 1, 0);
		for (int r = 0; r < comm.Size(); r++)
			Offsets[r + 1] = Offsets[r] + counts[r];
		Begin = Offsets[comm.Rank()];
		GlobalSize = Offsets.back();
		RecvStart = { 0 };
		SendStart = { 0 };
	}

	// The ownership of owners, with the ghosts ghosts (global numbers of entities owned by other processes, in any
	// order, possibly repeated), sorted. Collective.
	IndexMap(const IndexMap& owners, vector<BigNumber> ghosts)
		: IndexMap(owners.WithoutGhosts())
	{
		sort(ghosts.begin(), ghosts.end());
		ghosts.erase(unique(ghosts.begin(), ghosts.end()), ghosts.end());
		SetGhosts(std::move(ghosts));
	}

	// The ownership of owners, with the ghosts in the given order (no repetition): the local numbers of the caller.
	// Collective.
	static IndexMap InOrder(const IndexMap& owners, vector<BigNumber> ghosts)
	{
		IndexMap map = owners.WithoutGhosts();
		map.SetGhosts(std::move(ghosts));
		return map;
	}

	// The same ownership, without ghosts (no communication)
	IndexMap WithoutGhosts() const
	{
		IndexMap m;
		m.Comm = Comm;
		m.NOwned = NOwned;
		m.Begin = Begin;
		m.GlobalSize = GlobalSize;
		m.Offsets = Offsets;
		m.RecvStart = { 0 };
		m.SendStart = { 0 };
		return m;
	}

	int Size() const { return Comm->Size(); }
	BigNumber NGhosts() const { return Ghosts.size(); }
	BigNumber NLocal() const { return NOwned + NGhosts(); }
	bool IsOwned(BigNumber global) const { return global >= Begin && global < Begin + NOwned; }

	// The process that owns the entity of global number global
	int Owner(BigNumber global) const
	{
		assert(global >= 0 && global < GlobalSize);
		return (int)(upper_bound(Offsets.begin(), Offsets.end(), global) - Offsets.begin()) - 1;
	}

	BigNumber Global(BigNumber local) const
	{
		return local < NOwned ? Begin + local : Ghosts[local - NOwned];
	}

	// The local number of the entity global: owned or ghost; -1 if neither
	BigNumber Local(BigNumber global) const
	{
		if (IsOwned(global))
			return global - Begin;
		if (_sortedGhostPositions.empty())
		{
			auto it = lower_bound(Ghosts.begin(), Ghosts.end(), global);
			if (it != Ghosts.end() && *it == global)
				return NOwned + (it - Ghosts.begin());
			return -1;
		}
		auto it = lower_bound(_sortedGhostPositions.begin(), _sortedGhostPositions.end(), global, [this](BigNumber position, BigNumber g) { return Ghosts[position] < g; });
		if (it != _sortedGhostPositions.end() && Ghosts[*it] == global)
			return NOwned + *it;
		return -1;
	}

	//-------------------------------------------------------------------//
	//                          Halo exchanges                           //
	//-------------------------------------------------------------------//

	// The values of the ghosts (blocks of blockSize), from those of their owners: ghosts[g*blockSize...] (NGhosts()
	// blocks) from owned (NOwned blocks) on the owner. Collective on the processes of the halo.
	template <class T>
	void UpdateGhosts(const T* owned, T* ghosts, int blockSize) const
	{
		if (Size() == 1)
			return;
		vector<vector<T>> sendBuffers(SendRanks.size()), recvBuffers(RecvRanks.size());
		for (size_t k = 0; k < SendRanks.size(); k++)
		{
			vector<T>& buffer = sendBuffers[k];
			buffer.resize((SendStart[k + 1] - SendStart[k]) * blockSize);
			T* p = buffer.data();
			for (BigNumber s = SendStart[k]; s < SendStart[k + 1]; s++)
				for (int i = 0; i < blockSize; i++)
					*p++ = owned[SendIndices[s] * blockSize + i];
		}
		for (size_t k = 0; k < RecvRanks.size(); k++)
			recvBuffers[k].resize((RecvStart[k + 1] - RecvStart[k]) * blockSize);
		Comm->Neighbours(SendRanks, sendBuffers, RecvRanks, recvBuffers);
		for (size_t k = 0; k < RecvRanks.size(); k++)
		{
			const T* p = recvBuffers[k].data();
			for (BigNumber q = RecvStart[k]; q < RecvStart[k + 1]; q++)
				for (int i = 0; i < blockSize; i++)
					ghosts[RecvIndices[q] * blockSize + i] = *p++;
		}
	}

	// Vector of the local values (owned, then ghosts) from the owned ones
	Vector WithGhosts(const Vector& owned, int blockSize) const
	{
		Vector local(NLocal() * blockSize);
		local.head(owned.rows()) = owned;
		UpdateGhosts(owned.data(), local.data() + NOwned * blockSize, blockSize);
		return local;
	}

	// The values of the ghosts added to those of their owners: owned[i*blockSize...] += the values ghosts[g*blockSize...]
	// of the entity i on the processes where it is a ghost (in rank order). Collective on the processes of the halo.
	template <class T>
	void AddGhostsToOwners(const T* ghosts, T* owned, int blockSize) const
	{
		if (Size() == 1)
			return;
		// The reverse of UpdateGhosts(): the receivers send, the senders receive
		vector<vector<T>> sendBuffers(RecvRanks.size()), recvBuffers(SendRanks.size());
		for (size_t k = 0; k < RecvRanks.size(); k++)
			sendBuffers[k] = GhostValuesOf(k, ghosts, blockSize);
		for (size_t k = 0; k < SendRanks.size(); k++)
			recvBuffers[k].resize((SendStart[k + 1] - SendStart[k]) * blockSize);
		Comm->Neighbours(RecvRanks, sendBuffers, SendRanks, recvBuffers, ReverseTag);
		for (size_t k = 0; k < SendRanks.size(); k++)
		{
			const T* p = recvBuffers[k].data();
			for (BigNumber s = SendStart[k]; s < SendStart[k + 1]; s++)
				for (int i = 0; i < blockSize; i++)
					owned[SendIndices[s] * blockSize + i] += *p++;
		}
	}

	// The values of the owners sent to the processes where they are ghosts, without summing: for each owned entity, the
	// values ghosts[g*blockSize...] it has on each of these processes, by process (the result of a reverse exchange):
	// received[k] from SendRanks[k], in the order of SendIndices[SendStart[k]...]
	template <class T>
	vector<vector<T>> GatherGhostsOnOwners(const T* ghosts, int blockSize) const
	{
		vector<vector<T>> sendBuffers(RecvRanks.size()), recvBuffers(SendRanks.size());
		if (Size() == 1)
			return recvBuffers;
		for (size_t k = 0; k < RecvRanks.size(); k++)
			sendBuffers[k] = GhostValuesOf(k, ghosts, blockSize);
		for (size_t k = 0; k < SendRanks.size(); k++)
			recvBuffers[k].resize((SendStart[k + 1] - SendStart[k]) * blockSize);
		Comm->Neighbours(RecvRanks, sendBuffers, SendRanks, recvBuffers, ReverseTag);
		return recvBuffers;
	}

	// Data of any size along the halo, forward (from the owners to the processes where their entities are ghosts):
	// send[k] goes to SendRanks[k]; returns received[k] from RecvRanks[k]. Collective on the processes of the halo.
	template <class T>
	vector<vector<T>> ForwardExchange(const vector<vector<T>>& send) const
	{
		assert(send.size() == SendRanks.size());
		return ExchangeOfAnySize(SendRanks, send, RecvRanks, ForwardDataTag);
	}

	// Same, reverse (from the processes where the entities are ghosts to their owners): send[k] goes to RecvRanks[k];
	// returns received[k] from SendRanks[k]
	template <class T>
	vector<vector<T>> ReverseExchange(const vector<vector<T>>& send) const
	{
		assert(send.size() == RecvRanks.size());
		return ExchangeOfAnySize(RecvRanks, send, SendRanks, ReverseDataTag);
	}

	// The values of the ghosts of RecvRanks[k], in the order of its halo
	template <class T>
	vector<T> GhostValuesOf(size_t k, const T* ghosts, int blockSize) const
	{
		vector<T> values;
		values.reserve((RecvStart[k + 1] - RecvStart[k]) * blockSize);
		for (BigNumber q = RecvStart[k]; q < RecvStart[k + 1]; q++)
			for (int i = 0; i < blockSize; i++)
				values.push_back(ghosts[RecvIndices[q] * blockSize + i]);
		return values;
	}

private:
	static constexpr int ReverseTag = 7103;
	static constexpr int SizeTag = 7104;
	static constexpr int ForwardDataTag = 7105;
	static constexpr int ReverseDataTag = 7106;

	template <class T>
	vector<vector<T>> ExchangeOfAnySize(const vector<int>& sendRanks, const vector<vector<T>>& send, const vector<int>& recvRanks, int tag) const
	{
		vector<vector<T>> received(recvRanks.size());
		if (Size() == 1)
			return received;
		vector<vector<int64_t>> sendSizes(sendRanks.size()), recvSizes(recvRanks.size(), vector<int64_t>(1));
		for (size_t k = 0; k < sendRanks.size(); k++)
			sendSizes[k] = { (int64_t)send[k].size() };
		Comm->Neighbours(sendRanks, sendSizes, recvRanks, recvSizes, SizeTag);
		for (size_t k = 0; k < recvRanks.size(); k++)
			received[k].resize(recvSizes[k][0]);
		Comm->Neighbours(sendRanks, send, recvRanks, received, tag);
		return received;
	}

	// Ghosts given (unique global numbers of entities owned by other processes), in their local order
	void SetGhosts(vector<BigNumber> ghosts)
	{
		Ghosts = std::move(ghosts);
		vector<BigNumber> sorted(Ghosts.size());
		for (BigNumber g = 0; g < NGhosts(); g++)
			sorted[g] = g;
		bool increasing = is_sorted(Ghosts.begin(), Ghosts.end());
		if (!increasing)
			sort(sorted.begin(), sorted.end(), [this](BigNumber a, BigNumber b) { return Ghosts[a] < Ghosts[b]; });
		bool invalid = false;
		for (BigNumber q = 0; q < NGhosts(); q++)
		{
			BigNumber global = Ghosts[sorted[q]];
			if (global < 0 || global >= GlobalSize || IsOwned(global) || (q > 0 && global == Ghosts[sorted[q - 1]]))
				invalid = true;
		}
		if (Comm->Any(invalid))
			Utils::FatalError("IndexMap: the ghosts must be global numbers of entities of other processes (global size " + to_string(GlobalSize) + "), without repetition.");
		if (!increasing)
			_sortedGhostPositions = sorted;
		BuildHalo(sorted);
	}

	// sorted: the positions of the ghosts by increasing global number
	void BuildHalo(const vector<BigNumber>& sorted)
	{
		RecvRanks.clear();
		RecvStart = { 0 };
		RecvIndices = sorted;
		for (BigNumber q = 0; q < NGhosts(); q++)
		{
			int owner = Owner(Ghosts[sorted[q]]);
			if (RecvRanks.empty() || RecvRanks.back() != owner)
			{
				if (q > 0)
					RecvStart.push_back(q);
				RecvRanks.push_back(owner);
			}
		}
		if (NGhosts() > 0)
			RecvStart.push_back(NGhosts());

		// The owners learn which of their entities are ghosts where
		vector<vector<BigNumber>> requests(Size());
		for (size_t k = 0; k < RecvRanks.size(); k++)
			for (BigNumber q = RecvStart[k]; q < RecvStart[k + 1]; q++)
				requests[RecvRanks[k]].push_back(Ghosts[RecvIndices[q]]);
		vector<vector<BigNumber>> requested = Comm->AllToAll(requests);
		SendRanks.clear();
		SendStart = { 0 };
		SendIndices.clear();
		for (int r = 0; r < Size(); r++)
		{
			if (requested[r].empty())
				continue;
			SendRanks.push_back(r);
			for (BigNumber global : requested[r])
			{
				if (!IsOwned(global))
					Utils::FatalError("IndexMap: a process asks for the entity " + to_string(global) + ", not owned by the process " + to_string(Comm->Rank()) + ".");
				SendIndices.push_back(global - Begin);
			}
			SendStart.push_back(SendIndices.size());
		}
	}
};
