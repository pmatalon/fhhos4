#pragma once
#include <memory>
#include "CoarseningPass.h"
#include "../../../Utils/DistributedMatrix.h"
using namespace std;

// The parts of a coarsening pass of U-AMG (CoarseningPass) specific to several processes (UncondensedLevel::
// CoarsenInDistributedPasses()). Each process owns cells, which it aggregates alone (CoarseningPass::Aggregate()), and
// faces, numbered locally by FaceMap: its own faces, then the other faces of its cells (ghosts). The faces between
// cells of different processes are therefore kept, and collapsed with the others between the same two aggregates
// (CoarseningPass::CollapseInterfaces()), once the aggregates of the cells of the other processes are known. A coarse
// face between aggregates of two processes belongs to the process of the smaller aggregate (the smaller rank).
//
// The requirement on the faces: each face belongs to a process that owns one of its cells (e.g. the first one), and
// has at most two cells. The steps after CoarseningPass::Build(), Aggregate(), CollapseInterfaces(): FindRemoteCells(),
// AddRemoteAggregates(), NumberCoarseFaces(). All collective.
class DistributedPass
{
public:
	CoarseningPass& Pass;
	const IndexMap& FaceMap;
	IndexMap CellOwners;                     // the cells of the pass, all owned
	shared_ptr<const IndexMap> CellMap;      // the same, with the remote cells of the faces between processes as ghosts
	vector<BigNumber> FaceRemoteCell;        // for each local face, its cell on another process (local number in CellMap, >= NCells), -1 if none
	IndexMap AggregateOwners;                // after AddRemoteAggregates()
	vector<BigNumber> GhostAggregates;       // the aggregates of the remote cells (global numbers, increasing): numbered NAggregates + i in the pass
	shared_ptr<const IndexMap> CoarseFaceMap; // after NumberCoarseFaces(): the coarse faces of the aggregates of the process
	vector<BigNumber> CoarseFacePosition;    // the local number in CoarseFaceMap of each coarse face of the pass

	DistributedPass(CoarseningPass& pass, const IndexMap& faceMap)
		: Pass(pass), FaceMap(faceMap), CellOwners(*faceMap.Comm, pass.NCells)
	{
		assert(pass.NFaces == faceMap.NLocal());
		pass.NOwnedFaces = faceMap.NOwned;
	}

	const Communicator& Comm() const { return *FaceMap.Comm; }

	// After Pass.Build(), whose cells of the faces (FaceCells) are those of the process: the cells of the faces on other
	// processes (FaceRemoteCell, CellMap). The owner of a face gathers its cells, then sends them back.
	void FindRemoteCells()
	{
		const IndexMap& faces = FaceMap;
		BigNumber nOwned = faces.NOwned, nGhosts = faces.NGhosts();
		BigNumber cellBegin = CellOwners.Begin;
		const IndexLists& faceCells = Pass.FaceCells;
		bool invalid = false;

		// The ghost faces: their (single) cell of this process, if any (a ghost face may belong to no cell of the
		// process: a coupling of the operator, through a cell of another process)
		vector<BigNumber> ghostCell(nGhosts, -1);
		for (BigNumber g = 0; g < nGhosts; g++)
		{
			BigNumber f = nOwned + g;
			if (faceCells.Size(f) > 1)
				invalid = true; // both cells of a face of another process on this one
			else if (faceCells.Size(f) == 1)
				ghostCell[g] = cellBegin + faceCells[f][0];
		}
		// The owners gather the cells of their faces on the other processes
		vector<BigNumber> faceRemoteGlobal(faces.NLocal(), -1);
		vector<vector<BigNumber>> received = faces.GatherGhostsOnOwners(ghostCell.data(), 1);
		vector<int> nRemote(nOwned, 0);
		for (size_t k = 0; k < faces.SendRanks.size(); k++)
		{
			for (BigNumber s = faces.SendStart[k]; s < faces.SendStart[k + 1]; s++)
			{
				BigNumber f = faces.SendIndices[s];
				BigNumber cell = received[k][s - faces.SendStart[k]];
				if (cell < 0)
					continue;
				faceRemoteGlobal[f] = cell;
				nRemote[f]++;
			}
		}
		// The owners send their cell back
		vector<BigNumber> ownerCell(nOwned, -1), ghostOwnerCell(nGhosts, -1);
		for (BigNumber f = 0; f < nOwned; f++)
		{
			if (nRemote[f] == 0)
				continue;
			if (nRemote[f] > 1 || faceCells.Size(f) != 1)
				invalid = true; // more than two cells, or a face without cell on its process
			else
				ownerCell[f] = cellBegin + faceCells[f][0];
		}
		faces.UpdateGhosts(ownerCell.data(), ghostOwnerCell.data(), 1);
		for (BigNumber g = 0; g < nGhosts; g++)
		{
			if (ghostCell[g] < 0)
				continue; // not a face of the cells of the process
			if (ghostOwnerCell[g] < 0)
				invalid = true; // a face owned by a process that has none of its cells
			faceRemoteGlobal[nOwned + g] = ghostOwnerCell[g];
		}
		if (Comm().Any(invalid))
			Utils::FatalError("UncondensedAMG: on several processes, each face must belong to a process that owns one of its cells, and have at most two cells.");

		vector<BigNumber> remoteCells;
		for (BigNumber c : faceRemoteGlobal)
			if (c >= 0)
				remoteCells.push_back(c);
		CellMap = make_shared<const IndexMap>(CellOwners, std::move(remoteCells));
		FaceRemoteCell.assign(faces.NLocal(), -1);
		for (BigNumber f = 0; f < faces.NLocal(); f++)
			if (faceRemoteGlobal[f] >= 0)
				FaceRemoteCell[f] = CellMap->Local(faceRemoteGlobal[f]);
	}

	// After Pass.Aggregate(): the aggregates of the remote cells (GhostAggregates), and the remote cells added to the
	// cells of the faces (Pass.FaceCells: after the cells of the process) and to the cells of the aggregates
	// (Pass.CellAggregate: NAggregates + i for GhostAggregates[i])
	void AddRemoteAggregates()
	{
		CoarseningPass& pass = Pass;
		const IndexMap& cells = *CellMap;
		AggregateOwners = IndexMap(Comm(), pass.NAggregates);
		vector<BigNumber> cellAggregate(pass.NCells), ghostCellAggregate(cells.NGhosts());
		for (BigNumber e = 0; e < pass.NCells; e++)
			cellAggregate[e] = AggregateOwners.Begin + pass.CellAggregate[e];
		cells.UpdateGhosts(cellAggregate.data(), ghostCellAggregate.data(), 1);
		GhostAggregates = ghostCellAggregate;
		sort(GhostAggregates.begin(), GhostAggregates.end());
		GhostAggregates.erase(unique(GhostAggregates.begin(), GhostAggregates.end()), GhostAggregates.end());

		pass.CellAggregate.resize(pass.NCells + cells.NGhosts());
		for (BigNumber g = 0; g < cells.NGhosts(); g++)
			pass.CellAggregate[pass.NCells + g] = pass.NAggregates + (lower_bound(GhostAggregates.begin(), GhostAggregates.end(), ghostCellAggregate[g]) - GhostAggregates.begin());

		IndexLists faceCells;
		faceCells.Allocate(pass.NFaces, [&](BigNumber f) { return pass.FaceCells.Size(f) + (FaceRemoteCell[f] >= 0 ? 1 : 0); });
		for (BigNumber f = 0; f < pass.NFaces; f++)
		{
			copy(pass.FaceCells[f], pass.FaceCells[f] + pass.FaceCells.Size(f), faceCells[f]);
			if (FaceRemoteCell[f] >= 0)
				faceCells[f][pass.FaceCells.Size(f)] = FaceRemoteCell[f]; // after the cells of the process: still increasing
		}
		pass.FaceCells = std::move(faceCells);
	}

	// After Pass.CollapseInterfaces(): the global numbers of the coarse faces. A process numbers its coarse faces (those
	// of the boundaries of its aggregates included), in their order in the pass, and asks for the others (of the
	// processes of smaller rank) to their owners. CoarseFaceMap:
	// the coarse faces of the aggregates of the process; CoarseFacePosition: their local numbers.
	void NumberCoarseFaces()
	{
		CoarseningPass& pass = Pass;
		const Communicator& comm = Comm();
		BigNumber nAggregates = pass.NAggregates;
		auto remoteAggregate = [&](BigNumber a) { return GhostAggregates[a - nAggregates]; };
		vector<char> owned(pass.NCoarseFaces);
		BigNumber nOwned = 0;
		for (BigNumber cf = 0; cf < pass.NCoarseFaces; cf++)
		{
			BigNumber b = pass.CoarseFaceAggregates[cf].second;
			owned[cf] = b < nAggregates || AggregateOwners.Owner(remoteAggregate(b)) > comm.Rank();
			nOwned += owned[cf];
		}
		IndexMap coarseFaceOwners(comm, nOwned);
		vector<BigNumber> global(pass.NCoarseFaces, -1);
		BigNumber next = coarseFaceOwners.Begin;
		for (BigNumber cf = 0; cf < pass.NCoarseFaces; cf++)
			if (owned[cf])
				global[cf] = next++;

		// Requests for the others: (aggregate of the owner, aggregate of this process)
		vector<vector<BigNumber>> requests(comm.Size());
		vector<vector<BigNumber>> requested(comm.Size()); // the coarse faces of each request
		for (BigNumber cf = 0; cf < pass.NCoarseFaces; cf++)
		{
			if (owned[cf])
				continue;
			auto [a, b] = pass.CoarseFaceAggregates[cf];
			int owner = AggregateOwners.Owner(remoteAggregate(b));
			requests[owner].push_back(remoteAggregate(b));
			requests[owner].push_back(AggregateOwners.Begin + a);
			requested[owner].push_back(cf);
		}
		vector<vector<BigNumber>> received = comm.AllToAll(requests);
		vector<vector<BigNumber>> replies(comm.Size());
		bool notFound = false;
		for (int r = 0; r < comm.Size(); r++)
		{
			for (size_t q = 0; q < received[r].size(); q += 2)
			{
				BigNumber a = received[r][q] - AggregateOwners.Begin; // aggregate of this process
				BigNumber b = nAggregates + (lower_bound(GhostAggregates.begin(), GhostAggregates.end(), received[r][q + 1]) - GhostAggregates.begin());
				// The coarse faces of a, by increasing second aggregate
				const BigNumber* faces = pass.AggregateCoarseFaces[a];
				int n = pass.AggregateCoarseFaces.Size(a);
				BigNumber cf = -1;
				for (int i = 0; i < n; i++)
				{
					if (pass.CoarseFaceAggregates[faces[i]] == pair<BigNumber, BigNumber>(a, b))
					{
						cf = faces[i];
						break;
					}
				}
				if (cf < 0 || !owned[cf])
				{
					notFound = true;
					replies[r].push_back(-1);
				}
				else
					replies[r].push_back(global[cf]);
			}
		}
		vector<vector<BigNumber>> answers = comm.AllToAll(replies);
		if (comm.Any(notFound))
			Utils::FatalError("UncondensedAMG: a coarse face between two processes is unknown to its owner.");
		vector<BigNumber> ghosts;
		for (int r = 0; r < comm.Size(); r++)
		{
			for (size_t q = 0; q < requested[r].size(); q++)
			{
				global[requested[r][q]] = answers[r][q];
				ghosts.push_back(answers[r][q]);
			}
		}
		CoarseFaceMap = make_shared<const IndexMap>(coarseFaceOwners, std::move(ghosts));
		CoarseFacePosition.resize(pass.NCoarseFaces);
		for (BigNumber cf = 0; cf < pass.NCoarseFaces; cf++)
			CoarseFacePosition[cf] = CoarseFaceMap->Local(global[cf]);
	}
};
