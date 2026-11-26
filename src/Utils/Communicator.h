#pragma once
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <string>
#include <type_traits>
#include <vector>
#include <unistd.h>
#ifdef MPI_ENABLED
#include <mpi.h>
#endif
#include "Parallelism.h"
#include "Types.h"
#include "fhhos4_error.h"
using namespace std;

// The processes of a distributed computation (MPI): their number, and the reductions and exchanges between them.
// Without MPI (ENABLE_MPI=OFF), a single process: the reductions return their argument. The code calls these functions
// unconditionally: no #ifdef MPI_ENABLED outside this file.
//
// The distributed objects (IndexMap, DistributedMatrix) are only built on several processes: on one process, the
// sequential code runs, unchanged (the same results and performance).
class Communicator
{
private:
#ifdef MPI_ENABLED
	MPI_Comm _comm = MPI_COMM_SELF;
#endif
	int _rank = 0;
	int _size = 1;

public:
	// A single process (MPI_COMM_SELF with MPI)
	Communicator() {}

#ifdef MPI_ENABLED
	explicit Communicator(MPI_Comm comm)
		: _comm(comm)
	{
		MPI_Comm_rank(comm, &_rank);
		MPI_Comm_size(comm, &_size);
	}

	MPI_Comm Comm() const { return _comm; }
#endif

	// The Fortran handle of the communicator, for MUMPS; without MPI, MUMPS's USE_COMM_WORLD (its sequential version)
	int FortranHandle() const
	{
#ifdef MPI_ENABLED
		return MPI_Comm_c2f(_comm);
#else
		return -987654;
#endif
	}

	// Whether MPI is finalized (no MPI call after it); false without MPI
	static bool Finalized()
	{
		int finalized = 0;
#ifdef MPI_ENABLED
		MPI_Finalized(&finalized);
#endif
		return finalized != 0;
	}

	int Rank() const { return _rank; }
	int Size() const { return _size; }
	bool IsDistributed() const { return _size > 1; }
	bool IsRoot() const { return _rank == 0; }

	// The communicator of the reductions of the solvers (Parallel::Dot(), Norm()...), shared by all of them: the
	// program's MPI_COMM_WORLD, or, during a call of the library, the communicator of its fhhos4::Solver (see
	// CurrentCommunicator). A single process by default.
	static const Communicator*& CurrentPointer()
	{
		static const Communicator self;
		static const Communicator* current = &self;
		return current;
	}

	static const Communicator& Current()
	{
		return *CurrentPointer();
	}

	//-------------------------------------------------------------------//
	//                            Reductions                             //
	//-------------------------------------------------------------------//

	template <class T>
	T Sum(T x) const
	{
#ifdef MPI_ENABLED
		if (_size > 1)
			MPI_Allreduce(MPI_IN_PLACE, &x, 1, Type<T>(), MPI_SUM, _comm);
#endif
		return x;
	}

	template <class T>
	T Max(T x) const
	{
#ifdef MPI_ENABLED
		if (_size > 1)
			MPI_Allreduce(MPI_IN_PLACE, &x, 1, Type<T>(), MPI_MAX, _comm);
#endif
		return x;
	}

	template <class T>
	T Min(T x) const
	{
#ifdef MPI_ENABLED
		if (_size > 1)
			MPI_Allreduce(MPI_IN_PLACE, &x, 1, Type<T>(), MPI_MIN, _comm);
#endif
		return x;
	}

	// Whether x is true on at least one process (resp. on all of them)
	bool Any(bool x) const { return Max((int)x) != 0; }
	bool All(bool x) const { return Min((int)x) != 0; }

	// The sums of values[0..n) over the processes, in place
	template <class T>
	void SumInPlace(T* values, size_t n) const
	{
#ifdef MPI_ENABLED
		if (_size > 1 && n > 0)
			MPI_Allreduce(MPI_IN_PLACE, values, (int)n, Type<T>(), MPI_SUM, _comm);
#endif
	}

	// The sum of x over the processes of lower rank (0 on the first)
	template <class T>
	T ExclusiveScan(T x) const
	{
		T result = 0;
#ifdef MPI_ENABLED
		if (_size > 1)
		{
			MPI_Exscan(&x, &result, 1, Type<T>(), MPI_SUM, _comm);
			if (_rank == 0)
				result = 0; // undefined on rank 0
		}
#endif
		return result;
	}

	// The value x of each process, in rank order
	template <class T>
	vector<T> AllGather(T x) const
	{
		vector<T> values(_size, x);
#ifdef MPI_ENABLED
		if (_size > 1)
			MPI_Allgather(&x, 1, Type<T>(), values.data(), 1, Type<T>(), _comm);
#endif
		return values;
	}

	// The concatenation of the vectors of all the processes, in rank order (counts: the size of each, if not null)
	template <class T>
	vector<T> AllGatherV(const vector<T>& local, vector<int>* counts = nullptr) const
	{
		static_assert(is_trivially_copyable_v<T>);
		vector<int> sizes = AllGather((int)(local.size() * sizeof(T)));
		vector<int> displacements(_size + 1, 0);
		for (int r = 0; r < _size; r++)
			displacements[r + 1] = displacements[r] + sizes[r];
		vector<T> all(displacements[_size] / sizeof(T));
#ifdef MPI_ENABLED
		if (_size > 1)
			MPI_Allgatherv(local.data(), sizes[_rank], MPI_BYTE, all.data(), sizes.data(), displacements.data(), MPI_BYTE, _comm);
		else
#endif
			copy(local.begin(), local.end(), all.begin());
		if (counts)
		{
			counts->resize(_size);
			for (int r = 0; r < _size; r++)
				(*counts)[r] = sizes[r] / sizeof(T);
		}
		return all;
	}

	// The concatenation of the vectors of all the processes, in rank order, on the process root (empty on the others);
	// counts: the size of each, on root, if not null
	template <class T>
	vector<T> GatherV(const vector<T>& local, int root = 0, vector<int>* counts = nullptr) const
	{
		static_assert(is_trivially_copyable_v<T>);
#ifdef MPI_ENABLED
		if (_size > 1)
		{
			int bytes = Bytes<T>(local.size());
			vector<int> sizes(_rank == root ? _size : 0);
			MPI_Gather(&bytes, 1, MPI_INT, sizes.data(), 1, MPI_INT, root, _comm);
			vector<int> displacements(sizes.size() + 1, 0);
			for (size_t r = 0; r < sizes.size(); r++)
				displacements[r + 1] = displacements[r] + sizes[r];
			vector<T> all(displacements.back() / sizeof(T));
			MPI_Gatherv(local.data(), bytes, MPI_BYTE, all.data(), sizes.data(), displacements.data(), MPI_BYTE, root, _comm);
			if (counts)
			{
				counts->resize(sizes.size());
				for (size_t r = 0; r < sizes.size(); r++)
					(*counts)[r] = sizes[r] / sizeof(T);
			}
			return all;
		}
#endif
		if (counts)
			*counts = { (int)local.size() };
		return local;
	}

	// The other way around: counts[r] values of all, in rank order, to the process r, from the process root (all and
	// counts read there only); returns the count values of this process
	template <class T>
	vector<T> ScatterV(const vector<T>& all, const vector<int>& counts, size_t count, int root = 0) const
	{
		static_assert(is_trivially_copyable_v<T>);
		vector<T> local(count);
#ifdef MPI_ENABLED
		if (_size > 1)
		{
			vector<int> sizes, displacements;
			if (_rank == root)
			{
				sizes.resize(_size);
				displacements.assign(_size, 0);
				for (int r = 0; r < _size; r++)
				{
					sizes[r] = Bytes<T>(counts[r]);
					if (r > 0)
						displacements[r] = displacements[r - 1] + sizes[r - 1];
				}
			}
			MPI_Scatterv(all.data(), sizes.data(), displacements.data(), MPI_BYTE, local.data(), Bytes<T>(count), MPI_BYTE, root, _comm);
			return local;
		}
#endif
		copy(all.begin(), all.begin() + count, local.begin());
		return local;
	}

	// The number of processes of this communicator on the node of this one that run on the same CPUs (the same
	// affinity), among which its cores are divided (see Parallelism::SharingProcesses()). Collective.
	int SharingProcesses() const
	{
		int n = 1;
#if defined(MPI_ENABLED) && defined(__linux__)
		if (_size == 1)
			return 1;
		MPI_Comm node;
		MPI_Comm_split_type(_comm, MPI_COMM_TYPE_SHARED, _rank, MPI_INFO_NULL, &node);
		int nodeSize;
		MPI_Comm_size(node, &nodeSize);
		cpu_set_t cpus;
		CPU_ZERO(&cpus);
		sched_getaffinity(0, sizeof(cpus), &cpus);
		vector<cpu_set_t> all(nodeSize);
		MPI_Allgather(&cpus, sizeof(cpus), MPI_BYTE, all.data(), sizeof(cpus), MPI_BYTE, node);
		n = 0;
		for (const cpu_set_t& other : all)
			n += CPU_EQUAL(&other, &cpus) ? 1 : 0;
		MPI_Comm_free(&node);
#endif
		return n;
	}

	void Barrier() const
	{
#ifdef MPI_ENABLED
		if (_size > 1)
			MPI_Barrier(_comm);
#endif
	}

	// Sends the bytes of the root to all the processes (resized on the others)
	template <class T>
	void Broadcast(vector<T>& values, int root = 0) const
	{
		static_assert(is_trivially_copyable_v<T>);
#ifdef MPI_ENABLED
		if (_size == 1)
			return;
		size_t n = values.size();
		MPI_Bcast(&n, sizeof(n), MPI_BYTE, root, _comm);
		values.resize(n);
		if (n > 0)
			MPI_Bcast(values.data(), (int)(n * sizeof(T)), MPI_BYTE, root, _comm);
#endif
	}

	//-------------------------------------------------------------------//
	//                            Exchanges                              //
	//-------------------------------------------------------------------//

	// Personalized exchange: send[r] goes to the process r; returns what each process sent to this one (received[r]
	// from the process r). Collective: every process calls it, with send of size Size(). The sizes are exchanged with
	// MPI_Alltoall (O(Size()) per process): for the setups, not for the solves (see Neighbours below).
	template <class T>
	vector<vector<T>> AllToAll(const vector<vector<T>>& send) const
	{
		static_assert(is_trivially_copyable_v<T>);
		assert((int)send.size() == _size);
		vector<vector<T>> received(_size);
#ifdef MPI_ENABLED
		if (_size > 1)
		{
			vector<int64_t> sendCounts(_size), recvCounts(_size);
			for (int r = 0; r < _size; r++)
				sendCounts[r] = send[r].size();
			MPI_Alltoall(sendCounts.data(), 1, MPI_INT64_T, recvCounts.data(), 1, MPI_INT64_T, _comm);
			vector<MPI_Request> requests;
			requests.reserve(2 * _size);
			for (int r = 0; r < _size; r++)
			{
				received[r].resize(recvCounts[r]);
				if (recvCounts[r] > 0)
				{
					requests.emplace_back();
					MPI_Irecv(received[r].data(), Bytes<T>(recvCounts[r]), MPI_BYTE, r, ExchangeTag, _comm, &requests.back());
				}
			}
			for (int r = 0; r < _size; r++)
			{
				if (sendCounts[r] > 0)
				{
					requests.emplace_back();
					MPI_Isend(send[r].data(), Bytes<T>(sendCounts[r]), MPI_BYTE, r, ExchangeTag, _comm, &requests.back());
				}
			}
			MPI_Waitall((int)requests.size(), requests.data(), MPI_STATUSES_IGNORE);
			return received;
		}
#endif
		received[0] = send[0];
		return received;
	}

	// Exchange with known neighbours: sends sendBuffers[k] to sendRanks[k] and receives recvBuffers[k] (already sized)
	// from recvRanks[k]. Only the processes that exchange take part: for the solves (halo updates).
	template <class T>
	void Neighbours(const vector<int>& sendRanks, const vector<vector<T>>& sendBuffers, const vector<int>& recvRanks, vector<vector<T>>& recvBuffers, int tag = NeighbourTag) const
	{
		static_assert(is_trivially_copyable_v<T>);
#ifdef MPI_ENABLED
		vector<MPI_Request> requests;
		requests.reserve(sendRanks.size() + recvRanks.size());
		for (size_t k = 0; k < recvRanks.size(); k++)
		{
			if (recvBuffers[k].empty())
				continue;
			requests.emplace_back();
			MPI_Irecv(recvBuffers[k].data(), Bytes<T>(recvBuffers[k].size()), MPI_BYTE, recvRanks[k], tag, _comm, &requests.back());
		}
		for (size_t k = 0; k < sendRanks.size(); k++)
		{
			if (sendBuffers[k].empty())
				continue;
			requests.emplace_back();
			MPI_Isend(sendBuffers[k].data(), Bytes<T>(sendBuffers[k].size()), MPI_BYTE, sendRanks[k], tag, _comm, &requests.back());
		}
		MPI_Waitall((int)requests.size(), requests.data(), MPI_STATUSES_IGNORE);
#else
		assert(sendRanks.empty() && recvRanks.empty());
#endif
	}

private:
	static constexpr int ExchangeTag = 7101;
	static constexpr int NeighbourTag = 7102;

	template <class T>
	static int Bytes(size_t n)
	{
		size_t bytes = n * sizeof(T);
		assert(bytes <= (size_t)INT32_MAX); // MPI's counts are int
		return (int)bytes;
	}

#ifdef MPI_ENABLED
	template <class T>
	static MPI_Datatype Type()
	{
		if constexpr (is_same_v<T, double>)
			return MPI_DOUBLE;
		else if constexpr (is_same_v<T, float>)
			return MPI_FLOAT;
		else if constexpr (is_integral_v<T> && is_signed_v<T> && sizeof(T) == 4)
			return MPI_INT32_T;
		else if constexpr (is_integral_v<T> && is_signed_v<T> && sizeof(T) == 8)
			return MPI_INT64_T;
		else if constexpr (is_integral_v<T> && is_unsigned_v<T> && sizeof(T) == 8)
			return MPI_UINT64_T;
		else
			static_assert(sizeof(T) == 0, "Communicator: unsupported type");
	}
#endif
};

// Runs check() on all the processes of comm: if it throws on one of them, all of them throw (an error on one process
// alone would leave the others waiting in the next collective call), with the message of the first process where it
// threw
template <class Check>
void CollectiveCheck(const Communicator& comm, Check check)
{
	if (!comm.IsDistributed())
	{
		check();
		return;
	}
	string message;
	try
	{
		check();
	}
	catch (const fhhos4::Error& e)
	{
		message = e.what();
	}
	int failedRank = comm.Min(message.empty() ? comm.Size() : comm.Rank());
	if (failedRank == comm.Size())
		return;
	vector<char> text(message.begin(), message.end());
	comm.Broadcast(text, failedRank);
	throw fhhos4::Error(string(text.begin(), text.end()) + " (on the process " + to_string(failedRank) + ")");
}

// The MPI environment of an executable (the program, the tests): MPI_Init_thread() at the construction, MPI_Finalize()
// at the destruction or at exit() (e.g. after the help, or an error in the arguments), and the communicator of all the
// processes. Without MPI, a single process. OpenMP in each process: MPI is called by the main thread only
// (MPI_THREAD_FUNNELED).
class MPIEnvironment
{
private:
	Communicator _world;
public:
	MPIEnvironment(int& argc, char**& argv)
	{
#ifdef MPI_ENABLED
		int provided;
		MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
		InitializingProcess() = getpid();
		atexit(FinalizeIfNeeded);
		_world = Communicator(MPI_COMM_WORLD);
		// The processes that run on the same CPUs share them: by default (without OMP_NUM_THREADS), as many threads per
		// process as CPUs of the affinity divided among them, in the solvers as many as physical cores divided among them
		int sharing = _world.SharingProcesses();
		Parallelism::SharingProcesses() = sharing;
#ifdef _OPENMP
		if (sharing > 1 && !getenv("OMP_NUM_THREADS"))
			omp_set_num_threads(max(1, omp_get_max_threads() / sharing));
#endif
#endif
	}

	const Communicator& World() const { return _world; }

	// Ends all the processes (an error on some of them only would leave the others waiting)
	[[noreturn]] static void Abort(int exitCode)
	{
#ifdef MPI_ENABLED
		MPI_Abort(MPI_COMM_WORLD, exitCode);
#endif
		exit(exitCode);
	}

	~MPIEnvironment()
	{
		FinalizeIfNeeded();
	}

	MPIEnvironment(const MPIEnvironment&) = delete;
	MPIEnvironment& operator=(const MPIEnvironment&) = delete;

private:
	// The process that initialized MPI: not a child forked from it (e.g. by the death tests of GoogleTest), which must
	// not finalize MPI (MPI_Finalize() would wait for the other processes)
	static pid_t& InitializingProcess()
	{
		static pid_t pid = 0;
		return pid;
	}

	static void FinalizeIfNeeded()
	{
#ifdef MPI_ENABLED
		if (getpid() != InitializingProcess())
			return;
		int initialized = 0, finalized = 0;
		MPI_Initialized(&initialized);
		MPI_Finalized(&finalized);
		if (initialized && !finalized)
			MPI_Finalize();
#endif
	}
};

// While in scope: the communicator of the reductions of the solvers (Communicator::Current())
class CurrentCommunicator
{
private:
	const Communicator* _saved;
public:
	explicit CurrentCommunicator(const Communicator& comm)
	{
		_saved = Communicator::CurrentPointer();
		Communicator::CurrentPointer() = &comm;
	}
	~CurrentCommunicator()
	{
		Communicator::CurrentPointer() = _saved;
	}
	CurrentCommunicator(const CurrentCommunicator&) = delete;
	CurrentCommunicator& operator=(const CurrentCommunicator&) = delete;
};

// The reductions of the solvers on vectors distributed by rows (each process holds the rows it owns), over the
// processes of Communicator::Current(). On a process, by blocks of rows in parallel (RowBlocks): the results depend on
// the number of rows, not on the number of threads.
namespace Parallel
{
	// The blocks of rows of the vector operations of the solvers: each block computed by one thread (from 20000 rows, the
	// threshold of the other vector loops), its sums added to the others' in the order of the blocks. 4096 rows: the rows
	// of a few vectors stay in the cache of a core, so a pass can update a block then reduce it (SumsOverRowBlocks()); more
	// above 4M rows (at most MaxCount blocks, whose sums are on the stack). One block (at most 4096 rows): the sum of the
	// whole rows, e.g. Eigen's x.dot(y), bit for bit.
	struct RowBlocks
	{
		static constexpr BigNumber MaxCount = 1024;
		BigNumber Size;  // the rows of a block (the last one: the rest)
		BigNumber Count;

		explicit RowBlocks(BigNumber n)
			: Size(4096 * std::max<BigNumber>(1, (n + 4096 * MaxCount - 1) / (4096 * MaxCount))), Count((n + Size - 1) / Size)
		{}
	};

	// f(first, size) on the blocks of rows of a vector of n rows (RowBlocks), in parallel
	template <class F>
	void ForEachRowBlock(BigNumber n, F f)
	{
		RowBlocks blocks(n);
		#pragma omp parallel for schedule(static) if (n > 20000)
		for (BigNumber b = 0; b < blocks.Count; ++b)
			f(b * blocks.Size, std::min(blocks.Size, n - b * blocks.Size));
	}

	// The K sums of blockSums(first, size) (an std::array<double, K>: those of the rows [first, first + size), which it may
	// also update) over the blocks of rows of a vector of n rows, in the order of the blocks (RowBlocks); on this process
	template <int K, class BlockSums>
	std::array<double, K> SumsOverRowBlocks(BigNumber n, BlockSums blockSums)
	{
		RowBlocks blocks(n);
		std::array<double, K> partial[RowBlocks::MaxCount];
		ForEachRowBlock(n, [&](BigNumber first, BigNumber size) { partial[first / blocks.Size] = blockSums(first, size); });
		std::array<double, K> sums{};
		for (BigNumber b = 0; b < blocks.Count; ++b)
			for (int k = 0; k < K; ++k)
				sums[k] += partial[b][k];
		return sums;
	}

	// The parts of Dot() and SquaredNorm() on a block of rows [first, first + size): Eigen's; in single precision (the
	// cycle of U-AMG in single precision), Eigen's products of sub-blocks of 1024 coefficients, summed in double
	inline double DotOfRows(const Vector& x, const Vector& y, BigNumber first, BigNumber size)
	{
		return x.segment(first, size).dot(y.segment(first, size));
	}

	inline double DotOfRows(const FloatVector& x, const FloatVector& y, BigNumber first, BigNumber size)
	{
		double dot = 0;
		for (BigNumber i = first; i < first + size; i += 1024)
		{
			BigNumber m = std::min<BigNumber>(1024, first + size - i);
			dot += x.segment(i, m).dot(y.segment(i, m));
		}
		return dot;
	}

	inline double SquaredNormOfRows(const Vector& x, BigNumber first, BigNumber size)
	{
		return x.segment(first, size).squaredNorm();
	}

	inline double SquaredNormOfRows(const FloatVector& x, BigNumber first, BigNumber size)
	{
		return DotOfRows(x, x, first, size);
	}

	template <class V>
	double DotOf(const V& x, const V& y)
	{
		return Communicator::Current().Sum(SumsOverRowBlocks<1>(x.rows(), [&](BigNumber i, BigNumber m)
		{
			return std::array<double, 1>{ DotOfRows(x, y, i, m) };
		})[0]);
	}

	template <class V>
	double SquaredNormOf(const V& x)
	{
		return Communicator::Current().Sum(SumsOverRowBlocks<1>(x.rows(), [&](BigNumber i, BigNumber m)
		{
			return std::array<double, 1>{ SquaredNormOfRows(x, i, m) };
		})[0]);
	}

	inline double Dot(const Vector& x, const Vector& y)       { return DotOf(x, y); }
	inline double Dot(const FloatVector& x, const FloatVector& y) { return DotOf(x, y); }
	inline double SquaredNorm(const Vector& x)                { return SquaredNormOf(x); }
	inline double SquaredNorm(const FloatVector& x)           { return SquaredNormOf(x); }
	inline double Norm(const Vector& x)                       { return std::sqrt(SquaredNormOf(x)); }
	inline double Norm(const FloatVector& x)                  { return std::sqrt(SquaredNormOf(x)); }

	// The global number of rows of a vector distributed by rows
	inline GlobalNumber Rows(const Vector& x)
	{
		return Communicator::Current().Sum((GlobalNumber)x.rows());
	}
}
