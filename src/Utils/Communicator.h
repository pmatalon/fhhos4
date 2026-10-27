#pragma once
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <type_traits>
#include <vector>
#include <unistd.h>
#ifdef MPI_ENABLED
#include <mpi.h>
#endif
#include "Parallelism.h"
#include "Types.h"
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
// processes of Communicator::Current(). On one process, those of Eigen, bit for bit (x.norm() is sqrt(x.squaredNorm())).
namespace Parallel
{
	inline double Dot(const Vector& x, const Vector& y)
	{
		return Communicator::Current().Sum(x.dot(y));
	}

	inline double SquaredNorm(const Vector& x)
	{
		return Communicator::Current().Sum(x.squaredNorm());
	}

	inline double Norm(const Vector& x)
	{
		const Communicator& comm = Communicator::Current();
		if (!comm.IsDistributed())
			return x.norm();
		return std::sqrt(comm.Sum(x.squaredNorm()));
	}

	// The same for vectors in single precision (the cycle of U-AMG in single precision): Eigen's products of blocks of
	// 1024 coefficients, summed in double
	inline double Dot(const FloatVector& x, const FloatVector& y)
	{
		double dot = 0;
		for (BigNumber i = 0; i < x.rows(); i += 1024)
		{
			BigNumber m = std::min<BigNumber>(1024, x.rows() - i);
			dot += x.segment(i, m).dot(y.segment(i, m));
		}
		return Communicator::Current().Sum(dot);
	}

	inline double SquaredNorm(const FloatVector& x)
	{
		return Dot(x, x);
	}

	inline double Norm(const FloatVector& x)
	{
		return std::sqrt(SquaredNorm(x));
	}

	// The global number of rows of a vector distributed by rows
	inline BigNumber Rows(const Vector& x)
	{
		return Communicator::Current().Sum((BigNumber)x.rows());
	}
}
