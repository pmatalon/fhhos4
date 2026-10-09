#pragma once
#include <cassert>
#include <cstdlib>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef __linux__
#include <sched.h>
#include <fstream>
#include <set>
#include <string>
#endif
#include "NonZeroCoefficients.h"
using namespace std;

// The parallel loops are OpenMP loops:
//
//     #pragma omp parallel for
//     for (Element<Dim>* e : mesh->Elements)
//         ...
//
// Without OpenMP (ENABLE_OPENMP=OFF), the pragmas are ignored and the loops are executed sequentially.
// A parallel loop nested in another one is executed sequentially by the calling thread (OpenMP default).
namespace Parallelism
{
	// Number of physical cores among the CPUs the calling thread may run on (its affinity: e.g. those of the job, of the
	// MPI process), each core counted once whatever its number of hardware threads (hyper-threading): the number of
	// distinct sets of CPUs sharing a core (thread_siblings_list, in /sys on Linux, whatever the numbering of the cores
	// and packages). 0 if unknown: on other systems, if /sys can't be read, or beyond 1024 CPUs (CPU_SETSIZE).
	inline int PhysicalCores()
	{
#ifdef __linux__
		static const int cores = []()
		{
			cpu_set_t cpus;
			if (sched_getaffinity(0, sizeof(cpus), &cpus) != 0)
				return 0;
			set<string> siblingSets;
			for (int cpu = 0; cpu < CPU_SETSIZE; cpu++)
			{
				if (!CPU_ISSET(cpu, &cpus))
					continue;
				ifstream file("/sys/devices/system/cpu/cpu" + to_string(cpu) + "/topology/thread_siblings_list");
				string siblings;
				if (!getline(file, siblings) || siblings.empty())
					return 0;
				siblingSets.insert(siblings);
			}
			return (int)siblingSets.size();
		}();
		return cores;
#else
		return 0;
#endif
	}

	// The automatic number of threads of the solvers, from the OpenMP setting openMPThreads (omp_get_max_threads():
	// OMP_NUM_THREADS if set, otherwise the number of CPUs of the affinity): at most one thread per physical core,
	// without hyper-threading, unless the user set OpenMP's number of threads (OMP_NUM_THREADS) or its binding
	// (OMP_PROC_BIND, OMP_PLACES: the calling thread is then pinned to its place, whose CPUs aren't the process's).
	// The hyperthreads slow the solvers down (memory-bound), whereas they speed up the assembly (compute-bound):
	// PERFORMANCE.md.
	inline int WithoutHyperThreads(int openMPThreads)
	{
#ifdef _OPENMP
		if (getenv("OMP_NUM_THREADS") || omp_get_proc_bind() != omp_proc_bind_false)
			return openMPThreads;
#endif
		int cores = PhysicalCores();
		return cores > 0 ? min(openMPThreads, cores) : openMPThreads;
	}

	// The number of threads given to SetNThreads() (0: automatic)
	inline int& RequestedNThreads()
	{
		static int nThreads = 0;
		return nThreads;
	}

	// Sets the number of threads of the parallel loops (and of Eigen's internal parallelism).
	// 0: back to the automatic default, the OpenMP default (OMP_NUM_THREADS if set, all the logical CPUs otherwise), and
	// in the solvers (SolverThreads) at most one thread per physical core.
	inline void SetNThreads(int nThreads)
	{
		RequestedNThreads() = nThreads;
#ifdef _OPENMP
		static const int defaultNThreads = omp_get_max_threads(); // initialized at the first call, before any change
		omp_set_num_threads(nThreads > 0 ? nThreads : defaultNThreads);
#endif
	}

	// While in scope (or until End()): the number of threads of the solvers, WithoutHyperThreads() unless a number of
	// threads was given to SetNThreads() (-threads), which then applies to the solvers too
	class SolverThreads
	{
	private:
		int _savedNThreads = 0;
	public:
		SolverThreads()
		{
#ifdef _OPENMP
			int nThreads = omp_get_max_threads();
			int solverNThreads = RequestedNThreads() > 0 ? nThreads : WithoutHyperThreads(nThreads);
			if (solverNThreads != nThreads)
			{
				_savedNThreads = nThreads;
				omp_set_num_threads(solverNThreads);
			}
#endif
		}

		// Back to the number of threads before
		void End()
		{
#ifdef _OPENMP
			if (_savedNThreads > 0)
				omp_set_num_threads(_savedNThreads);
#endif
			_savedNThreads = 0;
		}

		~SolverThreads()
		{
			End();
		}

		SolverThreads(const SolverThreads&) = delete;
		SolverThreads& operator=(const SolverThreads&) = delete;
	};

	// Number of threads that would execute a parallel loop started by the calling thread
	inline int NThreads()
	{
#ifdef _OPENMP
		if (omp_get_active_level() >= omp_get_max_active_levels()) // nested: executed by the calling thread only
			return 1;
		return omp_get_max_threads();
#else
		return 1;
#endif
	}

	// Number of the calling thread in the current parallel loop (0 outside a parallel loop)
	inline int ThreadNumber()
	{
#ifdef _OPENMP
		return omp_get_thread_num();
#else
		return 0;
#endif
	}

	// Iterations [begin, end) of the calling thread when n iterations are split into contiguous chunks, one per thread of
	// the current parallel region, in thread order ([0, n) outside a parallel region)
	inline pair<BigNumber, BigNumber> ThreadChunk(BigNumber n)
	{
#ifdef _OPENMP
		BigNumber thread = omp_get_thread_num(), nThreads = omp_get_num_threads();
		return { n * thread / nThreads, n * (thread + 1) / nThreads };
#else
		return { 0, n };
#endif
	}

	// Max number of iterations processed by one thread in a parallel loop of n iterations (static schedule)
	inline BigNumber ChunkSize(BigNumber n)
	{
		BigNumber nThreads = NThreads();
		return (n + nThreads - 1) / nThreads;
	}
}

// One value of T per thread, to accumulate results in a parallel loop without synchronization:
//
//     ThreadLocal<vector<Element<Dim>*>> selected;
//     #pragma omp parallel for schedule(static)
//     for (Element<Dim>* e : mesh->Elements)
//         if (...)
//             selected.Local().push_back(e);
//     for (vector<Element<Dim>*>& list : selected) // in thread order
//         ...
//
// With schedule(static), each thread processes one contiguous chunk of iterations, and the chunks follow the
// thread order: the results, taken in thread order, follow the iteration order.
// Must be constructed outside the parallel loop.
template <class T>
class ThreadLocal
{
private:
	vector<T> _values;
public:
	// Each value is constructed with the arguments args
	template <class... Args>
	explicit ThreadLocal(const Args&... args)
	{
		int nThreads = Parallelism::NThreads();
		_values.reserve(nThreads);
		for (int i = 0; i < nThreads; i++)
			_values.emplace_back(args...);
	}

	ThreadLocal(const ThreadLocal&) = delete;
	ThreadLocal& operator=(const ThreadLocal&) = delete;

	// Value of the calling thread
	T& Local()
	{
		int threadNumber = Parallelism::ThreadNumber();
		assert(threadNumber < (int)_values.size());
		return _values[threadNumber];
	}

	typename vector<T>::iterator begin() { return _values.begin(); }
	typename vector<T>::iterator end()   { return _values.end(); }
};

// Non-zero coefficients of a sparse matrix, added in a parallel loop:
//
//     ThreadLocalCoeffs coeffs(mesh->Elements.size(), nnzPerElement);
//     #pragma omp parallel for
//     for (Element<Dim>* e : mesh->Elements)
//         coeffs.Local().Add(i, j, value);
//     coeffs.Fill(M);
class ThreadLocalCoeffs : public ThreadLocal<NonZeroCoefficients>
{
public:
	// nIterations, nnzPerIteration: to reserve the memory of each thread
	ThreadLocalCoeffs(BigNumber nIterations = 0, BigNumber nnzPerIteration = 0) :
		ThreadLocal<NonZeroCoefficients>(Parallelism::ChunkSize(nIterations) * nnzPerIteration)
	{}

	// Concatenation of the coefficients of all threads, in thread order
	NonZeroCoefficients Merge()
	{
		BigNumber size = 0;
		for (NonZeroCoefficients& coeffs : *this)
			size += coeffs.Size();
		NonZeroCoefficients global(size);
		for (NonZeroCoefficients& coeffs : *this)
			global.Add(coeffs);
		return global;
	}

	void Fill(SparseMatrix& m)
	{
		Merge().Fill(m);
	}
};
