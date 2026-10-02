#pragma once
#include <cassert>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
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
	// Sets the number of threads of the parallel loops (and of Eigen's internal parallelism).
	// 0: back to the OpenMP default (OMP_NUM_THREADS if set, the number of cores otherwise).
	inline void SetNThreads(int nThreads)
	{
#ifdef _OPENMP
		static const int defaultNThreads = omp_get_max_threads(); // initialized at the first call, before any change
		omp_set_num_threads(nThreads > 0 ? nThreads : defaultNThreads);
#endif
	}

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
