#pragma once
#include <cassert>
#include <thread>
#include <future>
#include <math.h>
#include <type_traits>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "NonZeroCoefficients.h"
using namespace std;

struct EmptyResultChunk
{};

class CoeffsChunk
{
public: NonZeroCoefficients Coeffs;
};

template <class ResultT = CoeffsChunk>
class ParallelChunk
{
public:
	int ThreadNumber;
	BigNumber Start; // starts at 0
	BigNumber End; // not included in the chunk
	BigNumber Size() { return End - Start; }
	std::future<void> ThreadFuture;
	ResultT Results;

	ParallelChunk(int threadNumber)
	{
		this->ThreadNumber = threadNumber;
	}
};

class BaseParallelLoop
{
protected:
	inline static unsigned int DefaultNThreads = std::thread::hardware_concurrency();

public:
	static void SetDefaultNThreads(unsigned int nThreads)
	{
		DefaultNThreads = nThreads;
		if (DefaultNThreads == 0)
		{
			DefaultNThreads = std::thread::hardware_concurrency();
			if (DefaultNThreads == 0)
			{
				cout << "Warning: std::thread::hardware_concurrency() returned 0. Falling down to sequential execution.";
				DefaultNThreads = 1;
			}
		}
#ifdef _OPENMP
		// Also used by Eigen's internal parallelism (sparse matrix-vector products, dense products)
		omp_set_num_threads(DefaultNThreads);
#endif
	}
	static unsigned int GetDefaultNThreads()
	{
		return DefaultNThreads;
	}
};

template <class ResultT>
class BaseChunksParallelLoop : public BaseParallelLoop
{
public:
	unsigned int NThreads;
	BigNumber ChunkMinSize;

	vector<ParallelChunk<ResultT>*> Chunks;

	BaseChunksParallelLoop(BigNumber loopSize, unsigned int nThreads)
	{
		NThreads = nThreads;
		if (NThreads == 0)
			NThreads = std::thread::hardware_concurrency();
		if (NThreads == 0) // hardware_concurrency() can return 0
			NThreads = 1;

		// Each thread must have at least 2 elements to process
		while (loopSize / NThreads < 2 && NThreads > 1)
			NThreads--;

		Chunks.resize(NThreads);
		ChunkMinSize = loopSize / NThreads;
		int rest = loopSize - NThreads * ChunkMinSize;

		BigNumber start = 0;
		for (unsigned int threadNumber = 0; threadNumber < NThreads; threadNumber++)
		{
			Chunks[threadNumber] = new ParallelChunk<ResultT>(threadNumber);
			Chunks[threadNumber]->Start = start;
			Chunks[threadNumber]->End = start + ChunkMinSize + (threadNumber < rest ? 1 : 0);
			start += Chunks[threadNumber]->Size();
		}
		assert(start == loopSize);
	}

	// The chunks are owned (deleted by the destructor): a copy would delete them twice.
	BaseChunksParallelLoop(const BaseChunksParallelLoop&) = delete;
	BaseChunksParallelLoop& operator=(const BaseChunksParallelLoop&) = delete;

	void InitChunks(function<void(ParallelChunk<ResultT>*)> functionInitChunks)
	{
		for (unsigned int threadNumber = 0; threadNumber < NThreads; threadNumber++)
			functionInitChunks(Chunks[threadNumber]);
	}

	void AggregateChunkResults(function<void(ResultT&)> aggregate)
	{
		for (unsigned int threadNumber = 0; threadNumber < NThreads; threadNumber++)
			aggregate(Chunks[threadNumber]->Results);
	}

	void Wait()
	{
		for (unsigned int threadNumber = 0; threadNumber < NThreads; threadNumber++)
			this->Chunks[threadNumber]->ThreadFuture.wait();
	}

	// Runs functionChunk on each chunk, one chunk per thread.
	// The chunk decomposition is static, so the results do not depend on the thread scheduling.
	template <class F>
	void ExecuteChunk(F&& functionChunk)
	{
		if (NThreads == 1)
			functionChunk(this->Chunks[0]);
		else
		{
#ifdef _OPENMP
			// The OpenMP threads are pooled: no thread creation cost at each call.
			// Nested calls (from inside a parallel loop) are executed sequentially by the calling thread.
			#pragma omp parallel for num_threads(NThreads) schedule(static, 1)
			for (int threadNumber = 0; threadNumber < (int)NThreads; threadNumber++)
				functionChunk(this->Chunks[threadNumber]);
#else
			for (unsigned int threadNumber = 0; threadNumber < NThreads; threadNumber++)
			{
				ParallelChunk<ResultT>* chunk = this->Chunks[threadNumber];
				chunk->ThreadFuture = std::async(std::launch::async, [chunk, &functionChunk]()
					{
						functionChunk(chunk);
					}
				);
			}
			this->Wait();
#endif
		}
	}

	virtual ~BaseChunksParallelLoop()
	{
		for (unsigned int threadNumber = 0; threadNumber < NThreads; threadNumber++)
		{
			ParallelChunk<ResultT>* chunk = Chunks[threadNumber];
			delete chunk;
		}
	}

	// Specialization when ResultT = CoeffsChunk
	void ReserveChunkCoeffsSize(BigNumber nnzForOneLoopIteration)
	{
		static_assert(std::is_same<ResultT, CoeffsChunk>::value, "Works only with CoeffsChunk!");
		for (unsigned int threadNumber = 0; threadNumber < this->NThreads; threadNumber++)
		{
			ParallelChunk<CoeffsChunk>* chunk = this->Chunks[threadNumber];
			chunk->Results.Coeffs = NonZeroCoefficients(chunk->Size() * nnzForOneLoopIteration);
		}
	}
	void Fill(SparseMatrix &m)
	{
		static_assert(std::is_same<ResultT, CoeffsChunk>::value, "Works only with CoeffsChunk!");
		NonZeroCoefficients global;
		for (unsigned int threadNumber = 0; threadNumber < this->NThreads; threadNumber++)
		{
			ParallelChunk<CoeffsChunk>* chunk = this->Chunks[threadNumber];
			global.Add(chunk->Results.Coeffs);
		}
		global.Fill(m);
	}
};

//----------------------------//
//     List parallel loop     //
//----------------------------//

template <class T, class ResultT = CoeffsChunk>//, typename enable_if<is_base_of<ParallelChunk, ChunkT>::value>::type = ParallelChunk >
class ParallelLoop : public BaseChunksParallelLoop<ResultT>
{
	//static_assert(std::is_base_of<ParallelChunk, ChunkT>::value, "ChunkT must inherit from ParallelChunk");
private:
	const vector<T>& _list;
public:

	ParallelLoop(const vector<T>& list) : 
		ParallelLoop(list, BaseParallelLoop::DefaultNThreads) {}
	
	ParallelLoop(const vector<T>& list, unsigned int nThreads) :
		BaseChunksParallelLoop<ResultT>(list.size(), nThreads),
		_list(list)
	{}

	// functionToExecute: void(T) or void(T, ParallelChunk<ResultT>*)
	template <class F>
	void Execute(F&& functionToExecute)
	{
		this->ExecuteChunk([this, &functionToExecute](ParallelChunk<ResultT>* chunk)
			{
				for (BigNumber i = chunk->Start; i < chunk->End; ++i)
				{
					if constexpr (std::is_invocable_v<F&, T, ParallelChunk<ResultT>*>)
						functionToExecute(this->_list[i], chunk);
					else
						functionToExecute(this->_list[i]);
				}
			});
	}

	static void Execute(const vector<T>& list, function<void(T)> functionToExecute)
	{
		ParallelLoop parallelLoop(list);
		parallelLoop.Execute(functionToExecute);
	}
};

//----------------------------//
//    Number parallel loop    //
//----------------------------//

template <class ResultT = CoeffsChunk>
class NumberParallelLoop : public BaseChunksParallelLoop<ResultT>
{
public:
	NumberParallelLoop(BigNumber endLoop) : NumberParallelLoop(endLoop, BaseParallelLoop::DefaultNThreads) {}

	NumberParallelLoop(BigNumber endLoop, unsigned int nThreads) :
		BaseChunksParallelLoop<ResultT>(endLoop, nThreads)
	{}

	// functionToExecute: void(BigNumber) or void(BigNumber, ParallelChunk<ResultT>*)
	template <class F>
	void Execute(F&& functionToExecute)
	{
		this->ExecuteChunk([&functionToExecute](ParallelChunk<ResultT>* chunk)
			{
				for (BigNumber i = chunk->Start; i < chunk->End; ++i)
				{
					if constexpr (std::is_invocable_v<F&, BigNumber, ParallelChunk<ResultT>*>)
						functionToExecute(i, chunk);
					else
						functionToExecute(i);
				}
			});
	}

	static void Execute(BigNumber endLoop, function<void(BigNumber)> functionToExecute)
	{
		NumberParallelLoop parallelLoop(endLoop);
		parallelLoop.Execute(functionToExecute);
	}
};