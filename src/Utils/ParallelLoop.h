#pragma once
#include <vector>
#include <thread>
#include <future>
#include <math.h>
#include <type_traits>
#include "Utils.h"
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
	static unsigned int DefaultNThreads;

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

		// Each thread must have at least 2 elements to process
		while (loopSize / NThreads < 2 && NThreads > 1)
			NThreads--;

		Chunks.reserve(NThreads);
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

	void InitChunkResults(function<void(ResultT&)> functionInitChunkResults)
	{
		for (unsigned int threadNumber = 0; threadNumber < NThreads; threadNumber++)
			functionInitChunkResults(Chunks[threadNumber]->Results);
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

	virtual ~BaseChunksParallelLoop()
	{
		for (unsigned int threadNumber = 0; threadNumber < NThreads; threadNumber++)
		{
			ParallelChunk<ResultT>* chunk = Chunks[threadNumber];
			delete chunk;
		}
	}
};

template <class T, class ResultT = CoeffsChunk>//, typename enable_if<is_base_of<ParallelChunk, ChunkT>::value>::type = ParallelChunk >
class ParallelLoop : public BaseChunksParallelLoop<ResultT>
{
	//static_assert(std::is_base_of<ParallelChunk, ChunkT>::value, "ChunkT must inherit from ParallelChunk");

public:
	vector<T> List;

	ParallelLoop(vector<T> list) :ParallelLoop(list, BaseParallelLoop::DefaultNThreads) {}
	
	ParallelLoop(vector<T> list, unsigned int nThreads) :
		BaseChunksParallelLoop<ResultT>(list.size(), nThreads)
	{
		List = list;
	}

	void Execute(function<void(T)> functionToExecute)
	{
		if (this->NThreads == 1)
		{
			ParallelChunk<ResultT>* chunk = this->Chunks[0];
			for (BigNumber i = chunk->Start; i < chunk->End; ++i)
				functionToExecute(this->List[i]);
		}
		else
		{
			for (unsigned int threadNumber = 0; threadNumber < this->NThreads; threadNumber++)
			{
				ParallelChunk<ResultT>* chunk = this->Chunks[threadNumber];
				chunk->ThreadFuture = std::async(std::launch::async, [chunk, this, functionToExecute]()
					{
						for (BigNumber i = chunk->Start; i < chunk->End; ++i)
							functionToExecute(this->List[i]);
					}
				);
			}
			this->Wait();
		}
	}

	void Execute(function<void(T, ParallelChunk<ResultT>*)> functionToExecute)
	{
		if (this->NThreads == 1)
		{
			ParallelChunk<ResultT>* chunk = this->Chunks[0];
			for (BigNumber i = chunk->Start; i < chunk->End; ++i)
				functionToExecute(this->List[i], chunk);
		}
		else
		{
			for (unsigned int threadNumber = 0; threadNumber < this->NThreads; threadNumber++)
			{
				ParallelChunk<ResultT>* chunk = this->Chunks[threadNumber];
				chunk->ThreadFuture = std::async(std::launch::async, [chunk, this, functionToExecute]()
					{
						for (BigNumber i = chunk->Start; i < chunk->End; ++i)
							functionToExecute(this->List[i], chunk);
					}
				);
			}
			this->Wait();
		}
	}

	static void Execute(vector<T> list, function<void(T)> functionToExecute)
	{
		ParallelLoop parallelLoop(list);
		parallelLoop.Execute(functionToExecute);
	}
};

template <class ResultT = CoeffsChunk>
class NumberParallelLoop : public BaseChunksParallelLoop<ResultT>
{
public:
	NumberParallelLoop(BigNumber endLoop) : NumberParallelLoop(endLoop, BaseParallelLoop::DefaultNThreads) {}

	NumberParallelLoop(BigNumber endLoop, unsigned int nThreads) :
		BaseChunksParallelLoop<ResultT>(endLoop, nThreads)
	{}

	void Execute(function<void(BigNumber)> functionToExecute)
	{
		if (this->NThreads == 1)
		{
			ParallelChunk<ResultT>* chunk = this->Chunks[0];
			for (BigNumber i = chunk->Start; i < chunk->End; ++i)
				functionToExecute(i);
		}
		else
		{
			for (unsigned int threadNumber = 0; threadNumber < this->NThreads; threadNumber++)
			{
				ParallelChunk<ResultT>* chunk = this->Chunks[threadNumber];
				chunk->ThreadFuture = std::async(std::launch::async, [chunk, this, functionToExecute]()
					{
						for (BigNumber i = chunk->Start; i < chunk->End; ++i)
							functionToExecute(i);
					}
				);
			}
			this->Wait();
		}
	}

	void Execute(function<void(BigNumber, ParallelChunk<ResultT>*)> functionToExecute)
	{
		if (this->NThreads == 1)
		{
			ParallelChunk<ResultT>* chunk = this->Chunks[0];
			for (BigNumber i = chunk->Start; i < chunk->End; ++i)
				functionToExecute(i, chunk);
		}
		else
		{
			for (unsigned int threadNumber = 0; threadNumber < this->NThreads; threadNumber++)
			{
				ParallelChunk<ResultT>* chunk = this->Chunks[threadNumber];
				chunk->ThreadFuture = std::async(std::launch::async, [chunk, this, functionToExecute]()
					{
						for (BigNumber i = chunk->Start; i < chunk->End; ++i)
							functionToExecute(i, chunk);
					}
				);
			}
			this->Wait();
		}
	}

	static void Execute(BigNumber endLoop, function<void(BigNumber)> functionToExecute)
	{
		NumberParallelLoop parallelLoop(endLoop);
		parallelLoop.Execute(functionToExecute);
	}
};

unsigned int BaseParallelLoop::DefaultNThreads = std::thread::hardware_concurrency();