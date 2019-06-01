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
	BigNumber Start;
	BigNumber End;
	BigNumber Size() { return End - Start; }
	std::future<void> ThreadFuture;
	ResultT Results;

	ParallelChunk(int threadNumber, BigNumber chunkMaxSize, BigNumber loopSize)
	{
		this->ThreadNumber = threadNumber;
		this->Start = threadNumber * chunkMaxSize;
		this->End = min(this->Start + chunkMaxSize, static_cast<long unsigned int>(loopSize));
	}
};

template <class T, class ResultT = CoeffsChunk>//, typename enable_if<is_base_of<ParallelChunk, ChunkT>::value>::type = ParallelChunk >
class ParallelLoop
{
	//static_assert(std::is_base_of<ParallelChunk, ChunkT>::value, "ChunkT must inherit from ParallelChunk");

public:
	vector<T> List;
	unsigned int NThreads;
	BigNumber ChunkMaxSize;

	vector<ParallelChunk<ResultT>*> Chunks;

	ParallelLoop(vector<T> list) :ParallelLoop(list, std::thread::hardware_concurrency()) {}
	
	ParallelLoop(vector<T> list, unsigned int nThreads)
	{
		List = list;
		BigNumber loopSize = list.size();
		NThreads = nThreads;
		if (NThreads == 0)
			NThreads = 1;
		if (NThreads > loopSize)
			NThreads = loopSize;

		Chunks.reserve(NThreads);
		ChunkMaxSize = (BigNumber)ceil(loopSize / (double)NThreads);

		for (unsigned int threadNumber = 0; threadNumber < NThreads; threadNumber++)
			Chunks[threadNumber] = new ParallelChunk<ResultT>(threadNumber, ChunkMaxSize, loopSize);
	}

	void InitChunkResults(function<void(ResultT)> functionInitChunkResults)
	{
		for (unsigned int threadNumber = 0; threadNumber < NThreads; threadNumber++)
		{
			ResultT results = Chunks[threadNumber]->Results;
			functionInitChunkResults(results);
		}
	}
	
	void AggregateChunkResults(function<void(ResultT)> aggregate)
	{
		for (unsigned int threadNumber = 0; threadNumber < NThreads; threadNumber++)
		{
			ResultT results = Chunks[threadNumber]->Results;
			aggregate(results);
		}
	}

	void Wait()
	{
		for (unsigned int threadNumber = 0; threadNumber < NThreads; threadNumber++)
			this->Chunks[threadNumber]->ThreadFuture.wait();
	}

	void Execute(function<void(T)> functionToExecute)
	{
		if (NThreads == 1)
		{
			ParallelChunk<ResultT>* chunk = Chunks[0];
			for (BigNumber i = chunk->Start; i < chunk->End; ++i)
				functionToExecute(this->List[i]);
		}
		else
		{
			for (unsigned int threadNumber = 0; threadNumber < NThreads; threadNumber++)
			{
				ParallelChunk<ResultT>* chunk = Chunks[threadNumber];
				chunk->ThreadFuture = std::async(std::launch::async, [chunk, this, functionToExecute]()
					{
						for (BigNumber i = chunk->Start; i < chunk->End; ++i)
							functionToExecute(this->List[i]);
					}
				);
			}
			Wait();
		}
	}

	void Execute(function<void(T, ParallelChunk<ResultT>*)> functionToExecute)
	{
		if (NThreads == 1)
		{
			ParallelChunk<ResultT>* chunk = Chunks[0];
			for (BigNumber i = chunk->Start; i < chunk->End; ++i)
				functionToExecute(this->List[i], chunk);
		}
		else
		{
			for (unsigned int threadNumber = 0; threadNumber < NThreads; threadNumber++)
			{
				ParallelChunk<ResultT>* chunk = Chunks[threadNumber];
				chunk->ThreadFuture = std::async(std::launch::async, [chunk, this, functionToExecute]()
					{
						for (BigNumber i = chunk->Start; i < chunk->End; ++i)
							functionToExecute(this->List[i], chunk);
					}
				);
			}
			Wait();
		}
	}

	static void Execute(vector<T> list, function<void(T)> functionToExecute)
	{
		ParallelLoop parallelLoop(list);
		parallelLoop.Execute(functionToExecute);
	}

	~ParallelLoop()
	{
		for (unsigned int threadNumber = 0; threadNumber < NThreads; threadNumber++)
		{
			ParallelChunk<ResultT>* chunk = Chunks[threadNumber];
			delete chunk;
		}
	}

};

/*template <>
void ParallelLoop<, CoeffsChunk>::Fill(Eigen::SparseMatrix<double> &m)
{
	NonZeroCoefficients global;
	for (unsigned int threadNumber = 0; threadNumber < NThreads; threadNumber++)
	{
		ParallelChunk<CoeffsChunk>* chunk = Chunks[threadNumber];
		global.Add(chunk->Results.Coeffs);
	}
	global.Fill(m);
}*/