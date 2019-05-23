#pragma once
#include <vector>
#include <thread>
#include <future>
#include <math.h>
#include "Utils.h"
using namespace std;

class ParallelChunk
{
public:
	int ThreadNumber;
	BigNumber Start;
	BigNumber End;
	BigNumber Size() { return End - Start; }
	std::future<void> ThreadFuture;
	NonZeroCoefficients Coeffs;

	ParallelChunk(int threadNumber, BigNumber chunkMaxSize, BigNumber loopSize)
	{
		this->ThreadNumber = threadNumber;
		this->Start = threadNumber * chunkMaxSize;
		this->End = min(this->Start + chunkMaxSize, static_cast<long unsigned int>(loopSize));
	}
};

template <class T>
class ParallelLoop
{
public:
	vector<T> List;
	unsigned int NThreads;
	BigNumber ChunkMaxSize;

	vector<ParallelChunk*> Chunks;

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
			Chunks[threadNumber] = new ParallelChunk(threadNumber, ChunkMaxSize, loopSize);
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
			ParallelChunk* chunk = Chunks[0];
			for (BigNumber i = chunk->Start; i < chunk->End; ++i)
				functionToExecute(this->List[i]);
		}
		else
		{
			for (unsigned int threadNumber = 0; threadNumber < NThreads; threadNumber++)
			{
				ParallelChunk* chunk = Chunks[threadNumber];
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

	void Execute(function<void(T, ParallelChunk*)> functionToExecute)
	{
		if (NThreads == 1)
		{
			ParallelChunk* chunk = Chunks[0];
			for (BigNumber i = chunk->Start; i < chunk->End; ++i)
				functionToExecute(this->List[i], chunk);
		}
		else
		{
			for (unsigned int threadNumber = 0; threadNumber < NThreads; threadNumber++)
			{
				ParallelChunk* chunk = Chunks[threadNumber];
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

	void ReserveChunkCoeffsSize(BigNumber nnzForOneLoopIterate)
	{
		for (unsigned int threadNumber = 0; threadNumber < NThreads; threadNumber++)
		{
			ParallelChunk* chunk = Chunks[threadNumber];
			chunk->Coeffs = NonZeroCoefficients(chunk->Size() * nnzForOneLoopIterate);
		}
	}

	inline void Fill(Eigen::SparseMatrix<double> &m)
	{
		NonZeroCoefficients global;
		for (unsigned int threadNumber = 0; threadNumber < NThreads; threadNumber++)
		{
			ParallelChunk* chunk = Chunks[threadNumber];
			global.Add(chunk->Coeffs);
		}
		global.Fill(m);
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
			ParallelChunk* chunk = Chunks[threadNumber];
			delete chunk;
		}
	}

};