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
	future<void> ThreadFuture;
	ParallelChunk(int threadNumber, BigNumber chunkMaxSize, BigNumber loopSize)
	{
		this->ThreadNumber = threadNumber;
		this->Start = threadNumber * chunkMaxSize;
		this->End = min(this->Start + chunkMaxSize, static_cast<long unsigned int>(loopSize));
	}
};

class ParallelLoop
{
public:
	unsigned int NThreads;
	BigNumber ChunkMaxSize;

	vector<ParallelChunk*> Chunks;

	ParallelLoop(BigNumber loopSize)
	{
		NThreads = std::thread::hardware_concurrency();
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

	template <class T>
	static void Execute(vector<T> list, function<void(T)> functionToExecute)
	{
		ParallelLoop parallelLoop(list.size());

		for (unsigned int threadNumber = 0; threadNumber < parallelLoop.NThreads; threadNumber++)
		{
			ParallelChunk* chunk = parallelLoop.Chunks[threadNumber];
			chunk->ThreadFuture = std::async([chunk, list, functionToExecute]()
				{
					for (BigNumber i = chunk->Start; i < chunk->End; ++i)
					{
						functionToExecute(list[i]);
					}
				}
			);
		}
		parallelLoop.Wait();
	}

	~ParallelLoop()
	{
		for (auto chunk : this->Chunks)
			delete chunk;
	}

};