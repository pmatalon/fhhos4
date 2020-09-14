#pragma once
#include "ParallelLoop.h"
#include "../Mesh/Element.h"
using namespace std;

template <class T>
class ParallelLoopFillingOneMatrix : public ParallelLoop<T, CoeffsChunk>
{
public:
	ParallelLoopFillingOneMatrix(const vector<T>& list) : ParallelLoop<T, CoeffsChunk>(list) {}
	ParallelLoopFillingOneMatrix(const vector<T>& list, unsigned int nThreads) : ParallelLoop<T, CoeffsChunk>(list, nThreads) {}

	void ReserveChunkCoeffsSize(BigNumber nnzForOneLoopIteration)
	{
		for (unsigned int threadNumber = 0; threadNumber < this->NThreads; threadNumber++)
		{
			ParallelChunk<CoeffsChunk>* chunk = this->Chunks[threadNumber];
			chunk->Results.Coeffs = NonZeroCoefficients(chunk->Size() * nnzForOneLoopIteration);
		}
	}

	void Fill(SparseMatrix &m)
	{
		NonZeroCoefficients global;
		for (unsigned int threadNumber = 0; threadNumber < this->NThreads; threadNumber++)
		{
			ParallelChunk<CoeffsChunk>* chunk = this->Chunks[threadNumber];
			global.Add(chunk->Results.Coeffs);
		}
		global.Fill(m);
	}
};

template <int Dim>
class ElementParallelLoop : public ParallelLoopFillingOneMatrix<Element<Dim>*>
{
public:
	ElementParallelLoop(const vector<Element<Dim>*>& list) : ParallelLoopFillingOneMatrix<Element<Dim>*>(list) {}
	ElementParallelLoop(const vector<Element<Dim>*>& list, unsigned int nThreads) : ParallelLoopFillingOneMatrix<Element<Dim>*>(list, nThreads) {}
};

template <int Dim>
class FaceParallelLoop : public ParallelLoopFillingOneMatrix<Face<Dim>*>
{
public:
	FaceParallelLoop(const vector<Face<Dim>*>& list) : ParallelLoopFillingOneMatrix<Face<Dim>*>(list) {}
	FaceParallelLoop(const vector<Face<Dim>*>& list, unsigned int nThreads) : ParallelLoopFillingOneMatrix<Face<Dim>*>(list, nThreads) {}
};