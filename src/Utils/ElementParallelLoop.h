#pragma once
#include "ParallelLoop.h"
#include "../Mesh/Element.h"
using namespace std;

template <int Dim>
class ElementParallelLoop : public ParallelLoop<Element<Dim>*, CoeffsChunk>
{
public:
	ElementParallelLoop(vector<Element<Dim>*> list) : ParallelLoop<Element<Dim>*, CoeffsChunk>(list) {}
	ElementParallelLoop(vector<Element<Dim>*> list, unsigned int nThreads) : ParallelLoop<Element<Dim>*, CoeffsChunk>(list, nThreads) {}

	void ReserveChunkCoeffsSize(BigNumber nnzForOneLoopIteration)
	{
		for (unsigned int threadNumber = 0; threadNumber < this->NThreads; threadNumber++)
		{
			ParallelChunk<CoeffsChunk>* chunk = this->Chunks[threadNumber];
			chunk->Results.Coeffs = NonZeroCoefficients(chunk->Size() * nnzForOneLoopIteration);
		}
	}

	void Fill(Eigen::SparseMatrix<double> &m)
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