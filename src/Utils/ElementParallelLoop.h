#pragma once
#include "ParallelLoop.h"
#include "../Mesh/Element.h"
using namespace std;

/*template <int Dim, class ResultT>
class ElementParallelLoop : public ParallelLoop<Element<Dim>*, ResultT>
{
public:
	ElementParallelLoop(vector<Element<Dim>*> list) : ParallelLoop<Element<Dim>*, ResultT>(list) {}
	ElementParallelLoop(vector<Element<Dim>*> list, unsigned int nThreads) : ParallelLoop<Element<Dim>*, ResultT>(list, nThreads) {}

	void ReserveChunkCoeffsSize(BigNumber nnzForOneLoopIteration)
	{	}

	void Fill(Eigen::SparseMatrix<double> &m)
	{	}
};

template <int Dim>
void ElementParallelLoop<Dim, CoeffsChunk>::ReserveChunkCoeffsSize(BigNumber nnzForOneLoopIteration)
{
	for (unsigned int threadNumber = 0; threadNumber < this->NThreads; threadNumber++)
	{
		ParallelChunk<CoeffsChunk>* chunk = this->Chunks[threadNumber];
		chunk->Results.Coeffs = NonZeroCoefficients(chunk->Size() * nnzForOneLoopIteration);
	}
}

template <int Dim>
void ElementParallelLoop<Dim, CoeffsChunk>::Fill(Eigen::SparseMatrix<double> &m)
{
	NonZeroCoefficients global;
	for (unsigned int threadNumber = 0; threadNumber < this->NThreads; threadNumber++)
	{
		ParallelChunk<CoeffsChunk>* chunk = this->Chunks[threadNumber];
		global.Add(chunk->Results.Coeffs);
	}
	global.Fill(m);
}*/

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
		//cout << "Fill" << endl;
		NonZeroCoefficients global;
		for (unsigned int threadNumber = 0; threadNumber < this->NThreads; threadNumber++)
		{
			ParallelChunk<CoeffsChunk>* chunk = this->Chunks[threadNumber];
			global.Add(chunk->Results.Coeffs);
		}
		global.Fill(m);
	}
};