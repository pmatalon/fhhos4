#pragma once
#include "HHOSkeletonSpace.h"

// Same as HHOSkeletonSpace<Dim>,
// except that we filter Dirichlet faces in every function.

template<int Dim>
class HHONonDirichletFaceSpace : public HHOSkeletonSpace<Dim>
{
public:
	HHONonDirichletFaceSpace()
	{}

	HHONonDirichletFaceSpace(const Mesh<Dim>* mesh, HHOParameters<Dim>* hho, vector<Diff_HHOFace<Dim>>& hhoFaces)
		: HHOSkeletonSpace<Dim>(mesh, hho, hhoFaces)
	{}

public:
	//--------------------------------------------//
	// Implementation of interface IDiscreteSpace //
	//--------------------------------------------//

	BigNumber Dimension() override
	{
		return this->HHO->nInteriorAndNeumannFaces * this->HHO->nFaceUnknowns;
	}

	double Measure() override
	{
		double measure = 0;
		for (Face<Dim>* f : this->ListFaces())
		{
			if (f->HasDirichletBC())
				continue;
			measure += f->Measure();
		}
		return measure;
	}

	Vector InnerProdWithBasis(DomFunction func) override
	{
		Vector innerProds = Vector(Dimension());
		ParallelLoop<Face<Dim>*>::Execute(this->ListFaces(), [this, &innerProds, func](Face<Dim>* f)
			{
				if (f->HasDirichletBC())
					return;
				BigNumber i = this->Number(f) * this->_nUnknowns;
				innerProds.segment(i, this->_nUnknowns) = this->HHOFace(f)->InnerProductWithBasis(func);
			}
		);
		return innerProds;
	}

	Vector ApplyMassMatrix(const Vector& v) override
	{
		assert(v.rows() == Dimension());
		Vector res(v.rows());
		ParallelLoop<Face<Dim>*>::Execute(this->ListFaces(), [this, &v, &res](Face<Dim>* f)
			{
				if (f->HasDirichletBC())
					return;
				BigNumber i = this->Number(f) * this->_nUnknowns;
				res.segment(i, this->_nUnknowns) = this->HHOFace(f)->ApplyMassMatrix(v.segment(i, this->_nUnknowns));
			});
		return res;
	}

	Vector SolveMassMatrix(const Vector& v) override
	{
		assert(v.rows() == Dimension());
		Vector res(v.rows());
		ParallelLoop<Face<Dim>*>::Execute(this->ListFaces(), [this, &v, &res](Face<Dim>* f)
			{
				if (f->HasDirichletBC())
					return;
				BigNumber i = this->Number(f) * this->_nUnknowns;
				res.segment(i, this->_nUnknowns) = this->HHOFace(f)->SolveMassMatrix(v.segment(i, this->_nUnknowns));
			});
		return res;
	}

	Vector Project(DomFunction func) override
	{
		Vector vectorOfDoFs = Vector(Dimension());
		ParallelLoop<Face<Dim>*>::Execute(this->ListFaces(), [this, &vectorOfDoFs, func](Face<Dim>* f)
			{
				if (f->HasDirichletBC())
					return;
				BigNumber i = this->Number(f) * this->_nUnknowns;
				vectorOfDoFs.segment(i, this->_nUnknowns) = this->HHOFace(f)->ProjectOnBasis(func);
			}
		);
		return vectorOfDoFs;
	}
	
	double L2InnerProd(const Vector& v1, const Vector& v2) override
	{
		assert(v1.rows() == Dimension());
		assert(v2.rows() == Dimension());

		struct ChunkResult { double total = 0; };

		ParallelLoop<Face<Dim>*, ChunkResult> parallelLoop(this->ListFaces());
		parallelLoop.Execute([this, &v1, &v2](Face<Dim>* f, ParallelChunk<ChunkResult>* chunk)
			{
				if (f->HasDirichletBC())
					return;
				BigNumber i = this->Number(f) * this->_nUnknowns;
				chunk->Results.total += this->HHOFace(f)->InnerProd(v1.segment(i, this->_nUnknowns), v2.segment(i, this->_nUnknowns));
			});

		double total = 0;
		parallelLoop.AggregateChunkResults([&total](ChunkResult chunkResult)
			{
				total += chunkResult.total;
			});
		return total;
	}

	double Integral(const Vector& faceCoeffs) override
	{
		assert(faceCoeffs.rows() == Dimension());

		struct ChunkResult { double total = 0; };

		ParallelLoop<Face<Dim>*, ChunkResult> parallelLoop(this->ListFaces());
		parallelLoop.Execute([this, &faceCoeffs](Face<Dim>* f, ParallelChunk<ChunkResult>* chunk)
			{
				if (f->HasDirichletBC())
					return;
				auto i = this->Number(f) * this->_nUnknowns;
				chunk->Results.total += this->HHOFace(f)->Integral(faceCoeffs.segment(i, this->_nUnknowns));
			});

		double total = 0;
		parallelLoop.AggregateChunkResults([&total](ChunkResult chunkResult)
			{
				total += chunkResult.total;
			});
		return total;
	}

	double Integral(DomFunction func) override
	{
		struct ChunkResult { double total = 0; };

		ParallelLoop<Face<Dim>*, ChunkResult> parallelLoop(this->ListFaces());
		parallelLoop.Execute([this, func](Face<Dim>* f, ParallelChunk<ChunkResult>* chunk)
			{
				if (f->HasDirichletBC())
					return;
				chunk->Results.total += f->Integral(func);
			});

		double total = 0;
		parallelLoop.AggregateChunkResults([&total](ChunkResult chunkResult)
			{
				total += chunkResult.total;
			});
		return total;
	}
};
