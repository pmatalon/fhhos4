#pragma once
#include "IDiscreteSpace.h"
#include "../Diff_HHOElement.h"
#include "../../../Utils/ElementParallelLoop.h"


template<int Dim>
class HHOFaceListSpace : public IDiscreteSpace
{
protected:
	int _nUnknowns;

	HHOFaceListSpace() {}

	HHOFaceListSpace(int nUnknowns)
	{
		_nUnknowns = nUnknowns;
	}
private:
	virtual const vector<Face<Dim>*>& ListFaces() = 0;

	virtual Diff_HHOFace<Dim>* HHOFace(Face<Dim>* f) = 0;

	virtual BigNumber Number(Face<Dim>* f) = 0;

public:
	//--------------------------------------------//
	// Implementation of interface IDiscreteSpace //
	//--------------------------------------------//

	BigNumber Dimension() override
	{
		return ListFaces().size() * _nUnknowns;
	}

	Vector InnerProdWithBasis(DomFunction func) override
	{
		Vector innerProds = Vector(Dimension());
		ParallelLoop<Face<Dim>*>::Execute(ListFaces(), [this, &innerProds, func](Face<Dim>* f)
			{
				Diff_HHOFace<Dim>* face = HHOFace(f);
				BigNumber i = Number(f) * _nUnknowns;

				innerProds.segment(i, _nUnknowns) = face->InnerProductWithBasis(func);
			}
		);
		return innerProds;
	}

	Vector ApplyMassMatrix(const Vector& v) override
	{
		assert(v.rows() == Dimension());
		Vector res(v.rows());
		ParallelLoop<Face<Dim>*>::Execute(ListFaces(), [this, &v, &res](Face<Dim>* f)
			{
				Diff_HHOFace<Dim>* face = HHOFace(f);
				BigNumber i = Number(f) * _nUnknowns;

				res.segment(i, _nUnknowns) = face->ApplyMassMatrix(v.segment(i, _nUnknowns));
			});
		return res;
	}

	Vector SolveMassMatrix(const Vector& v) override
	{
		assert(v.rows() == Dimension());
		Vector res(v.rows());
		ParallelLoop<Face<Dim>*>::Execute(ListFaces(), [this, &v, &res](Face<Dim>* f)
			{
				Diff_HHOFace<Dim>* face = HHOFace(f);
				BigNumber i = Number(f) * _nUnknowns;

				res.segment(i, _nUnknowns) = face->SolveMassMatrix(v.segment(i, _nUnknowns));
			});
		return res;
	}

	Vector Project(DomFunction func) override
	{
		Vector vectorOfDoFs = Vector(Dimension());
		ParallelLoop<Face<Dim>*>::Execute(ListFaces(), [this, &vectorOfDoFs, func](Face<Dim>* f)
			{
				Diff_HHOFace<Dim>* face = HHOFace(f);
				BigNumber i = Number(f) * _nUnknowns;

				vectorOfDoFs.segment(i, _nUnknowns) = face->ProjectOnBasis(func);
			}
		);
		return vectorOfDoFs;
	}

	double L2InnerProd(const Vector& v1, const Vector& v2) override
	{
		assert(v1.rows() == Dimension());
		assert(v2.rows() == Dimension());

		struct ChunkResult { double total = 0; };

		ParallelLoop<Face<Dim>*, ChunkResult> parallelLoop(ListFaces());
		parallelLoop.Execute([this, &v1, &v2](Face<Dim>* f, ParallelChunk<ChunkResult>* chunk)
			{
				Diff_HHOFace<Dim>* face = HHOFace(f);
				BigNumber i = Number(f) * _nUnknowns;

				chunk->Results.total += face->InnerProd(v1.segment(i, _nUnknowns), v2.segment(i, _nUnknowns));
			});

		double total = 0;
		parallelLoop.AggregateChunkResults([&total](ChunkResult chunkResult)
			{
				total += chunkResult.total;
			});
		return total;
	}

	double Integral(const Vector& v) override
	{
		assert(v.rows() == Dimension());

		struct ChunkResult { double total = 0; };

		ParallelLoop<Face<Dim>*, ChunkResult> parallelLoop(ListFaces());
		parallelLoop.Execute([this, &v](Face<Dim>* f, ParallelChunk<ChunkResult>* chunk)
			{
				Diff_HHOFace<Dim>* hhoFace = HHOFace(f);
				auto i = Number(f) * _nUnknowns;
				chunk->Results.total += hhoFace->Integral(v.segment(i, _nUnknowns));
			});

		double total = 0;
		parallelLoop.AggregateChunkResults([&total](ChunkResult chunkResult)
			{
				total += chunkResult.total;
			});
		return total;
	}
};
