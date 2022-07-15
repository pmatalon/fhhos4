#pragma once
#include "IDiscreteSpace.h"
#include "../Diff_HHOElement.h"
#include "../../../Utils/ElementParallelLoop.h"


template<int Dim>
class HHOSkeletonSpace : public IDiscreteSpace
{
private:
	const Mesh<Dim>* _mesh = nullptr;
	vector<Diff_HHOFace<Dim>>* _hhoFaces = nullptr;
	HHOParameters<Dim>* HHO = nullptr;

public:
	HHOSkeletonSpace()
	{}

	HHOSkeletonSpace(const Mesh<Dim>* mesh, HHOParameters<Dim>* hho, vector<Diff_HHOFace<Dim>>& hhoFaces)
	{
		HHO = hho;
		_mesh = mesh;
		_hhoFaces = &hhoFaces;
	}

private:
	const vector<Face<Dim>*>& ListFaces()
	{
		return _mesh->Faces;
	}

	Diff_HHOFace<Dim>* HHOFace(Face<Dim>* f)
	{
		return &(*_hhoFaces)[f->Number];
	}

	BigNumber Number(Face<Dim>* f)
	{
		return f->Number;
	}

public:
	//--------------------------------------------//
	// Implementation of interface IDiscreteSpace //
	//--------------------------------------------//

	BigNumber Dimension() override
	{
		return _mesh->Faces.size() * HHO->nFaceUnknowns;
	}

	double Measure() override
	{
		return _mesh->SkeletonMeasure();
	}

	Vector InnerProdWithBasis(DomFunction func) override
	{
		Vector innerProds = Vector(Dimension());
		ParallelLoop<Face<Dim>*>::Execute(ListFaces(), [this, &innerProds, func](Face<Dim>* f)
			{
				Diff_HHOFace<Dim>* face = HHOFace(f);
				BigNumber i = Number(f) * HHO->nFaceUnknowns;

				innerProds.segment(i, HHO->nFaceUnknowns) = face->InnerProductWithBasis(func);
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
				BigNumber i = Number(f) * HHO->nFaceUnknowns;

				res.segment(i, HHO->nFaceUnknowns) = face->ApplyMassMatrix(v.segment(i, HHO->nFaceUnknowns));
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
				BigNumber i = Number(f) * HHO->nFaceUnknowns;

				res.segment(i, HHO->nFaceUnknowns) = face->SolveMassMatrix(v.segment(i, HHO->nFaceUnknowns));
			});
		return res;
	}

	Vector Project(DomFunction func) override
	{
		Vector vectorOfDoFs = Vector(HHO->nTotalFaceUnknowns);
		ParallelLoop<Face<Dim>*>::Execute(ListFaces(), [this, &vectorOfDoFs, func](Face<Dim>* f)
			{
				if (f->HasDirichletBC()) // CAREFUL: just the non-Dirichlet faces
					return;
				Diff_HHOFace<Dim>* face = HHOFace(f);
				BigNumber i = Number(f) * HHO->nFaceUnknowns;

				vectorOfDoFs.segment(i, HHO->nFaceUnknowns) = face->ProjectOnBasis(func);
			}
		);
		return vectorOfDoFs;
	}
	
	double L2InnerProd(const Vector& v1, const Vector& v2) override
	{
		Utils::FatalError("To be implemented");
		return 0;
	}

	double Integral(const Vector& faceCoeffs) override
	{
		assert(faceCoeffs.rows() == Dimension());

		struct ChunkResult { double total = 0; };

		ParallelLoop<Face<Dim>*, ChunkResult> parallelLoop(ListFaces());
		parallelLoop.Execute([this, &faceCoeffs](Face<Dim>* f, ParallelChunk<ChunkResult>* chunk)
			{
				Diff_HHOFace<Dim>* hhoFace = HHOFace(f);
				auto i = Number(f) * HHO->nFaceUnknowns;
				chunk->Results.total += hhoFace->Integral(faceCoeffs.segment(i, HHO->nFaceUnknowns));
			});

		double total = 0;
		parallelLoop.AggregateChunkResults([&total](ChunkResult chunkResult)
			{
				total += chunkResult.total;
			});
		return total;
	}
};
