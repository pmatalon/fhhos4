#pragma once
#include "IDiscreteSpace.h"
#include "../Diff_HHOElement.h"
#include "../../../Utils/ElementParallelLoop.h"


template<int Dim>
class HHOReconstructSpace : public IDiscreteSpace
{
private:
	const Mesh<Dim>* _mesh = nullptr;
	vector<Diff_HHOElement<Dim>>* _hhoElements = nullptr;
	HHOParameters<Dim>* HHO = nullptr;

public:
	HHOReconstructSpace()
	{}

	HHOReconstructSpace(const Mesh<Dim>* mesh, HHOParameters<Dim>* hho, vector<Diff_HHOElement<Dim>>& hhoElements)
	{
		HHO = hho;
		_mesh = mesh;
		_hhoElements = &hhoElements;
	}

	Diff_HHOElement<Dim>* HHOElement(Element<Dim>* e)
	{
		return &(*_hhoElements)[e->Number];
	}

	//--------------------------------------------//
	// Implementation of interface IDiscreteSpace //
	//--------------------------------------------//
	
	BigNumber Dimension() override
	{
		return _mesh->Elements.size() * HHO->nReconstructUnknowns;
	}

	double Measure() override
	{
		return _mesh->Measure();
	}

	Vector InnerProdWithBasis(DomFunction func) override
	{
		Vector innerProds = Vector(Dimension());
		ParallelLoop<Element<Dim>*>::Execute(this->_mesh->Elements, [this, &innerProds, func](Element<Dim>* e)
			{
				Diff_HHOElement<Dim>* elem = HHOElement(e);
				BigNumber i = e->Number * HHO->nReconstructUnknowns;

				innerProds.segment(i, HHO->nReconstructUnknowns) = elem->InnerProductWithBasis(elem->ReconstructionBasis, func);
			}
		);
		return innerProds;
	}

	Vector ApplyMassMatrix(const Vector& v) override
	{
		assert(v.rows() == Dimension());

		if (HHO->OrthonormalizeElemBases())
			return v;

		Vector res(v.rows());
		ParallelLoop<Element<Dim>*>::Execute(this->_mesh->Elements, [this, &v, &res](Element<Dim>* e)
			{
				BigNumber i = e->Number * HHO->nReconstructUnknowns;
				res.segment(i, HHO->nReconstructUnknowns) = HHOElement(e)->ApplyReconstructMassMatrix(v.segment(i, HHO->nReconstructUnknowns));
			});
		return res;
	}

	Vector SolveMassMatrix(const Vector& v) override
	{
		assert(v.rows() == Dimension());

		if (HHO->OrthonormalizeElemBases())
			return v;

		Vector res(v.rows());
		ParallelLoop<Element<Dim>*>::Execute(this->_mesh->Elements, [this, &v, &res](Element<Dim>* e)
			{
				BigNumber i = e->Number * HHO->nReconstructUnknowns;
				res.segment(i, HHO->nReconstructUnknowns) = HHOElement(e)->SolveReconstructMassMatrix(v.segment(i, HHO->nReconstructUnknowns));
			});
		return res;
	}
	
	Vector Project(DomFunction func) override
	{
		Vector vectorOfDoFs = Vector(Dimension());
		ParallelLoop<Element<Dim>*>::Execute(this->_mesh->Elements, [this, &vectorOfDoFs, func](Element<Dim>* e)
			{
				BigNumber i = e->Number * HHO->nReconstructUnknowns;
				vectorOfDoFs.segment(i, HHO->nReconstructUnknowns) = HHOElement(e)->ProjectOnReconstructBasis(func);
			}
		);
		return vectorOfDoFs;
	}

	double L2InnerProd(const Vector& v1, const Vector& v2) override
	{
		assert(v1.rows() == Dimension());
		assert(v2.rows() == Dimension());

		if (HHO->OrthonormalizeElemBases())
			return v1.dot(v2);

		struct ChunkResult { double total = 0; };

		ParallelLoop<Element<Dim>*, ChunkResult> parallelLoop(_mesh->Elements);
		parallelLoop.Execute([this, &v1, &v2](Element<Dim>* e, ParallelChunk<ChunkResult>* chunk)
			{
				BigNumber i = e->Number * HHO->nReconstructUnknowns;
				chunk->Results.total += v1.segment(i, HHO->nReconstructUnknowns).dot(HHOElement(e)->ApplyReconstructMassMatrix(v2.segment(i, HHO->nReconstructUnknowns)));
			});

		double total = 0;
		parallelLoop.AggregateChunkResults([&total](ChunkResult chunkResult)
			{
				total += chunkResult.total;
			});
		return total;
	}

	double Integral(const Vector& reconstructedCoeffs) override
	{
		assert(reconstructedCoeffs.rows() == Dimension());

		struct ChunkResult { double total = 0; };

		ParallelLoop<Element<Dim>*, ChunkResult> parallelLoop(_mesh->Elements);
		parallelLoop.Execute([this, &reconstructedCoeffs](Element<Dim>* e, ParallelChunk<ChunkResult>* chunk)
			{
				auto i = e->Number * HHO->nReconstructUnknowns;
				chunk->Results.total += HHOElement(e)->IntegralReconstruct(reconstructedCoeffs.segment(i, HHO->nReconstructUnknowns));
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

		ParallelLoop<Element<Dim>*, ChunkResult> parallelLoop(_mesh->Elements);
		parallelLoop.Execute([this, func](Element<Dim>* e, ParallelChunk<ChunkResult>* chunk)
			{
				chunk->Results.total += e->Integral(func);
			});

		double total = 0;
		parallelLoop.AggregateChunkResults([&total](ChunkResult chunkResult)
			{
				total += chunkResult.total;
			});
		return total;
	}

	/*
	SparseMatrix StiffnessMatrix()
	{
		ElementParallelLoop<Dim> parallelLoop(_mesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(HHO->nReconstructUnknowns * HHO->nReconstructUnknowns);

		parallelLoop.Execute([this](Element<Dim>* e, ParallelChunk<CoeffsChunk>* chunk)
			{
				Diff_HHOElement<Dim>* elem = this->HHOElement(e);

				chunk->Results.Coeffs.Add(e->Number * HHO->nReconstructUnknowns, e->Number * HHO->nReconstructUnknowns, e->IntegralGradGradMatrix(elem->ReconstructionBasis));
			});

		SparseMatrix stiff(HHO->nTotalReconstructUnknowns, HHO->nTotalReconstructUnknowns);
		parallelLoop.Fill(stiff);
		return stiff;
	}*/
};
