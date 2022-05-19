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
		Vector innerProds = Vector(HHO->nTotalReconstructUnknowns);
		ParallelLoop<Element<Dim>*>::Execute(this->_mesh->Elements, [this, &innerProds, func](Element<Dim>* e)
			{
				Diff_HHOElement<Dim>* elem = HHOElement(e);
				BigNumber i = e->Number * HHO->nReconstructUnknowns;

				innerProds.segment(i, HHO->nReconstructUnknowns) = elem->InnerProductWithBasis(elem->ReconstructionBasis, func);
			}
		);
		return innerProds;
	}

	Vector SolveMassMatrix(const Vector& v) override
	{
		assert(v.rows() == HHO->nTotalReconstructUnknowns);
		Vector res(v.rows());
		ParallelLoop<Element<Dim>*>::Execute(this->_mesh->Elements, [this, &v, &res](Element<Dim>* e)
			{
				Diff_HHOElement<Dim>* elem = HHOElement(e);
				BigNumber i = e->Number * HHO->nReconstructUnknowns;

				res.segment(i, HHO->nReconstructUnknowns) = elem->SolveReconstructMassMatrix(v.segment(i, HHO->nReconstructUnknowns));
			});
		return res;
	}
	
	Vector Project(DomFunction func) override
	{
		Vector vectorOfDoFs = Vector(HHO->nTotalReconstructUnknowns);
		ParallelLoop<Element<Dim>*>::Execute(this->_mesh->Elements, [this, &vectorOfDoFs, func](Element<Dim>* e)
			{
				Diff_HHOElement<Dim>* elem = HHOElement(e);
				vectorOfDoFs.segment(e->Number * HHO->nReconstructUnknowns, HHO->nReconstructUnknowns) = elem->ProjectOnReconstructBasis(func);
			}
		);
		return vectorOfDoFs;
	}

	double L2InnerProd(const Vector& v1, const Vector& v2) override
	{
		assert(v1.rows() == _mesh->Elements.size() * HHO->nReconstructUnknowns);
		assert(v2.rows() == _mesh->Elements.size() * HHO->nReconstructUnknowns);

		struct ChunkResult { double total = 0; };

		ParallelLoop<Element<Dim>*, ChunkResult> parallelLoop(_mesh->Elements);
		parallelLoop.Execute([this, &v1, &v2](Element<Dim>* e, ParallelChunk<ChunkResult>* chunk)
			{
				Diff_HHOElement<Dim>* elem = HHOElement(e);
				BigNumber i = elem->Number() * HHO->nReconstructUnknowns;

				chunk->Results.total += v1.segment(i, HHO->nReconstructUnknowns).dot(elem->ApplyReconstructMassMatrix(v2.segment(i, HHO->nReconstructUnknowns)));
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
		assert(reconstructedCoeffs.rows() == HHO->nTotalReconstructUnknowns);

		struct ChunkResult { double total = 0; };

		ParallelLoop<Element<Dim>*, ChunkResult> parallelLoop(_mesh->Elements);
		parallelLoop.Execute([this, &reconstructedCoeffs](Element<Dim>* e, ParallelChunk<ChunkResult>* chunk)
			{
				Diff_HHOElement<Dim>* hhoElem = HHOElement(e);
				auto i = e->Number * HHO->nReconstructUnknowns;
				chunk->Results.total += hhoElem->IntegralReconstruct(reconstructedCoeffs.segment(i, HHO->nReconstructUnknowns));
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
