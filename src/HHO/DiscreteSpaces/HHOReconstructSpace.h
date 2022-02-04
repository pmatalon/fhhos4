#pragma once
#include "IDiscreteSpace.h"
#include "../Diff_HHOElement.h"
#include "../../Utils/ElementParallelLoop.h"


template<int Dim>
class HHOReconstructSpace : public IDiscreteSpace
{
private:
	Mesh<Dim>* _mesh = nullptr;
	vector<Diff_HHOElement<Dim>>* _hhoElements = nullptr;
	HHOParameters<Dim>* HHO = nullptr;

public:
	HHOReconstructSpace()
	{}

	HHOReconstructSpace(Mesh<Dim>* mesh, HHOParameters<Dim>* hho, vector<Diff_HHOElement<Dim>>& hhoElements)
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

				innerProds.segment(i, HHO->nReconstructUnknowns) = elem->InnerProductWithReconstructBasis(func);
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
		Utils::FatalError("To be implemented");
	}

	double L2InnerProd(const Vector& v1, const Vector& v2) override
	{
		Utils::FatalError("To be implemented");
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
};
