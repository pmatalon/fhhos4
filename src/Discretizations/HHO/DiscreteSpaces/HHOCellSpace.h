#pragma once
#include "IDiscreteSpace.h"
#include "../Diff_HHOElement.h"
#include "../../../Utils/ElementParallelLoop.h"


template<int Dim>
class HHOCellSpace : public IDiscreteSpace
{
private:
	const Mesh<Dim>* _mesh = nullptr;
	vector<Diff_HHOElement<Dim>>* _hhoElements = nullptr;
	HHOParameters<Dim>* HHO = nullptr;

public:
	HHOCellSpace()
	{}

	HHOCellSpace(const Mesh<Dim>* mesh, HHOParameters<Dim>* hho, vector<Diff_HHOElement<Dim>>& hhoElements)
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
		return _mesh->Elements.size() * HHO->nCellUnknowns;
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
				BigNumber i = e->Number * HHO->nCellUnknowns;

				innerProds.segment(i, HHO->nCellUnknowns) = elem->InnerProductWithBasis(elem->CellBasis, func);
			}
		);
		return innerProds;
	}

	Vector ApplyMassMatrix(const Vector& v) override
	{
		assert(v.rows() == Dimension());
		Vector res(v.rows());
		ParallelLoop<Element<Dim>*>::Execute(this->_mesh->Elements, [this, &v, &res](Element<Dim>* e)
			{
				Diff_HHOElement<Dim>* elem = HHOElement(e);
				BigNumber i = e->Number * HHO->nCellUnknowns;

				res.segment(i, HHO->nCellUnknowns) = elem->ApplyCellMassMatrix(v.segment(i, HHO->nCellUnknowns));
			});
		return res;
	}

	Vector SolveMassMatrix(const Vector& v) override
	{
		assert(v.rows() == Dimension());
		Vector res(v.rows());
		ParallelLoop<Element<Dim>*>::Execute(this->_mesh->Elements, [this, &v, &res](Element<Dim>* e)
			{
				Diff_HHOElement<Dim>* elem = HHOElement(e);
				BigNumber i = e->Number * HHO->nCellUnknowns;

				res.segment(i, HHO->nCellUnknowns) = elem->SolveCellMassMatrix(v.segment(i, HHO->nCellUnknowns));
			});
		return res;
	}
	
	Vector Project(DomFunction func) override
	{
		Vector vectorOfDoFs = Vector(Dimension());
		ParallelLoop<Element<Dim>*>::Execute(this->_mesh->Elements, [this, &vectorOfDoFs, func](Element<Dim>* e)
			{
				Diff_HHOElement<Dim>* elem = HHOElement(e);
				vectorOfDoFs.segment(e->Number * HHO->nCellUnknowns, HHO->nCellUnknowns) = elem->ProjectOnCellBasis(func);
			}
		);
		return vectorOfDoFs;
	}

	double L2InnerProd(const Vector& v1, const Vector& v2) override
	{
		assert(v1.rows() == Dimension());
		assert(v2.rows() == Dimension());

		struct ChunkResult { double total = 0; };

		ParallelLoop<Element<Dim>*, ChunkResult> parallelLoop(_mesh->Elements);
		parallelLoop.Execute([this, &v1, &v2](Element<Dim>* e, ParallelChunk<ChunkResult>* chunk)
			{
				Diff_HHOElement<Dim>* elem = HHOElement(e);
				BigNumber i = elem->Number() * HHO->nCellUnknowns;

				chunk->Results.total += v1.segment(i, HHO->nCellUnknowns).dot(elem->ApplyCellMassMatrix(v2.segment(i, HHO->nCellUnknowns)));
			});

		double total = 0;
		parallelLoop.AggregateChunkResults([&total](ChunkResult chunkResult)
			{
				total += chunkResult.total;
			});
		return total;
	}

	double Integral(const Vector& cellCoeffs) override
	{
		assert(cellCoeffs.rows() == Dimension());

		struct ChunkResult { double total = 0; };

		ParallelLoop<Element<Dim>*, ChunkResult> parallelLoop(_mesh->Elements);
		parallelLoop.Execute([this, &cellCoeffs](Element<Dim>* e, ParallelChunk<ChunkResult>* chunk)
			{
				Diff_HHOElement<Dim>* hhoElem = HHOElement(e);
				auto i = e->Number * HHO->nCellUnknowns;
				chunk->Results.total += hhoElem->IntegralCell(cellCoeffs.segment(i, HHO->nCellUnknowns));
			});

		double total = 0;
		parallelLoop.AggregateChunkResults([&total](ChunkResult chunkResult)
			{
				total += chunkResult.total;
			});
		return total;
	}
};
