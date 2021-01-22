#pragma once
#include "../Utils/Utils.h"
#include <mutex>
using namespace std;

struct ElementAggregate;
struct AlgebraicElement
{
	BigNumber Number;
	vector<AlgebraicElement*> Neighbours;
	bool IsAggregated = false;
	ElementAggregate* CoarseElement = nullptr;
};


struct ElementAggregate
{
	BigNumber Number;
	vector<AlgebraicElement*> FineElements;

	ElementAggregate(BigNumber number, vector<AlgebraicElement*> elements)
		: Number(number), FineElements(elements) 
	{}
};

class AlgebraicMesh
{
private:
	int _blockSize;

	const SparseMatrix* A;
public:
	vector<AlgebraicElement> _elements;
	vector<ElementAggregate> _coarseElements;

public:
	AlgebraicMesh(int blockSize)
	{
		this->_blockSize = blockSize;
	}

	void Build(const SparseMatrix& A)
	{
		this->A = &A;

		BigNumber nElements = A.rows() / _blockSize;
		this->_elements = vector<AlgebraicElement>(nElements);

		NumberParallelLoop<EmptyResultChunk> parallelLoopElem(_elements.size());
		parallelLoopElem.Execute([this, &A](BigNumber elemNumber, ParallelChunk<EmptyResultChunk>* chunk)
			{
				AlgebraicElement& elem = _elements[elemNumber];
				elem.Number = elemNumber;

				// RowMajor --> the following line iterates over the non-zeros of the elemNumber-th row.
				for (SparseMatrix::InnerIterator it(A, elemNumber*_blockSize); it; ++it)
				{
					BigNumber neighbourNumber = it.col() / _blockSize;
					if (neighbourNumber != elemNumber)
					{
						AlgebraicElement* neighbour = &_elements[neighbourNumber];
						elem.Neighbours.push_back(neighbour);
					}
				}
			});
	}

	void PairWiseAggregate(bool& coarsestPossibleMeshReached)
	{
		// Element aggregation
		_coarseElements.reserve(_elements.size() / 2);

		for (BigNumber i = 0; i < _elements.size(); i++)
		{
			AlgebraicElement& elem = _elements[i];
			if (elem.IsAggregated)
				continue;

			if (elem.Neighbours.empty())
			{
				coarsestPossibleMeshReached = true;
				return;
			}

			AlgebraicElement* strongestNeighbour = StrongestNeighbour(elem, true);
			// If no neighbour available, get the already aggregated strongest neighbour
			if (!strongestNeighbour)
				strongestNeighbour = StrongestNeighbour(elem, false);
			
			if (!strongestNeighbour)
			{
				coarsestPossibleMeshReached = true;
				return;
			}

			Aggregate(elem, *strongestNeighbour);
		}
	}

private:
	AlgebraicElement* StrongestNeighbour(const AlgebraicElement& e, bool checkAvailability)
	{
		AlgebraicElement* strongestNeighbour = nullptr;
		double strongestNegativeCoupling = 0;
		int smallestAggregateSize = 1000000;
		//DenseMatrix elemBlock = A->block(e.Number*_blockSize, e.Number*_blockSize, _blockSize, _blockSize);
		//double elemKappa = elemBlock(0, 0);
		for (AlgebraicElement* n : e.Neighbours)
		{
			if (n == &e || (checkAvailability && n->IsAggregated))
				continue;

			/*DenseMatrix neighourBlock = A->block(n->Number*_blockSize, n->Number*_blockSize, _blockSize, _blockSize);
			double neighbourKappa = neighourBlock(0, 0);
			double minKappa = min(elemKappa, neighbourKappa);
			double maxKappa = max(elemKappa, neighbourKappa);
			if (maxKappa/minKappa > 10)
				continue;*/

			DenseMatrix couplingBlock = A->block(e.Number*_blockSize, n->Number*_blockSize, _blockSize, _blockSize);
			double coupling = couplingBlock(0, 0);
			if (coupling < strongestNegativeCoupling)
			{
				if (checkAvailability)
				{
					strongestNegativeCoupling = coupling;
					strongestNeighbour = n;
				}
				else
				{
					int aggregateSize = n->CoarseElement->FineElements.size();
					if (aggregateSize < smallestAggregateSize)
					{
						strongestNegativeCoupling = coupling;
						strongestNeighbour = n;
						smallestAggregateSize = aggregateSize;
					}
				}
			}
		}
		return strongestNeighbour;
	}

	
	void Aggregate(AlgebraicElement& e1, AlgebraicElement& e2)
	{
		assert(!(e1.IsAggregated && e2.IsAggregated));

		if (!e1.IsAggregated && !e2.IsAggregated)
		{
			_coarseElements.push_back(ElementAggregate(_coarseElements.size(), { &e1, &e2 }));
			e1.IsAggregated = true;
			e2.IsAggregated = true;
			e1.CoarseElement = &_coarseElements.back();
			e2.CoarseElement = &_coarseElements.back();
		}
		else if (e1.IsAggregated)
		{
			e1.CoarseElement->FineElements.push_back(&e2);
			e2.IsAggregated = true;
			e2.CoarseElement = e1.CoarseElement;
		}
		else
			Aggregate(e2, e1);
	}
};