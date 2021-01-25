#pragma once
#include "../Utils/Utils.h"
#include <mutex>
using namespace std;

struct ElementAggregate;
struct AlgebraicElement
{
	BigNumber Number;
	vector<pair<AlgebraicElement*, double>> Neighbours;
	int NStrongNeighbours = 0;
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
	double _strongCouplingThreshold;

	const SparseMatrix* A;
public:
	vector<AlgebraicElement> _elements;
	vector<ElementAggregate> _coarseElements;

public:
	AlgebraicMesh(int blockSize, double strongCouplingThreshold)
	{
		this->_blockSize = blockSize;
		this->_strongCouplingThreshold = strongCouplingThreshold;
	}

	void Build(const SparseMatrix& A)
	{
		this->A = &A;

		BigNumber nElements = A.rows() / _blockSize;
		this->_elements = vector<AlgebraicElement>(nElements);

		// Numbering
		NumberParallelLoop<EmptyResultChunk> parallelLoop(_elements.size());
		parallelLoop.Execute([this](BigNumber elemNumber, ParallelChunk<EmptyResultChunk>* chunk)
			{
				AlgebraicElement& elem = _elements[elemNumber];
				elem.Number = elemNumber;
			});

		// Neighbours and coupling
		NumberParallelLoop<EmptyResultChunk> parallelLoop2(_elements.size());
		parallelLoop2.Execute([this, &A](BigNumber elemNumber, ParallelChunk<EmptyResultChunk>* chunk)
			{
				AlgebraicElement& elem = _elements[elemNumber];
				for (int k = 0; k < _blockSize; k++)
				{
					// RowMajor --> the following line iterates over the non-zeros of the elemNumber-th row.
					for (SparseMatrix::InnerIterator it(A, elemNumber*_blockSize + k); it; ++it)
					{
						BigNumber neighbourNumber = it.col() / _blockSize;
						if (neighbourNumber != elemNumber)
						{
							AlgebraicElement* neighbour = &_elements[neighbourNumber];
							if (find_if(elem.Neighbours.begin(), elem.Neighbours.end(), [neighbour](const pair<AlgebraicElement*, double>& n) { return n.first == neighbour; }) == elem.Neighbours.end())
							{
								double coupling = this->CouplingValue(elem, *neighbour);
								elem.Neighbours.push_back({ neighbour, coupling });
							}
						}
					}
				}

				sort(elem.Neighbours.begin(), elem.Neighbours.end(),
					[](const pair<AlgebraicElement*, double>& n1, const pair<AlgebraicElement*, double>& n2)
					{
						return n1.second < n2.second; // Sort by ascending coupling
					});

				for (auto it = elem.Neighbours.begin(); it != elem.Neighbours.end(); ++it)
				{
					double coupling = it->second;
					if (IsStronglyCoupled(elem, coupling))
						elem.NStrongNeighbours++;
					else
						break;
				}
			});
	}

	void PairWiseAggregate(bool& coarsestPossibleMeshReached)
	{
		_coarseElements.reserve(_elements.size() / 2);

		vector<AlgebraicElement*> sortedElements;
		for (AlgebraicElement& e : _elements)
			sortedElements.push_back(&e);
		sort(sortedElements.begin(), sortedElements.end(),
			[](AlgebraicElement* e1, AlgebraicElement* e2)
			{
				return e1->NStrongNeighbours < e2->NStrongNeighbours; // Sort by ascending number of strong neighbours
			});

		// Element aggregation
		for (BigNumber i = 0; i < sortedElements.size(); i++)
		{
			AlgebraicElement* elem = sortedElements[i];
			if (elem->IsAggregated)
				continue;

			if (elem->Neighbours.empty())
			{
				coarsestPossibleMeshReached = true;
				return;
			}

			AlgebraicElement* strongestNeighbour = StrongestNeighbour(*elem, true);
			// If no neighbour available, get the already aggregated strongest neighbour
			if (!strongestNeighbour)
				strongestNeighbour = StrongestNeighbour(*elem, false);
			
			if (!strongestNeighbour)
			{
				coarsestPossibleMeshReached = true;
				return;
			}

			Aggregate(*elem, *strongestNeighbour);
		}
	}

private:
	AlgebraicElement* StrongestNeighbour(const AlgebraicElement& e, bool checkAvailability)
	{
		AlgebraicElement* strongestNeighbour = nullptr;
		double strongestNegativeCoupling = 0;
		int smallestAggregateSize = 1000000;
		for (auto it = e.Neighbours.begin(); it != e.Neighbours.end(); ++it)
		{
			AlgebraicElement* n = it->first;
			double coupling = it->second;

			if (checkAvailability && n->IsAggregated)
				continue;

			if (!IsStronglyCoupled(e, coupling))
				break; // the neighbours are sorted, so if this one is not strongly coupled, the next ones won't be.

			if (checkAvailability)
			{
				strongestNeighbour = n;
				break;
			}
			else
			{
				if (coupling <= strongestNegativeCoupling)
				{
					strongestNegativeCoupling = coupling;
					int aggregateSize = n->CoarseElement->FineElements.size();
					if (aggregateSize < smallestAggregateSize)
					{
						strongestNeighbour = n;
						smallestAggregateSize = aggregateSize;
					}
				}
				else
					break;
			}
		}
		return strongestNeighbour;
	}

	double CouplingValue(const AlgebraicElement& e1, const AlgebraicElement& e2)
	{
		DenseMatrix couplingBlock = A->block(e1.Number*_blockSize, e2.Number*_blockSize, _blockSize, _blockSize);
		//double coupling = couplingBlock(0, 0);
		//Eigen::VectorXd eigenVals = (couplingBlock.transpose()*couplingBlock).selfadjointView().eigenvalues();
		//cout << "Trace = " << endl << couplingBlock.trace() << endl;
		/*Eigen::SelfAdjointEigenSolver<DenseMatrix> es(couplingBlock.transpose()*couplingBlock);
		Eigen::VectorXd eigenVals = es.eigenvalues();
		cout << "Eigenvalues = " << endl << eigenVals << endl;
		cout << "sqrt(Eigenvalues) = " << endl << eigenVals.cwiseSqrt() << endl;
		//cout << eigenVals.cwiseSqrt().minCoeff() << endl;
		cout << endl;
		double coupling = eigenVals.cwiseSqrt().minCoeff();*/
		double coupling = couplingBlock.trace();
		return coupling;
	}

	bool IsStronglyCoupled(const AlgebraicElement& e, double coupling)
	{
		return coupling < 0 && coupling < _strongCouplingThreshold * e.Neighbours[0].second;
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