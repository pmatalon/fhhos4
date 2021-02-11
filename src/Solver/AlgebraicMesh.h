#pragma once
#include "../Utils/Utils.h"
#include <mutex>
using namespace std;

struct ElementAggregate;
struct AlgebraicElement
{
	BigNumber Number;
	mutex Mutex;
	vector<pair<AlgebraicElement*, double>> Neighbours;
	vector<AlgebraicElement*> StrongNeighbours; // sorted by descending strength (the strongest neighbour is first)
	int NElementsIAmStrongNeighbourOf = 0;
	ElementAggregate* CoarseElement = nullptr;
	ElementAggregate* FinalAggregate = nullptr;
	BigNumber FinalAggregateNumber;
};

struct ElementAggregate
{
	BigNumber Number;
	vector<AlgebraicElement*> FineElements;

	ElementAggregate(BigNumber number)
		: Number(number)
	{}

	ElementAggregate()
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

	AlgebraicMesh* FinerMesh = nullptr;

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
					AlgebraicElement* neighbour = it->first;
					double coupling = it->second;
					if (IsStronglyCoupled(elem, coupling))
					{
						elem.StrongNeighbours.push_back(neighbour);
						neighbour->Mutex.lock();
						neighbour->NElementsIAmStrongNeighbourOf++;
						neighbour->Mutex.unlock();
					}
					else
						break;
				}
			});
	}

	void PairWiseAggregate(bool& coarsestPossibleMeshReached)
	{
		_coarseElements.reserve(_elements.size()); // we need to be sure that the vector won't be resized

		int nSingletons = 0;
		int nPairs = 0;

		BigNumber remainingElements = _elements.size();
		AlgebraicElement* nextElement = &_elements[0];
		int min = _elements[0].NElementsIAmStrongNeighbourOf;
		for (AlgebraicElement& e : _elements)
		{
			if (e.NElementsIAmStrongNeighbourOf < min)
			{
				nextElement = &e;
				min = e.NElementsIAmStrongNeighbourOf;
			}
		}

		// Element aggregation
		while (nextElement)
		{
			AlgebraicElement* elem = nextElement;
			nextElement = nullptr;
			int currentMin = elem->NElementsIAmStrongNeighbourOf;

			assert(!elem->CoarseElement);

			// New aggregate
			_coarseElements.emplace_back(_coarseElements.size());
			ElementAggregate* aggregate = &_coarseElements.back();

			// Add elem
			AddToAggregate(*elem, *aggregate);
			remainingElements--;

			AlgebraicElement* neighbour = StrongestAvailableNeighbour(*elem);
			if (neighbour)
			{
				// Add neighbour?
				//bool strongConnectionIsReciprocal = find(neighbour->StrongNeighbours.begin(), neighbour->StrongNeighbours.end(), elem) != neighbour->StrongNeighbours.end();
				//if (strongConnectionIsReciprocal)
				//{
					AddToAggregate(*neighbour, *aggregate);
					remainingElements--;
					nPairs++;
				//}
				//else
					//nSingletons++;
			}
			else
				nSingletons++;

			nextElement = NextElementInTheNeighbourhood(*aggregate);

			if (!nextElement && remainingElements > 0)
			{
				// Browse the previous aggregates in the reverse order to find the next element
				BigNumber i = aggregate->Number - 1;
				while (!nextElement && i >= 0)
				//for (BigNumber i = aggregate->Number - 1; i >= 0; i--)
				{
					const ElementAggregate& previousAgg = _coarseElements[i];
					nextElement = NextElementInTheNeighbourhood(previousAgg);
					i--;
				}

				if (!nextElement)
				{
					auto it = find_if(_elements.begin(), _elements.end(), [](const AlgebraicElement& e) { return !e.CoarseElement; });
					if (it != _elements.end())
						nextElement = &*it;
				}
				/*min = 100000;
				for (AlgebraicElement& e : _elements)
				{
					if (!e.CoarseElement && e.NElementsIAmStrongNeighbourOf < min)
					{
						nextElement = &e;
						min = e.NElementsIAmStrongNeighbourOf;
					}
				}*/
			}
		}

		//cout << "Singletons: " << nSingletons << ", pairs: " << nPairs << endl;
	}

private:
	static bool CompareNElementsIAmStrongNeighbourOf(AlgebraicElement* e1, AlgebraicElement* e2)
	{
		return e1->NElementsIAmStrongNeighbourOf < e2->NElementsIAmStrongNeighbourOf; // Sort by ascending number
	}

	AlgebraicElement* StrongestAvailableNeighbour(const AlgebraicElement& e)
	{
		auto it = find_if(e.StrongNeighbours.begin(), e.StrongNeighbours.end(), [](AlgebraicElement* n) { return !n->CoarseElement; });
		if (it != e.StrongNeighbours.end())
			return *it;
		return nullptr;
		/*double bestCoupling = e.Neighbours[0].second;

		for (int i = 0; i< e.Neighbours.size(); i++)
		{
			const pair<AlgebraicElement*, double>& p = e.Neighbours[i];
			AlgebraicElement* n = p.first;
			if (n->CoarseElement)
				continue;

			if (i == 0)
				return n;

			double coupling = p.second;
			if (abs(coupling - bestCoupling) < Utils::Eps * bestCoupling)
				return n;
			else
				return nullptr;
		}
		return nullptr;*/
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

	void AddToAggregate(AlgebraicElement& e, ElementAggregate& aggregate)
	{
		aggregate.FineElements.push_back(&e);
		e.CoarseElement = &aggregate;

		for (AlgebraicElement* n : e.StrongNeighbours)
		{
			if (!n->CoarseElement)
				n->NElementsIAmStrongNeighbourOf--;
		}
	}

	AlgebraicElement* NextElementInTheNeighbourhood(const ElementAggregate& aggregate)
	{
		AlgebraicElement* nextElement = nullptr;
		int min = 1000000000;
		for (AlgebraicElement* e : aggregate.FineElements)
		{
			for (pair<AlgebraicElement*, double>& p : e->Neighbours)
			{
				AlgebraicElement* n = p.first;
				if (!n->CoarseElement && n->NElementsIAmStrongNeighbourOf < min)
				{
					nextElement = n;
					min = n->NElementsIAmStrongNeighbourOf;
				}
			}
		}
		return nextElement;
	}
};