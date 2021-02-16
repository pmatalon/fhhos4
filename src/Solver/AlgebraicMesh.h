#pragma once
#include "../Utils/Utils.h"
#include "PairwiseAggregation.h"
#include "AllNeighbourAggregation.h"
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
		PairwiseAggregation<AlgebraicElement, ElementAggregate> aggregProcess;
		_coarseElements = aggregProcess.Perform(_elements, coarsestPossibleMeshReached);
	}

	void AllNeighbourAggregate(bool& coarsestPossibleMeshReached)
	{
		AllNeighbourAggregation<AlgebraicElement, ElementAggregate> aggregProcess;
		_coarseElements = aggregProcess.Perform(_elements, coarsestPossibleMeshReached);
	}

private:
	static bool CompareNElementsIAmStrongNeighbourOf(AlgebraicElement* e1, AlgebraicElement* e2)
	{
		return e1->NElementsIAmStrongNeighbourOf < e2->NElementsIAmStrongNeighbourOf; // Sort by ascending number
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
};