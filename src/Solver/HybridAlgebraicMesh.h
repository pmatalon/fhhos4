#pragma once
#include "../Utils/Utils.h"
#include <mutex>
using namespace std;

struct HybridAlgebraicFace;
struct HybridElementAggregate;
struct HybridAlgebraicElement
{
	BigNumber Number;
	vector<HybridAlgebraicFace*> Faces;
	vector<pair<HybridAlgebraicElement*, double>> Neighbours;
	vector<HybridAlgebraicElement*> StrongNeighbours; // sorted by descending strength (the strongest neighbour is first)
	HybridElementAggregate* CoarseElement = nullptr;
	int NElementsIAmStrongNeighbourOf = 0;
};

struct HybridElementAggregate;
struct HybridFaceAggregate;
struct HybridAlgebraicFace
{
	BigNumber Number;
	vector<HybridAlgebraicElement*> Elements;
	bool IsRemovedOnCoarseMesh = false;
	vector<HybridElementAggregate*> CoarseElements;
	HybridFaceAggregate* CoarseFace = nullptr;
	mutex Mutex;
};

struct HybridElementAggregate
{
	BigNumber Number;
	vector<HybridAlgebraicElement*> FineElements;
	vector<HybridAlgebraicFace*> FineFaces;
	vector<HybridAlgebraicFace*> RemovedFineFaces;
	map<HybridElementAggregate*, vector<HybridAlgebraicFace*>> Neighbours;
	vector<HybridFaceAggregate*> CoarseFaces;

	HybridElementAggregate(BigNumber number)
		: Number(number)
	{}
};

struct HybridFaceAggregate
{
	BigNumber Number;
	vector<HybridAlgebraicFace*> FineFaces;
	
	HybridFaceAggregate(BigNumber number, vector<HybridAlgebraicFace*> faces)
		: Number(number), FineFaces(faces)
	{}
};

class HybridAlgebraicMesh
{
private:
	int _cellBlockSize;
	int _faceBlockSize;
	double _strongCouplingThreshold;

	const SparseMatrix* A_T_T;
	const SparseMatrix* A_T_F;
	const SparseMatrix* A_F_F;
public:
	vector<HybridAlgebraicElement> _elements;
	vector<HybridAlgebraicFace> _faces;
	vector<HybridElementAggregate> _coarseElements;
	vector<HybridFaceAggregate> _coarseFaces;

public:
	HybridAlgebraicMesh(int cellBlockSize, int faceBlockSize, double strongCouplingThreshold)
	{
		this->_cellBlockSize = cellBlockSize;
		this->_faceBlockSize = faceBlockSize;
		this->_strongCouplingThreshold = strongCouplingThreshold;
	}

	void Build(const SparseMatrix& A_T_T, const SparseMatrix& A_T_F, const SparseMatrix& A_F_F)
	{
		assert(A_T_T.rows() > 0 && A_T_F.rows() > 0);
		assert(A_T_T.rows() == A_T_F.rows());

		this->A_T_T = &A_T_T;
		this->A_T_F = &A_T_F;
		this->A_F_F = &A_F_F;

		if (!A_T_F.IsRowMajor)
			assert("A_T_F must be row-major");

		BigNumber nElements = A_T_F.rows() / _cellBlockSize;
		this->_elements = vector<HybridAlgebraicElement>(nElements);

		BigNumber nFaces = A_T_F.cols() / _faceBlockSize;
		this->_faces = vector<HybridAlgebraicFace>(nFaces);

		//-------------------------//
		// Filling elements' faces //
		//-------------------------//

		NumberParallelLoop<EmptyResultChunk> parallelLoopElem(_elements.size());
		parallelLoopElem.Execute([this, &A_T_F](BigNumber elemNumber, ParallelChunk<EmptyResultChunk>* chunk)
			{
				HybridAlgebraicElement& elem = _elements[elemNumber];
				elem.Number = elemNumber;

				for (int k = 0; k < _cellBlockSize; k++)
				{
					// RowMajor --> the following line iterates over the non-zeros of the elemNumber-th row.
					for (SparseMatrix::InnerIterator it(A_T_F, elemNumber*_cellBlockSize + k); it; ++it)
					{
						BigNumber faceNumber = it.col() / _faceBlockSize;
						HybridAlgebraicFace* face = &_faces[faceNumber];
						if (find(elem.Faces.begin(), elem.Faces.end(), face) == elem.Faces.end())
							elem.Faces.push_back(face);
					}
				}
			});

		//-------------------------//
		// Filling faces' elements //
		//-------------------------//

		ColMajorSparseMatrix A_T_F_ColMajor = A_T_F;

		NumberParallelLoop<EmptyResultChunk> parallelLoopFace(_faces.size());
		parallelLoopFace.Execute([this, &A_T_F_ColMajor](BigNumber faceNumber)
			{
				HybridAlgebraicFace& face = _faces[faceNumber];
				face.Number = faceNumber;

				// ColMajor --> the following line iterates over the non-zeros of the faceNumber-th col.
				for (ColMajorSparseMatrix::InnerIterator it(A_T_F_ColMajor, faceNumber*_faceBlockSize); it; ++it)
				{
					assert(it.col() / _faceBlockSize == faceNumber);
					BigNumber elemNumber = it.row() / _cellBlockSize;
					HybridAlgebraicElement* elem = &_elements[elemNumber];
					if (find(face.Elements.begin(), face.Elements.end(), elem) == face.Elements.end())
						face.Elements.push_back(elem);
				}
			});

		//--------------------//
		// Element neighbours //
		//--------------------//

		parallelLoopElem = NumberParallelLoop<EmptyResultChunk>(_elements.size());
		parallelLoopElem.Execute([this](BigNumber elemNumber)
			{
				HybridAlgebraicElement& elem = _elements[elemNumber];
				for (HybridAlgebraicFace* face : elem.Faces)
				{
					for (HybridAlgebraicElement* neighbour : face->Elements)
					{
						if (neighbour->Number != elem.Number)
						{
							double coupling = this->CouplingValue(elem, *face);
							elem.Neighbours.push_back({ neighbour, coupling });
						}
					}
				}

				sort(elem.Neighbours.begin(), elem.Neighbours.end(),
					[](const pair<HybridAlgebraicElement*, double>& n1, const pair<HybridAlgebraicElement*, double>& n2)
					{
						return n1.second < n2.second; // Sort by ascending coupling
					});

				double elemKappa = this->DiffusionCoeff(elem);
				for (auto it = elem.Neighbours.begin(); it != elem.Neighbours.end(); ++it)
				{
					HybridAlgebraicElement* neighbour = it->first;
					double coupling = it->second;
					if (IsStronglyCoupled(elem, coupling))
					{
						double neighbourKappa = this->DiffusionCoeff(*neighbour);
						double minKappa = min(elemKappa, neighbourKappa);
						double maxKappa = max(elemKappa, neighbourKappa);
						if (maxKappa / minKappa < 10)
						{
							elem.StrongNeighbours.push_back(neighbour);
							neighbour->NElementsIAmStrongNeighbourOf++;
						}
					}
					else
						break;
				}
			});
	}

	void PairWiseAggregate(CoarseningStrategy coarseningStgy, bool& coarsestPossibleMeshReached)
	{
		PairWiseElementAggregate(coarsestPossibleMeshReached);
		
		// Removal of faces shared by elements in the same aggregate.
		// Determination of the faces of the aggregates.
		NumberParallelLoop<EmptyResultChunk> parallelLoopCE(_coarseElements.size());
		parallelLoopCE.Execute([this](BigNumber coarseElemNumber)
			{
				HybridElementAggregate& coarseElem = _coarseElements[coarseElemNumber];
				for (int i = 0; i < coarseElem.FineElements.size(); i++)
				{
					HybridAlgebraicElement* elem1 = coarseElem.FineElements[i];
					for (HybridAlgebraicFace* face : elem1->Faces)
					{
						for (int j = i + 1; j < coarseElem.FineElements.size(); j++)
						{
							HybridAlgebraicElement* elem2 = coarseElem.FineElements[j];
							if (find(elem2->Faces.begin(), elem2->Faces.end(), face) != elem2->Faces.end())
							{
								// This face is shared by elem1 and elem2, so we remove it on the coarse grid
								face->IsRemovedOnCoarseMesh = true;
								coarseElem.RemovedFineFaces.push_back(face);
								//face->CoarseElements.push_back(&coarseElem);
								break;
							}
						}

						if (!face->IsRemovedOnCoarseMesh)
						{
							// If it is not an inner face, then it's a face of the aggregate
							face->Mutex.lock();
							face->CoarseElements.push_back(&coarseElem);
							face->Mutex.unlock();
							coarseElem.FineFaces.push_back(face);
						}
					}
				}
			});

		// Computation of the aggregates' neighbours
		parallelLoopCE = NumberParallelLoop<EmptyResultChunk>(_coarseElements.size());
		parallelLoopCE.Execute([this](BigNumber coarseElemNumber)
			{
				HybridElementAggregate* coarseElem = &_coarseElements[coarseElemNumber];
				set<HybridElementAggregate*> neighbours;
				for (HybridAlgebraicFace* face : coarseElem->FineFaces)
				{
					for (HybridElementAggregate* neighbour : face->CoarseElements)
					{
						if (neighbour != coarseElem)
						{
							auto it = coarseElem->Neighbours.find(neighbour);
							if (it == coarseElem->Neighbours.end())
								coarseElem->Neighbours.insert({ neighbour, {face} });
							else
								it->second.push_back(face);
						}
					}
				}
			});

		// Face aggregation by collapsing multiple interfaces
		this->_coarseFaces.reserve(_faces.size()); // must reserve sufficient space
		if (coarseningStgy == CoarseningStrategy::CAMGCollapseElementInterfaces ||
			coarseningStgy == CoarseningStrategy::CAMGCollapseElementInterfacesAndTryAggregInteriorToBoundaries)
		{
			// Collapse faces interfacing two element aggregates.
			for (HybridElementAggregate& coarseElem : _coarseElements)
			{
				for (auto it = coarseElem.Neighbours.begin(); it != coarseElem.Neighbours.end(); it++)
				{
					HybridElementAggregate* neighbour = it->first;
					vector<HybridAlgebraicFace*> fineFaces = it->second;
					if (!fineFaces.front()->CoarseFace)
					{
						_coarseFaces.emplace_back(_coarseFaces.size(), fineFaces);
						HybridFaceAggregate* coarseFace = &_coarseFaces.back(); // if _coarseFaces is reallocated, all the pointers already taken are invalid
						for (HybridAlgebraicFace* fineFace : fineFaces)
							fineFace->CoarseFace = coarseFace;
						coarseElem.CoarseFaces.push_back(coarseFace);
						neighbour->CoarseFaces.push_back(coarseFace);
					}
				}
			}
		}

		if (coarseningStgy == CoarseningStrategy::CAMGCollapseElementInterfacesAndTryAggregInteriorToBoundaries)
		{
			// Agglomerate removed faces
			for (HybridAlgebraicFace& face : _faces)
			{
				if (!face.IsRemovedOnCoarseMesh)
					continue;

				assert(!face.CoarseFace);

				HybridAlgebraicFace* strongestNeighbour = StrongestNeighbour(face);
				if (strongestNeighbour)
				{
					face.CoarseFace = strongestNeighbour->CoarseFace;
					strongestNeighbour->CoarseFace->FineFaces.push_back(&face);
				}
				else
				{
					//cout << "Face " << face.Number << " does not have any strong neighbour" << endl;
					/*for (SparseMatrix::InnerIterator it(*A_F_F, face.Number); it; ++it)
						cout << "(" << it.row() << ", " << it.col() << ") " << it.value() << endl;
					coarsestPossibleMeshReached = true;
					return;*/
				}
			}
		}
	}

private:
	void PairWiseElementAggregate(bool& coarsestPossibleMeshReached)
	{
		_coarseElements.reserve(_elements.size()); // we need to be sure that the vector won't be resized

		HybridAlgebraicElement* nextElement = &_elements[0];

		// Element aggregation
		while (nextElement)
		{
			HybridAlgebraicElement* elem = nextElement;
			nextElement = nullptr;
			int currentMin = elem->NElementsIAmStrongNeighbourOf;

			assert(!elem->CoarseElement);

			// New aggregate
			_coarseElements.emplace_back(_coarseElements.size());
			HybridElementAggregate* aggregate = &_coarseElements.back();

			// Add elem
			AddToAggregate(*elem, *aggregate);

			HybridAlgebraicElement* neighbour = StrongestAvailableNeighbour(*elem);
			if (neighbour)
			{
				// Add neighbour?
				bool strongConnectionIsReciprocal = find(neighbour->StrongNeighbours.begin(), neighbour->StrongNeighbours.end(), elem) != neighbour->StrongNeighbours.end();
				if (strongConnectionIsReciprocal)
					AddToAggregate(*neighbour, *aggregate);
			}

			nextElement = NextElementInTheNeighbourhood(*aggregate);

			if (!nextElement)
			{
				auto it = find_if(_elements.begin(), _elements.end(), [](const HybridAlgebraicElement& e) { return !e.CoarseElement; });
				if (it != _elements.end())
					nextElement = &*it;
			}
		}
	}

	static bool CompareNElementsIAmStrongNeighbourOf(HybridAlgebraicElement* e1, HybridAlgebraicElement* e2)
	{
		return e1->NElementsIAmStrongNeighbourOf < e2->NElementsIAmStrongNeighbourOf; // Sort by ascending number
	}

	HybridAlgebraicElement* StrongestAvailableNeighbour(const HybridAlgebraicElement& e)
	{
		auto it = find_if(e.StrongNeighbours.begin(), e.StrongNeighbours.end(), [](HybridAlgebraicElement* n) { return !n->CoarseElement; });
		if (it != e.StrongNeighbours.end())
			return *it;
		return nullptr;
	}

	double CouplingValue(const HybridAlgebraicElement& e, const HybridAlgebraicFace& f)
	{
		DenseMatrix couplingBlock = A_T_F->block(e.Number*_cellBlockSize, f.Number*_faceBlockSize, _cellBlockSize, _faceBlockSize);
		double coupling = couplingBlock.trace();
		return coupling;
	}

	bool IsStronglyCoupled(const HybridAlgebraicElement& e, double coupling)
	{
		return coupling < 0 && coupling < _strongCouplingThreshold * e.Neighbours[0].second;
	}

	double DiffusionCoeff(const HybridAlgebraicElement& e)
	{
		DenseMatrix elemBlock = A_T_T->block(e.Number*_cellBlockSize, e.Number*_cellBlockSize, _cellBlockSize, _cellBlockSize);
		double kappa = elemBlock(0, 0);
		return kappa;
	}

	HybridAlgebraicFace* StrongestNeighbour(const HybridAlgebraicFace& f, vector<const HybridAlgebraicFace*> tabooList = {})
	{
		HybridAlgebraicFace* strongestNeighbour = nullptr;
		double strongestNegativeCoupling = 0;
		int smallestAggregateSize = 1000000;
		for (HybridAlgebraicElement* e : f.Elements)
		{
			for (HybridAlgebraicFace* n : e->Faces)
			{
				if (n == &f)
					continue;

				if (find(tabooList.begin(), tabooList.end(), n) != tabooList.end())
					continue;

				double coupling = this->CouplingValue(f, *n);
				if (coupling < strongestNegativeCoupling)
				{
					if (!n->CoarseFace)
					{
						tabooList.push_back(&f);
						n = StrongestNeighbour(*n, tabooList);
						if (!n)
							continue;
						assert(n->CoarseFace);
					}

					assert(n->CoarseFace);
					int aggregateSize = n->CoarseFace->FineFaces.size();
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

	double CouplingValue(const HybridAlgebraicFace& f1, const HybridAlgebraicFace& f2)
	{
		DenseMatrix couplingBlock = A_F_F->block(f1.Number*_faceBlockSize, f2.Number*_faceBlockSize, _faceBlockSize, _faceBlockSize);
		double coupling = couplingBlock.trace();
		return coupling;
	}

	void AddToAggregate(HybridAlgebraicElement& e, HybridElementAggregate& aggregate)
	{
		aggregate.FineElements.push_back(&e);
		e.CoarseElement = &aggregate;

		for (HybridAlgebraicElement* n : e.StrongNeighbours)
		{
			if (!n->CoarseElement)
				n->NElementsIAmStrongNeighbourOf--;
		}
	}

	HybridAlgebraicElement* NextElementInTheNeighbourhood(const HybridElementAggregate& aggregate)
	{
		HybridAlgebraicElement* nextElement = nullptr;
		int min = 1000000000;
		for (HybridAlgebraicElement* e : aggregate.FineElements)
		{
			for (HybridAlgebraicElement* n : e->StrongNeighbours)
			{
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