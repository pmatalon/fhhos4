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
	vector<HybridAlgebraicElement*> Neighbours;
	bool IsAggregated = false;
	HybridElementAggregate* CoarseElement = nullptr;
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

	HybridElementAggregate(BigNumber number, vector<HybridAlgebraicElement*> elements)
		: Number(number), FineElements(elements) 
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

	const SparseMatrix* A_T_T;
	const SparseMatrix* A_T_F;
	const SparseMatrix* A_F_F;
public:
	vector<HybridAlgebraicElement> _elements;
	vector<HybridAlgebraicFace> _faces;
	vector<HybridElementAggregate> _coarseElements;
	vector<HybridFaceAggregate> _coarseFaces;

public:
	HybridAlgebraicMesh(int cellBlockSize, int faceBlockSize)
	{
		this->_cellBlockSize = cellBlockSize;
		this->_faceBlockSize = faceBlockSize;
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
						{
							//face->Number = faceNumber;
							elem.Faces.push_back(face);
							//face->Elements.push_back(&elem);
						}
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
							elem.Neighbours.push_back(neighbour);
					}
				}
			});
	}

	void PairWiseAggregate(CoarseningStrategy coarseningStgy, bool& coarsestPossibleMeshReached)
	{
		// Element aggregation
		_coarseElements.reserve(_elements.size() / 2);

		for (BigNumber i = 0; i < _elements.size(); i++)
		{
			HybridAlgebraicElement& elem = _elements[i];
			if (elem.IsAggregated)
				continue;

			if (elem.Neighbours.empty())
			{
				coarsestPossibleMeshReached = true;
				return;
			}

			HybridAlgebraicElement* strongestNeighbour = StrongestNeighbour(elem, true);
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

		// Face aggregation:

		this->_coarseFaces.reserve(_faces.size()); // must reserve sufficient space
		if (coarseningStgy == CoarseningStrategy::CAMGCollapseElementInterfaces ||
			coarseningStgy == CoarseningStrategy::CAMGCollapseElementInterfacesAndTyrAggregInteriorToBoundaries)
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


		if (coarseningStgy == CoarseningStrategy::CAMGCollapseElementInterfacesAndTyrAggregInteriorToBoundaries)
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
	HybridAlgebraicElement* StrongestNeighbour(const HybridAlgebraicElement& e, bool checkAvailability)
	{
		HybridAlgebraicElement* strongestNeighbour = nullptr;
		double strongestNegativeCoupling = 0;
		int smallestAggregateSize = 1000000;
		DenseMatrix elemBlock = A_T_T->block(e.Number*_cellBlockSize, e.Number*_cellBlockSize, _cellBlockSize, _cellBlockSize);
		double elemKappa = elemBlock(0, 0);
		for (HybridAlgebraicFace* f : e.Faces)
		{
			for (HybridAlgebraicElement* n : f->Elements)
			{
				if (n == &e || (checkAvailability && n->IsAggregated))
					continue;

				DenseMatrix neighourBlock = A_T_T->block(n->Number*_cellBlockSize, n->Number*_cellBlockSize, _cellBlockSize, _cellBlockSize);
				double neighbourKappa = neighourBlock(0, 0);
				double minKappa = min(elemKappa, neighbourKappa);
				double maxKappa = max(elemKappa, neighbourKappa);
				if (maxKappa/minKappa > 10)
					continue;

				DenseMatrix couplingElemFace = A_T_F->block(e.Number*_cellBlockSize, f->Number*_faceBlockSize, _cellBlockSize, _faceBlockSize);
				double coupling = couplingElemFace(0, 0);
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
		}
		return strongestNeighbour;
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

				//if (checkNeighbourCoarseFace && !n->CoarseFace)
					//continue;
				if (find(tabooList.begin(), tabooList.end(), n) != tabooList.end())
					continue;

				DenseMatrix couplingFaces = A_F_F->block(f.Number*_faceBlockSize, n->Number*_faceBlockSize, _faceBlockSize, _faceBlockSize);
				double coupling = couplingFaces(0, 0);
				if (/*!strongestNeighbour ||*/ coupling < strongestNegativeCoupling)
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

	void Aggregate(HybridAlgebraicElement& e1, HybridAlgebraicElement& e2)
	{
		assert(!(e1.IsAggregated && e2.IsAggregated));

		if (!e1.IsAggregated && !e2.IsAggregated)
		{
			_coarseElements.push_back(HybridElementAggregate(_coarseElements.size(), { &e1, &e2 }));
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