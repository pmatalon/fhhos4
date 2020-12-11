#pragma once
#include "Multigrid.h"
#include <mutex>
using namespace std;

struct AlgebraicFace;
struct ElementAggregate;
struct AlgebraicElement
{
	BigNumber Number;
	vector<AlgebraicFace*> Faces;
	vector<AlgebraicElement*> Neighbours;
	bool IsAggregated = false;
	ElementAggregate* CoarseElement;
};

struct ElementAggregate;
struct FaceAggregate;
struct AlgebraicFace
{
	BigNumber Number;
	vector<AlgebraicElement*> Elements;
	bool IsRemovedOnCoarseMesh = false;
	vector<ElementAggregate*> CoarseElements;
	FaceAggregate* CoarseFace;
	mutex Mutex;
};

struct ElementAggregate
{
	BigNumber Number;
	vector<AlgebraicElement*> FineElements;
	vector<AlgebraicFace*> FineFaces;
	map<ElementAggregate*, vector<AlgebraicFace*>> Neighbours;

	ElementAggregate(BigNumber number, vector<AlgebraicElement*> elements)
		: Number(number), FineElements(elements) 
	{}
};

struct FaceAggregate
{
	BigNumber Number;
	vector<AlgebraicFace*> FineFaces;
	
	FaceAggregate(BigNumber number, vector<AlgebraicFace*> faces)
		: Number(number), FineFaces(faces)
	{}
};

class CondensedAlgebraicMesh
{
private:
	int _cellBlockSize;
	int _faceBlockSize;

	const SparseMatrix* A_T_T;
	const SparseMatrix* A_T_F;
public:
	vector<AlgebraicElement> _elements;
	vector<AlgebraicFace> _faces;
	vector<ElementAggregate> _coarseElements;
	vector<FaceAggregate> _coarseFaces;

public:
	CondensedAlgebraicMesh(int cellBlockSize, int faceBlockSize)
	{
		this->_cellBlockSize = cellBlockSize;
		this->_faceBlockSize = faceBlockSize;
	}

	void Build(const SparseMatrix& A_T_T, const SparseMatrix& A_T_F)
	{
		assert(A_T_T.rows() > 0 && A_T_F.rows() > 0);
		assert(A_T_T.rows() == A_T_F.rows());

		this->A_T_T = &A_T_T;
		this->A_T_F = &A_T_F;

		if (!A_T_F.IsRowMajor)
			assert("A_T_F must be row-major");

		BigNumber nElements = A_T_T.rows() / _cellBlockSize;
		this->_elements = vector<AlgebraicElement>(nElements);

		BigNumber nFaces = A_T_F.cols() / _faceBlockSize;
		this->_faces = vector<AlgebraicFace>(nFaces);

		//-------------------------//
		// Filling elements' faces //
		//-------------------------//

		NumberParallelLoop<EmptyResultChunk> parallelLoopElem(_elements.size());
		parallelLoopElem.Execute([this, &A_T_F](BigNumber elemNumber, ParallelChunk<EmptyResultChunk>* chunk)
			{
				AlgebraicElement& elem = _elements[elemNumber];
				elem.Number = elemNumber;

				// RowMajor --> the following line iterates over the non-zeros of the elemNumber-th row.
				for (SparseMatrix::InnerIterator it(A_T_F, elemNumber*_cellBlockSize); it; ++it)
				{
					BigNumber faceNumber = it.col() / _faceBlockSize;
					AlgebraicFace* face = &_faces[faceNumber];
					if (find(elem.Faces.begin(), elem.Faces.end(), face) == elem.Faces.end())
					{
						//face->Number = faceNumber;
						elem.Faces.push_back(face);
						//face->Elements.push_back(&elem);
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
				AlgebraicFace& face = _faces[faceNumber];
				face.Number = faceNumber;

				// ColMajor --> the following line iterates over the non-zeros of the faceNumber-th col.
				for (ColMajorSparseMatrix::InnerIterator it(A_T_F_ColMajor, faceNumber*_faceBlockSize); it; ++it)
				{
					assert(it.col() / _faceBlockSize == faceNumber);
					BigNumber elemNumber = it.row() / _cellBlockSize;
					AlgebraicElement* elem = &_elements[elemNumber];
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
				AlgebraicElement& elem = _elements[elemNumber];
				for (AlgebraicFace* face : elem.Faces)
				{
					for (AlgebraicElement* neighbour : face->Elements)
					{
						if (neighbour->Number != elem.Number)
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
				coarsestPossibleMeshReached = true;

			AlgebraicElement* strongestNeighbour;
			for (AlgebraicElement* neighbour : elem.Neighbours)
			{
				//DenseMatrix blockI = A_T_F->block(i*_cellBlockSize, i*_cellBlockSize, _cellBlockSize, _cellBlockSize);
				if (!neighbour->IsAggregated)
				{
					Aggregate(elem, *neighbour);
					break;
				}
			}

			// If no neighbour available
			if (!elem.IsAggregated)
			{
				// Aggregate with first neighbour
				AlgebraicElement* neighbour = elem.Neighbours[0];
				Aggregate(elem, *neighbour);
			}
		}

		// Removal of faces shared by elements in the same aggregate.
		// Determination of the faces of the aggregates.
		NumberParallelLoop<EmptyResultChunk> parallelLoopCE(_coarseElements.size());
		parallelLoopCE.Execute([this](BigNumber coarseElemNumber)
			{
				ElementAggregate& coarseElem = _coarseElements[coarseElemNumber];
				for (int i = 0; i < coarseElem.FineElements.size(); i++)
				{
					AlgebraicElement* elem1 = coarseElem.FineElements[i];
					for (AlgebraicFace* face : elem1->Faces)
					{
						for (int j = i + 1; j < coarseElem.FineElements.size(); j++)
						{
							AlgebraicElement* elem2 = coarseElem.FineElements[j];
							if (find(elem2->Faces.begin(), elem2->Faces.end(), face) != elem2->Faces.end())
							{
								// This face is shared by elem1 and elem2, so we remove it on the coarse grid
								face->IsRemovedOnCoarseMesh = true;
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
				ElementAggregate* coarseElem = &_coarseElements[coarseElemNumber];
				set<ElementAggregate*> neighbours;
				for (AlgebraicFace* face : coarseElem->FineFaces)
				{
					for (ElementAggregate* neighbour : face->CoarseElements)
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
		// Collapse faces interfacing two element aggregates.
		_coarseFaces.reserve(_faces.size() / 3);
		for (ElementAggregate& coarseElem : _coarseElements)
		{
			for (auto it = coarseElem.Neighbours.begin(); it != coarseElem.Neighbours.end(); it++)
			{
				ElementAggregate* neighbour = it->first;
				vector<AlgebraicFace*> fineFaces = it->second;
				if (!fineFaces.front()->CoarseFace)
				{
					_coarseFaces.emplace_back(_coarseFaces.size(), fineFaces);
					for (AlgebraicFace* fineFace : fineFaces)
						fineFace->CoarseFace = &_coarseFaces.back();
				}
			}
		}
	}

private:
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