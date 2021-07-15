#pragma once
#include "../../../Utils/Utils.h"
#include "../AggregAMG/AlgebraicMesh.h"
using namespace std;

struct HybridAlgebraicFace;
struct HybridElementAggregate;
struct HybridAlgebraicElement
{
	BigNumber Number;
	mutex Mutex;
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

	int LocalFineElementNumber(HybridAlgebraicElement* fe) const
	{
		for (int i = 0; i < FineElements.size(); ++i)
		{
			if (fe == FineElements[i])
				return i;
		}
		assert(false);
	}

	int LocalFineFaceNumber(HybridAlgebraicFace* ff) const
	{
		for (int i = 0; i < FineFaces.size(); ++i)
		{
			if (ff == FineFaces[i])
				return i;
		}
		assert(false);
	}

	int LocalRemovedFineFaceNumber(HybridAlgebraicFace* ff) const
	{
		for (int i = 0; i < RemovedFineFaces.size(); ++i)
		{
			if (ff == RemovedFineFaces[i])
				return i;
		}
		assert(false);
	}

	int LocalCoarseFaceNumber(HybridFaceAggregate* cf) const
	{
		for (int i = 0; i < CoarseFaces.size(); ++i)
		{
			if (cf == CoarseFaces[i])
				return i;
		}
		assert(false);
	}
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

public:
	const SparseMatrix* A_T_T;
	const SparseMatrix* A_T_F;
	const SparseMatrix* A_F_F;

	vector<HybridAlgebraicElement> Elements;
	vector<HybridAlgebraicFace> Faces;
	vector<HybridElementAggregate> CoarseElements;
	vector<HybridFaceAggregate> CoarseFaces;

public:
	HybridAlgebraicMesh(const SparseMatrix* A_T_T, const SparseMatrix* A_T_F, const SparseMatrix* A_F_F, int cellBlockSize, int faceBlockSize, double strongCouplingThreshold)
	{
		assert(A_T_T->rows() > 0 && A_T_F->rows() > 0);
		assert(A_T_T->rows() == A_T_F->rows());

		this->A_T_T = A_T_T;
		this->A_T_F = A_T_F;
		this->A_F_F = A_F_F;

		if (!A_T_F->IsRowMajor)
			assert("A_T_F must be row-major");

		this->_cellBlockSize = cellBlockSize;
		this->_faceBlockSize = faceBlockSize;
		this->_strongCouplingThreshold = strongCouplingThreshold;
	}

	void Build()
	{
		const SparseMatrix& A_T_T = *this->A_T_T;
		const SparseMatrix& A_T_F = *this->A_T_F;

		BigNumber nElements = A_T_F.rows() / _cellBlockSize;
		this->Elements = vector<HybridAlgebraicElement>(nElements);

		BigNumber nFaces = A_T_F.cols() / _faceBlockSize;
		this->Faces = vector<HybridAlgebraicFace>(nFaces);

		//cout << nElements << " elements, " << nFaces << " faces found" << endl;

		//-------------------------//
		// Filling elements' faces //
		//-------------------------//

		NumberParallelLoop<EmptyResultChunk> parallelLoopElem(Elements.size());
		parallelLoopElem.Execute([this, &A_T_F](BigNumber elemNumber, ParallelChunk<EmptyResultChunk>* chunk)
			{
				HybridAlgebraicElement& elem = Elements[elemNumber];
				elem.Number = elemNumber;

				for (int k = 0; k < _cellBlockSize; k++)
				{
					// RowMajor --> the following line iterates over the non-zeros of the elemNumber-th row.
					for (SparseMatrix::InnerIterator it(A_T_F, elemNumber*_cellBlockSize + k); it; ++it)
					{
						BigNumber faceNumber = it.col() / _faceBlockSize;
						HybridAlgebraicFace* face = &Faces[faceNumber];
						if (find(elem.Faces.begin(), elem.Faces.end(), face) == elem.Faces.end())
							elem.Faces.push_back(face);
					}
				}

				if (elem.Faces.empty())
					Utils::Error("Element " + to_string(elemNumber) + " has no face (no non-zero coefficient in row " + to_string(elemNumber) + " of A_TF)");
			});

		//-------------------------//
		// Filling faces' elements //
		//-------------------------//

		ColMajorSparseMatrix A_T_F_ColMajor = A_T_F;

		NumberParallelLoop<EmptyResultChunk> parallelLoopFace(Faces.size());
		parallelLoopFace.Execute([this, &A_T_F_ColMajor](BigNumber faceNumber)
			{
				HybridAlgebraicFace& face = Faces[faceNumber];
				face.Number = faceNumber;

				// ColMajor --> the following line iterates over the non-zeros of the faceNumber-th col.
				for (ColMajorSparseMatrix::InnerIterator it(A_T_F_ColMajor, faceNumber*_faceBlockSize); it; ++it)
				{
					assert(it.col() / _faceBlockSize == faceNumber);
					BigNumber elemNumber = it.row() / _cellBlockSize;
					HybridAlgebraicElement* elem = &Elements[elemNumber];
					if (find(face.Elements.begin(), face.Elements.end(), elem) == face.Elements.end())
						face.Elements.push_back(elem);
				}
			});

		//--------------------//
		// Element neighbours //
		//--------------------//

		parallelLoopElem = NumberParallelLoop<EmptyResultChunk>(Elements.size());
		parallelLoopElem.Execute([this](BigNumber elemNumber)
			{
				HybridAlgebraicElement& elem = Elements[elemNumber];
				for (HybridAlgebraicFace* face : elem.Faces)
				{
					for (HybridAlgebraicElement* neighbour : face->Elements)
					{
						if (neighbour->Number != elem.Number)
						{
							double coupling = this->CouplingValue(elem, *neighbour, *face);
							elem.Neighbours.push_back({ neighbour, coupling });
						}
					}
				}

				sort(elem.Neighbours.begin(), elem.Neighbours.end(),
					[](const pair<HybridAlgebraicElement*, double>& n1, const pair<HybridAlgebraicElement*, double>& n2)
					{
						return n1.second < n2.second; // Sort by ascending coupling
					});

				for (auto it = elem.Neighbours.begin(); it != elem.Neighbours.end(); ++it)
				{
					HybridAlgebraicElement* neighbour = it->first;
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

	void Coarsen(H_CoarsStgy elemCoarseningStgy, FaceCoarseningStrategy faceCoarseningStgy, bool& coarsestPossibleMeshReached)
	{
		coarsestPossibleMeshReached = false;

		if (elemCoarseningStgy == H_CoarsStgy::DoublePairwiseAggregation || elemCoarseningStgy == H_CoarsStgy::MultiplePairwiseAggregation)
		{
			cout << "\tElement pairwise aggregation" << endl;
			PairwiseAggregation<HybridAlgebraicElement, HybridElementAggregate> aggregProcess;
			CoarseElements = aggregProcess.Perform(Elements, coarsestPossibleMeshReached);
		}
		else if (elemCoarseningStgy == H_CoarsStgy::AgglomerationCoarseningByFaceNeighbours || elemCoarseningStgy == H_CoarsStgy::MultipleAgglomerationCoarseningByFaceNeighbours)
		{
			cout << "\tElement agglomeration" << endl;
			AllNeighbourAggregation<HybridAlgebraicElement, HybridElementAggregate> aggregProcess;
			CoarseElements = aggregProcess.Perform(Elements, coarsestPossibleMeshReached);
		}
		else
			Utils::FatalError("Unmanaged element coarsening strategy");

		if (coarsestPossibleMeshReached)
			return;
		
		// Removal of faces shared by elements in the same aggregate.
		// Determination of the faces of the aggregates.
		NumberParallelLoop<EmptyResultChunk> parallelLoopCE(CoarseElements.size());
		parallelLoopCE.Execute([this](BigNumber coarseElemNumber)
			{
				HybridElementAggregate& coarseElem = CoarseElements[coarseElemNumber];
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
								face->CoarseElements.push_back(&coarseElem);
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
		parallelLoopCE = NumberParallelLoop<EmptyResultChunk>(CoarseElements.size());
		parallelLoopCE.Execute([this](BigNumber coarseElemNumber)
			{
				HybridElementAggregate* coarseElem = &CoarseElements[coarseElemNumber];
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
		this->CoarseFaces.reserve(Faces.size()); // must reserve sufficient space
		if (faceCoarseningStgy == FaceCoarseningStrategy::InterfaceCollapsing ||
			faceCoarseningStgy == FaceCoarseningStrategy::InterfaceCollapsingAndTryAggregInteriorToInterfaces)
		{
			cout << "\tInterface collapsing" << endl;

			// Collapse faces interfacing two element aggregates.
			for (HybridElementAggregate& coarseElem : CoarseElements)
			{
				for (auto it = coarseElem.Neighbours.begin(); it != coarseElem.Neighbours.end(); it++)
				{
					HybridElementAggregate* neighbour = it->first;
					vector<HybridAlgebraicFace*> fineFaces = it->second;
					if (!fineFaces.front()->CoarseFace)
					{
						CoarseFaces.emplace_back(CoarseFaces.size(), fineFaces);
						HybridFaceAggregate* coarseFace = &CoarseFaces.back(); // if CoarseFaces is reallocated, all the pointers already taken are invalid
						for (HybridAlgebraicFace* fineFace : fineFaces)
							fineFace->CoarseFace = coarseFace;
						coarseElem.CoarseFaces.push_back(coarseFace);
						neighbour->CoarseFaces.push_back(coarseFace);
					}
				}
			}
		}

		if (faceCoarseningStgy == FaceCoarseningStrategy::InterfaceCollapsingAndTryAggregInteriorToInterfaces)
		{
			cout << "\tInterior faces agglomeration" << endl;

			// Agglomerate removed faces
			for (HybridAlgebraicFace& face : Faces)
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

		//cout << CoarseElements.size() << " coarse elements, " << CoarseFaces.size() << " coarse faces" << endl;
		if (CoarseFaces.empty())
			coarsestPossibleMeshReached = true;
	}

private:
	static bool CompareNElementsIAmStrongNeighbourOf(HybridAlgebraicElement* e1, HybridAlgebraicElement* e2)
	{
		return e1->NElementsIAmStrongNeighbourOf < e2->NElementsIAmStrongNeighbourOf; // Sort by ascending number
	}

	double CouplingValue(const HybridAlgebraicElement& e1, const HybridAlgebraicElement& e2, const HybridAlgebraicFace& f)
	{
		DenseMatrix couplingBlock = A_T_F->block(e1.Number*_cellBlockSize, f.Number*_faceBlockSize, _cellBlockSize, _faceBlockSize);
		double coupling = couplingBlock.trace();

		//double e1Kappa = this->DiffusionCoeff(e1);
		//double e2Kappa = this->DiffusionCoeff(e2);
		double e1Kappa = abs(coupling);
		DenseMatrix couplingBlock2 = A_T_F->block(e2.Number*_cellBlockSize, f.Number*_faceBlockSize, _cellBlockSize, _faceBlockSize);
		double e2Kappa = abs(couplingBlock2.trace());

		double minKappa = min(e1Kappa, e2Kappa);
		double maxKappa = max(e1Kappa, e2Kappa);

		return coupling * (minKappa / maxKappa);
	}

	bool IsStronglyCoupled(const HybridAlgebraicElement& e, double coupling)
	{
		return coupling < 0 && coupling < _strongCouplingThreshold * e.Neighbours[0].second;
	}

	/*double DiffusionCoeff(const HybridAlgebraicElement& e)
	{
		DenseMatrix elemBlock = A_T_T->block(e.Number*_cellBlockSize, e.Number*_cellBlockSize, _cellBlockSize, _cellBlockSize);
		double kappa = elemBlock(0, 0);
		return kappa;
	}*/

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
};