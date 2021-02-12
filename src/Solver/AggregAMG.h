#pragma once
#include <mutex>
#include "AlgebraicMesh.h"
#include "Multigrid.h"
using namespace std;

class AggregLevel : public Level
{
private:
	int _blockSize;
	double _strongCouplingThreshold;
	AlgebraicMesh _mesh;
	vector<ElementAggregate> _aggregates;
	vector<vector<BigNumber>> _numberAggregates;
	double _restrictionCost = 0;
public:
	SparseMatrix Ac;

public:
	AggregLevel(int number, int blockSize, double strongCouplingThreshold)
		: Level(number), _mesh(blockSize, strongCouplingThreshold)
	{
		this->_blockSize = blockSize;
		this->_strongCouplingThreshold = strongCouplingThreshold;
	}

	void ExportVector(const Vector& v, string suffix, int levelNumber) override
	{ }

	void ExportMatrix(const SparseMatrix& M, string suffix, int levelNumber) override
	{
		string file = Utils::ProgramArgs.OutputDirectory + "/" + suffix;
		if (levelNumber > 0)
			file += "_" + to_string(levelNumber);
		file += ".dat";
		Eigen::saveMarket(M, file);
	}

	void CoarsenMesh(CoarseningStrategy coarseningStgy, int coarseningFactor, bool& noCoarserMeshProvided, bool& coarsestPossibleMeshReached) override
	{
		cout << "\tBuild algebraic mesh" << endl;

		//ExportMatrix(*this->OperatorMatrix, "A", 0);

		_mesh.Build(*this->OperatorMatrix);

		//----------------------------//
		// First pairwise aggregation //
		//----------------------------//

		cout << "\tPairwise aggregation 1" << endl;
		_mesh.PairWiseAggregate(coarsestPossibleMeshReached);
		if (coarsestPossibleMeshReached)
			return;

		/*// Cell-prolongation operator with only one 1 coefficient per row
		SparseMatrix Q_T1 = BuildQ_T(_mesh);
		SparseMatrix A1_old = Q_T1.transpose()* *this->OperatorMatrix * Q_T1;*/

		// Intermediate coarse operator
		SparseMatrix A1 = GalerkinOperator(*this->OperatorMatrix, _mesh);

		//-----------------------------//
		// Second pairwise aggregation //
		//-----------------------------//

		cout << "\tBuild intermediary coarse mesh" << endl;

		AlgebraicMesh coarseMesh(_blockSize, _strongCouplingThreshold);
		coarseMesh.FinerMesh = &_mesh;
		coarseMesh.Build(A1);

		cout << "\tPairwise aggregation 2" << endl;

		coarseMesh.PairWiseAggregate(coarsestPossibleMeshReached);
		if (coarsestPossibleMeshReached)
			return;

		/*SparseMatrix Q_T2 = BuildQ_T(coarseMesh);
		SparseMatrix Ac_old = Q_T2.transpose()* A1_old * Q_T2;
		//this->P = Q_T1 * Q_T2;*/

		this->Ac = GalerkinOperator(A1, coarseMesh);

		//-----------------------------------//
		// Assembly of the double-aggregates //
		//-----------------------------------//

		_aggregates = vector<ElementAggregate>(coarseMesh._coarseElements.size());
		_numberAggregates = vector<vector<BigNumber>>(coarseMesh._coarseElements.size());
		NumberParallelLoop<EmptyResultChunk> parallelLoop(coarseMesh._coarseElements.size());
		parallelLoop.Execute([this, &coarseMesh](BigNumber finalAggregNumber)
			{
				ElementAggregate& finalAggreg = _aggregates[finalAggregNumber];
				vector<BigNumber>& finalNumberAggreg = _numberAggregates[finalAggregNumber];
				finalAggreg.Number = finalAggregNumber;
				finalAggreg.FineElements = GetFineElements(coarseMesh._coarseElements[finalAggregNumber], coarseMesh);
				for (AlgebraicElement* e : finalAggreg.FineElements)
				{
					e->FinalAggregate = &finalAggreg;
					e->FinalAggregateNumber = finalAggreg.Number;
					finalNumberAggreg.push_back(e->Number);
				}
			});
		
		_restrictionCost = 0;
		for (ElementAggregate& a : _aggregates)
			_restrictionCost += a.FineElements.size();
	}

private:
	vector<AlgebraicElement*> GetFineElements(const ElementAggregate& aggreg, const AlgebraicMesh& mesh)
	{
		if (!mesh.FinerMesh)
			return aggreg.FineElements;
		
		vector<AlgebraicElement*> fineElements;
		for (AlgebraicElement* e : aggreg.FineElements)
		{
			ElementAggregate& finerAggreg = mesh.FinerMesh->_coarseElements[e->Number];
			fineElements = Utils::Join(fineElements, GetFineElements(finerAggreg, *mesh.FinerMesh));
		}
		return fineElements;
	}

	SparseMatrix GalerkinOperator(const SparseMatrix& A, const AlgebraicMesh& mesh)
	{
		// Formula (2.3) of Notay's 2010 paper
		NonZeroCoefficients coeffs(A.nonZeros());
		for (BigNumber k = 0; k < A.rows(); k++)
		{
			for (RowMajorSparseMatrix::InnerIterator it(A, k); it; ++it)
			{
				auto l = it.col();
				double a_kl = it.value();

				auto i = mesh._elements[k].CoarseElement->Number;
				auto j = mesh._elements[l].CoarseElement->Number;

				coeffs.Add(i, j, a_kl);
			}
		}

		SparseMatrix Ac(mesh._coarseElements.size(), mesh._coarseElements.size());
		coeffs.Fill(Ac);
		return Ac;
	}


	// Cell prolongation Q_T with only one 1 coefficient per row
	SparseMatrix BuildQ_T(const AlgebraicMesh& mesh)
	{
		DenseMatrix Id = DenseMatrix::Identity(_blockSize, _blockSize);
		NumberParallelLoop<CoeffsChunk> parallelLoopQ_T(mesh._elements.size());
		parallelLoopQ_T.Execute([this, &mesh, &Id](BigNumber elemNumber, ParallelChunk<CoeffsChunk>* chunk)
			{
				const AlgebraicElement& elem = mesh._elements[elemNumber];
				chunk->Results.Coeffs.Add(elem.Number*_blockSize, elem.CoarseElement->Number*_blockSize, Id);
			});
		SparseMatrix Q_T = SparseMatrix(mesh._elements.size()*_blockSize, mesh._coarseElements.size()*_blockSize);
		parallelLoopQ_T.Fill(Q_T);
		return Q_T;
	}

public:
	void OnStartSetup() override
	{
		cout << "\t\tMesh                : " << this->NUnknowns() / _blockSize << " elements";
		if (!this->IsFinestLevel())
		{
			double nFine = this->FinerLevel->NUnknowns();
			double nCoarse = this->NUnknowns();
			cout << ", coarsening factor = " << (nFine/nCoarse);
		}
		cout << endl;
	}

	Vector Prolong(Vector& coarseV) override
	{
		if (_blockSize == 1)
		{
			Vector fineV(_mesh._elements.size());
			/*for (AlgebraicElement& e : _mesh._elements)
				fineV[e.Number] = coarseV[e.FinalAggregate->Number];
			return fineV;*/
			for (AlgebraicElement& e : _mesh._elements)
				fineV[e.Number] = coarseV[e.FinalAggregateNumber];
			return fineV;
		}
		else
		{
			Vector fineV = Vector::Zero(_mesh._elements.size() * _blockSize);
			for (AlgebraicElement& e : _mesh._elements)
				fineV[e.Number * _blockSize] = coarseV[e.FinalAggregate->Number * _blockSize];
			return fineV;
		}
	}

	double ProlongCost() override
	{
		return 0;
	}

	Vector Restrict(Vector& fineV) override
	{
		if (_blockSize == 1)
		{
			/*Vector coarseV(_aggregates.size());
			for (ElementAggregate& a : _aggregates)
			{
				coarseV[a.Number] = 0;
				for (AlgebraicElement* e : a.FineElements)
					coarseV[a.Number] += fineV[e->Number];
			}
			return coarseV;*/
			Vector coarseV(_numberAggregates.size());
			for (BigNumber i = 0; i < _numberAggregates.size(); i++)
			{
				vector<BigNumber> fineElements = _numberAggregates[i];
				coarseV[i] = 0;
				for (BigNumber fineI : fineElements)
					coarseV[i] += fineV[fineI];
			}
			return coarseV;

		}
		else
		{
			Vector coarseV = Vector::Zero(_aggregates.size() * _blockSize);
			for (ElementAggregate& a : _aggregates)
			{
				coarseV[a.Number * _blockSize] = 0;
				for (AlgebraicElement* e : a.FineElements)
					coarseV[a.Number * _blockSize] += fineV[e->Number * _blockSize];
			}
			return coarseV;
		}
	}

	double RestrictCost() override
	{
		return _restrictionCost;
	}
};

class AggregAMG : public Multigrid
{
private:
	int _blockSize;
	double _strongCouplingThreshold;
public:

	AggregAMG(int blockSize, double strongCouplingThreshold, int nLevels = 0)
		: Multigrid(nLevels)
	{
		this->_blockSize = blockSize;
		this->_strongCouplingThreshold = strongCouplingThreshold;
		this->BlockSizeForBlockSmoothers = blockSize;
		this->UseGalerkinOperator = true;
		this->_fineLevel = new AggregLevel(0, blockSize, strongCouplingThreshold);
	}

	void BeginSerialize(ostream& os) const override
	{
		os << "AggregAMG" << endl;
	}

	void EndSerialize(ostream& os) const override
	{
	}

	Vector Solve(const Vector& b, string initialGuessCode) override
	{
		if (initialGuessCode.compare("smooth") == 0)
			Utils::Warning("Smooth initial guess unmanaged in AMG.");
		return Multigrid::Solve(b, initialGuessCode);
	}
	
protected:
	Level* CreateCoarseLevel(Level* fineLevel) override
	{
		AggregLevel* coarseLevel = new AggregLevel(fineLevel->Number + 1, _blockSize, _strongCouplingThreshold);
		coarseLevel->OperatorMatrix = &dynamic_cast<AggregLevel*>(fineLevel)->Ac;
		return coarseLevel;
	}
};