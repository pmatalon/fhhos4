#pragma once
#include <mutex>
#include "AlgebraicMesh.h"
#include "Multigrid.h"
using namespace std;

class AggregLevel : public Level
{
private:
	int _blockSize;
	AlgebraicMesh _mesh;
//public:
	//const SparseMatrix* A;
private:
	SparseMatrix Q_T;
	SparseMatrix Q_F;

	//SparseMatrix Ac;

public:
	AggregLevel(int number, int blockSize)
		: Level(number), _mesh(blockSize)
	{
		this->_blockSize = blockSize;
	}

	/*BigNumber NUnknowns() override
	{
		return A->rows();
	}*/

	void ExportVector(const Vector& v, string suffix, int levelNumber) override
	{ }

	void ExportMatrix(const SparseMatrix& M, string suffix, int levelNumber) override
	{ 
		Eigen::saveMarket(M, Utils::ProgramArgs.OutputDirectory + "/" + suffix + ".dat");
	}

	void CoarsenMesh(CoarseningStrategy coarseningStgy, int coarseningFactor, bool& noCoarserMeshProvided, bool& coarsestPossibleMeshReached) override
	{
		if (!this->OperatorMatrix)
			ComputeGalerkinOperator();

		cout << "\tDouble pairwise aggregation" << endl;

		_mesh.Build(*this->OperatorMatrix);

		//ExportMatrix(*A_F_F, "A_F_F", 0);

		//----------------------------//
		// First pairwise aggregation //
		//----------------------------//

		_mesh.PairWiseAggregate(coarsestPossibleMeshReached);
		if (coarsestPossibleMeshReached)
			return;

		/*cout << "------------- Elem" << endl;
		for (HybridElementAggregate& agg : _mesh._coarseElements)
		{
			cout << "(";
			for (auto e : agg.FineElements)
				cout << e->Number << ", ";
			cout << ")" << endl;
		}

		cout << "------------- Faces" << endl;
		for (HybridFaceAggregate& agg : _mesh._coarseFaces)
		{
			cout << "(";
			for (auto f : agg.FineFaces)
				cout << f->Number << ", ";
			cout << ")" << endl;
		}*/

		// Cell-prolongation operator with only one 1 coefficient per row
		SparseMatrix Q_T1 = BuildQ_T(_mesh);
		
		// Intermediate coarse operators
		SparseMatrix A1 = Q_T1.transpose() * (*this->OperatorMatrix) * Q_T1;

		//-----------------------------//
		// Second pairwise aggregation //
		//-----------------------------//

		AlgebraicMesh coarseMesh(_blockSize);
		coarseMesh.Build(A1);
		coarseMesh.PairWiseAggregate(coarsestPossibleMeshReached);
		if (coarsestPossibleMeshReached)
			return;

		SparseMatrix Q_T2 = BuildQ_T(coarseMesh);



//		this->Ac = Q_T2.transpose() * A1 * Q_T2;


		this->P = Q_T1 * Q_T2;
	}

private:
	// Cell prolongation Q_T with only one 1 coefficient per row
	SparseMatrix BuildQ_T(const AlgebraicMesh& mesh)
	{
		DenseMatrix Id = DenseMatrix::Identity(_blockSize, _blockSize);
		NumberParallelLoop<CoeffsChunk> parallelLoopQ_T(mesh._elements.size());
		parallelLoopQ_T.Execute([this, &mesh, &Id](BigNumber elemNumber, ParallelChunk<CoeffsChunk>* chunk)
			{
				const AlgebraicElement& elem = mesh._elements[elemNumber];
				chunk->Results.Coeffs.Add(elem.Number, elem.CoarseElement->Number, Id);
			});
		SparseMatrix Q_T = SparseMatrix(mesh._elements.size()*_blockSize, mesh._coarseElements.size()*_blockSize);
		parallelLoopQ_T.Fill(Q_T);
		return Q_T;
	}

public:
	/*void SetupDiscretizedOperator() override 
	{
		//SparseMatrix* schur = new SparseMatrix(*A_F_F - (A_T_F->transpose()) * (*inv_A_T_T) * (*A_T_F));
		//this->OperatorMatrix = schur;
	}*/

	void OnStartSetup() override
	{
		cout << "\t\tMesh                : " << this->NUnknowns() / _blockSize << " elements";
		if (!this->IsFinestLevel())
		{
			AggregLevel* fine = dynamic_cast<AggregLevel*>(this->FinerLevel);
			double nFine = fine->NUnknowns();
			double nCoarse = this->NUnknowns();
			cout << ", coarsening factor = " << (nFine/nCoarse);
		}
		cout << endl;
	}

	void SetupProlongation() override
	{}

	void SetupRestriction() override
	{
		double scalingFactor = 1.0;
		R = scalingFactor * P.transpose();
	}

	/*void OnEndSetup() override
	{
		if (this->CoarserLevel)
		{
			AggregLevel* coarse = dynamic_cast<AggregLevel*>(this->CoarserLevel);
			coarse->A = &Ac;
		}
	}

	~AggregLevel()
	{
		if (!this->IsFinestLevel() && !this->UseGalerkinOperator)
			delete OperatorMatrix;
	}*/
};

class AggregAMG : public Multigrid
{
private:
	int _blockSize;
public:

	AggregAMG(int blockSize, int nLevels = 0)
		: Multigrid(nLevels)
	{
		this->_blockSize = blockSize;
		this->BlockSizeForBlockSmoothers = blockSize;
		this->UseGalerkinOperator = true;
		this->_fineLevel = new AggregLevel(0, blockSize);
	}

	void BeginSerialize(ostream& os) const override
	{
		os << "AggregAMG" << endl;
	}

	void EndSerialize(ostream& os) const override
	{
	}

	/*void Setup(const SparseMatrix& A) override
	{
		assert(false && "This Setup method cannot be used in this solver.");
	}

	void Setup(const SparseMatrix& A) override
	{
		AggregLevel* fine = dynamic_cast<AggregLevel*>(this->_fineLevel);
		fine->A = &A;
		Multigrid::Setup(A);
	}*/

	Vector Solve(const Vector& b, string initialGuessCode) override
	{
		if (initialGuessCode.compare("smooth") == 0)
			Utils::Warning("Smooth initial guess unmanaged in AMG.");
		return Multigrid::Solve(b, initialGuessCode);
	}
	
protected:
	Level* CreateCoarseLevel(Level* fineLevel) override
	{
		AggregLevel* coarseLevel = new AggregLevel(fineLevel->Number + 1, _blockSize);
		return coarseLevel;
	}
};