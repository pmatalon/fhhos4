#pragma once
#include <mutex>
#include "CondensedAlgebraicMesh.h"
#include "Multigrid.h"
using namespace std;

class CondensedLevel : public Level
{
private:
	int _cellBlockSize;
	int _faceBlockSize;
	CondensedAlgebraicMesh _mesh;
public:
	const SparseMatrix* A_T_T;
	const SparseMatrix* A_T_F;
	const SparseMatrix* A_F_F;
	const SparseMatrix* inv_A_T_T;
private:
	SparseMatrix Q_T;
	SparseMatrix Q_F;

	SparseMatrix A_T_Tc;
	SparseMatrix A_T_Fc;
	SparseMatrix A_F_Fc;
	SparseMatrix* inv_A_T_Tc;

public:
	CondensedLevel(int number, int cellBlockSize, int faceBlockSize)
		: Level(number), _mesh(cellBlockSize, faceBlockSize)
	{
		this->_cellBlockSize = cellBlockSize;
		this->_faceBlockSize = faceBlockSize;
	}

	BigNumber NUnknowns() override
	{
		return A_F_F->rows();
	}

	void ExportVector(const Vector& v, string suffix, int levelNumber) override
	{ }
	void ExportMatrix(const SparseMatrix& M, string suffix, int levelNumber) override
	{ }

	void CoarsenMesh(CoarseningStrategy coarseningStgy, int coarseningFactor, bool& noCoarserMeshProvided, bool& coarsestPossibleMeshReached) override
	{
		_mesh.Build(*A_T_T, *A_T_F);

		//----------------------------//
		// First pairwise aggregation //
		//----------------------------//

		_mesh.PairWiseAggregate(coarsestPossibleMeshReached);
		if (coarsestPossibleMeshReached)
			return;

		//-------------------------------//
		// Intermediate coarse operators //
		//-------------------------------//

		// Prolongation operator with only one 1 coefficient per row
		SparseMatrix Q_Ti = BuildQ_T(_mesh);
		SparseMatrix Q_Fi = BuildQ_F(_mesh);

		// Intermediate coarse operators
		SparseMatrix A_T_Ti = Q_Ti.transpose() * (*A_T_T) * Q_Ti;
		SparseMatrix A_T_Fi = Q_Ti.transpose() * (*A_T_F) * Q_Fi;

		//-----------------------------//
		// Second pairwise aggregation //
		//-----------------------------//

		CondensedAlgebraicMesh coarseMesh(_cellBlockSize, _faceBlockSize);
		coarseMesh.Build(A_T_Ti, A_T_Fi);
		coarseMesh.PairWiseAggregate(coarsestPossibleMeshReached);
		if (coarsestPossibleMeshReached)
			return;

		// Assembly of the double-aggregates
		vector<vector<AlgebraicElement*>> doubleElemAggregates(coarseMesh._coarseElements.size());
		NumberParallelLoop<EmptyResultChunk> parallelLoopE(coarseMesh._coarseElements.size());
		parallelLoopE.Execute([this, &coarseMesh, &doubleElemAggregates](BigNumber secondAggregNumber)
			{
				ElementAggregate& agg2 = coarseMesh._coarseElements[secondAggregNumber];
				vector<AlgebraicElement*> doubleAggregate;
				for (AlgebraicElement* ce1 : agg2.FineElements)
				{
					ElementAggregate& agg1 = _mesh._coarseElements[ce1->Number];
					for (AlgebraicElement* fe : agg1.FineElements)
						doubleAggregate.push_back(fe);
				}
				doubleElemAggregates[secondAggregNumber] = doubleAggregate;
			});

		vector<vector<AlgebraicFace*>> doubleFaceAggregates(coarseMesh._coarseFaces.size());
		NumberParallelLoop<EmptyResultChunk> parallelLoopF(coarseMesh._coarseFaces.size());
		parallelLoopF.Execute([this, &coarseMesh, &doubleFaceAggregates](BigNumber secondAggregNumber)
			{
				FaceAggregate& agg2 = coarseMesh._coarseFaces[secondAggregNumber];
				vector<AlgebraicFace*> doubleAggregate;
				for (AlgebraicFace* cf1 : agg2.FineFaces)
				{
					FaceAggregate& agg1 = _mesh._coarseFaces[cf1->Number];
					for (AlgebraicFace* fe : agg1.FineFaces)
						doubleAggregate.push_back(fe);
				}
				doubleFaceAggregates[secondAggregNumber] = doubleAggregate;
			});

		//------------------//
		// Coarse operators //
		//------------------//

		SparseMatrix Q_T = BuildQ_T(doubleElemAggregates);
		SparseMatrix Q_F = BuildQ_F(doubleFaceAggregates);

		this->A_T_Tc = Q_T.transpose() * (*A_T_T) * Q_T;
		this->A_T_Fc = Q_T.transpose() * (*A_T_F) * Q_F;

		//-----------------------//
		// Prolongation operator //
		//-----------------------//

		// Theta: reconstruction from the coarse faces to the coarse cells
		this->inv_A_T_Tc = new SparseMatrix(Utils::InvertBlockDiagMatrix(A_T_Tc, _cellBlockSize));
		SparseMatrix Theta = -(*inv_A_T_Tc) * A_T_Fc;

		// Pi: average on both sides of each face
		DenseMatrix traceOfConstant = DenseMatrix::Zero(_faceBlockSize, _cellBlockSize);
		traceOfConstant(0, 0) = 1;
		NumberParallelLoop<CoeffsChunk> parallelLoopPi(_mesh._faces.size());
		parallelLoopPi.Execute([this, &traceOfConstant](BigNumber faceNumber, ParallelChunk<CoeffsChunk>* chunk)
			{
				AlgebraicFace& face = _mesh._faces[faceNumber];
				for (AlgebraicElement* elem : face.Elements)
					chunk->Results.Coeffs.Add(faceNumber, elem->Number, 1.0 / face.Elements.size()*traceOfConstant);
			});
		SparseMatrix Pi = SparseMatrix(_mesh._faces.size()*_faceBlockSize, _mesh._elements.size()*_cellBlockSize);
		parallelLoopPi.Fill(Pi);

		// Prolongation
		this->P = Pi * Q_T * Theta;

		if (true)
			this->A_F_Fc = Q_F.transpose() * (*A_F_F) * Q_F;
		else
			this->A_F_Fc = P.transpose() * (*A_F_F) * P;
	}

private:
	// Cell prolongation Q_T with only one 1 coefficient per row
	SparseMatrix BuildQ_T(const CondensedAlgebraicMesh& mesh)
	{
		DenseMatrix Id = DenseMatrix::Identity(_cellBlockSize, _cellBlockSize);
		NumberParallelLoop<CoeffsChunk> parallelLoopQ_T(mesh._elements.size());
		parallelLoopQ_T.Execute([this, &mesh, &Id](BigNumber elemNumber, ParallelChunk<CoeffsChunk>* chunk)
			{
				const AlgebraicElement& elem = mesh._elements[elemNumber];
				chunk->Results.Coeffs.Add(elem.Number, elem.CoarseElement->Number, Id);
			});
		SparseMatrix Q_T = SparseMatrix(mesh._elements.size()*_cellBlockSize, mesh._coarseElements.size()*_cellBlockSize);
		parallelLoopQ_T.Fill(Q_T);
		return Q_T;
	}

	// Face prolongation Q_F with only one 1 coefficient per row
	SparseMatrix BuildQ_F(const CondensedAlgebraicMesh& mesh)
	{
		DenseMatrix Id = DenseMatrix::Identity(_faceBlockSize, _faceBlockSize);
		NumberParallelLoop<CoeffsChunk> parallelLoopQ_F(mesh._coarseFaces.size());
		parallelLoopQ_F.Execute([this, &mesh, &Id](BigNumber faceAggNumber, ParallelChunk<CoeffsChunk>* chunk)
			{
				const FaceAggregate* agg = &mesh._coarseFaces[faceAggNumber];
				for (AlgebraicFace* face : agg->FineFaces)
					chunk->Results.Coeffs.Add(face->Number, agg->Number, Id);
			});
		SparseMatrix Q_F = SparseMatrix(mesh._faces.size()*_faceBlockSize, mesh._coarseFaces.size()*_faceBlockSize);
		parallelLoopQ_F.Fill(Q_F);
		return Q_F;
	}


	// Cell prolongation Q_T with only one 1 coefficient per row
	SparseMatrix BuildQ_T(const vector<vector<AlgebraicElement*>>& aggregates)
	{
		DenseMatrix Id = DenseMatrix::Identity(_cellBlockSize, _cellBlockSize);
		NumberParallelLoop<CoeffsChunk> parallelLoop(aggregates.size());
		parallelLoop.Execute([this, &aggregates, &Id](BigNumber aggregNumber, ParallelChunk<CoeffsChunk>* chunk)
			{
				const vector<AlgebraicElement*>& aggreg = aggregates[aggregNumber];
				for (AlgebraicElement* e : aggreg)
					chunk->Results.Coeffs.Add(e->Number, aggregNumber, Id);
			});
		SparseMatrix Q_T = SparseMatrix(_mesh._elements.size()*_cellBlockSize, aggregates.size()*_cellBlockSize);
		parallelLoop.Fill(Q_T);
		return Q_T;
	}

	// Face prolongation Q_F with only one 1 coefficient per row
	SparseMatrix BuildQ_F(vector<vector<AlgebraicFace*>> aggregates)
	{
		DenseMatrix Id = DenseMatrix::Identity(_faceBlockSize, _faceBlockSize);
		NumberParallelLoop<CoeffsChunk> parallelLoop(aggregates.size());
		parallelLoop.Execute([this, &aggregates, &Id](BigNumber aggregNumber, ParallelChunk<CoeffsChunk>* chunk)
			{
				const vector<AlgebraicFace*>& aggreg = aggregates[aggregNumber];
				for (AlgebraicFace* face : aggreg)
					chunk->Results.Coeffs.Add(face->Number, aggregNumber, Id);
			});
		SparseMatrix Q_F = SparseMatrix(_mesh._faces.size()*_faceBlockSize, aggregates.size()*_faceBlockSize);
		parallelLoop.Fill(Q_F);
		return Q_F;
	}
public:
	void SetupDiscretizedOperator() override 
	{
		SparseMatrix* schur = new SparseMatrix(*A_F_F - (A_T_F->transpose()) * (*inv_A_T_T) * (*A_T_F));
		this->OperatorMatrix = schur;
	}

	void OnStartSetup() override
	{}

	void SetupProlongation() override
	{}

	void SetupRestriction() override
	{
		R = P.transpose();
	}

	void OnEndSetup() override
	{
		if (this->CoarserLevel)
		{
			CondensedLevel* coarse = dynamic_cast<CondensedLevel*>(this->CoarserLevel);
			coarse->A_T_T = &A_T_Tc;
			coarse->A_T_F = &A_T_Fc;
			coarse->A_F_F = &A_F_Fc;
			coarse->inv_A_T_T = inv_A_T_Tc;
		}
	}

	~CondensedLevel()
	{
		if (!this->IsFinestLevel())
			delete OperatorMatrix;
	}
};

class CondensedAMG : public Multigrid
{
private:
	int _cellBlockSize;
	int _faceBlockSize;
public:

	CondensedAMG(int cellBlockSize, int faceBlockSize, int nLevels = 0)
		: Multigrid(nLevels)
	{
		this->_cellBlockSize = cellBlockSize;
		this->_faceBlockSize = faceBlockSize;
		this->BlockSizeForBlockSmoothers = faceBlockSize;
		this->UseGalerkinOperator = true;
		this->_fineLevel = new CondensedLevel(0, cellBlockSize, faceBlockSize);
	}

	void BeginSerialize(ostream& os) const override
	{
		os << "CondensedAMG" << endl;
	}

	void EndSerialize(ostream& os) const override
	{
	}

	void Setup(const SparseMatrix& A) override
	{
		assert(false && "This Setup method cannot be used in this solver.");
	}

	void Setup(const SparseMatrix& A, const SparseMatrix& A_T_T, const SparseMatrix& A_T_F, const SparseMatrix& A_F_F) override
	{
		CondensedLevel* fine = dynamic_cast<CondensedLevel*>(this->_fineLevel);
		fine->A_T_T = &A_T_T;
		fine->A_T_F = &A_T_F;
		fine->A_F_F = &A_F_F;
		SparseMatrix* inv_A_T_T = new SparseMatrix(Utils::InvertBlockDiagMatrix(A_T_T, _cellBlockSize));
		fine->inv_A_T_T = inv_A_T_T;
		Multigrid::Setup(A);
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
		CondensedLevel* coarseLevel = new CondensedLevel(fineLevel->Number + 1, this->_cellBlockSize, this->_faceBlockSize);
		return coarseLevel;
	}
};