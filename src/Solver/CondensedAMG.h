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

	/*SparseMatrix A_T_Tc;
	SparseMatrix A_T_Fc;
	SparseMatrix A_F_Fc;
	SparseMatrix* inv_A_T_Tc;*/

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
		_mesh.Build(*A_T_T, *A_T_F, *A_F_F);
		_mesh.PairWiseAggregate(coarsestPossibleMeshReached);
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
	{
		_mesh.ComputeProlongation();
		this->P = _mesh.P;
	}

	void SetupRestriction() override
	{
		R = P.transpose();
	}

	void OnEndSetup() override
	{
		if (this->CoarserLevel)
		{
			CondensedLevel* coarse = dynamic_cast<CondensedLevel*>(this->CoarserLevel);
			coarse->A_T_T = &_mesh.A_T_Tc;
			coarse->A_T_F = &_mesh.A_T_Fc;
			coarse->A_F_F = &_mesh.A_F_Fc;
			coarse->inv_A_T_T = _mesh.inv_A_T_Tc;
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
		//SparseMatrix* schur = new SparseMatrix(A_F_F - A_T_F.transpose() * (*inv_A_T_T) * A_T_F);
		Multigrid::Setup(A);
	}

	Vector Solve(const Vector& b, string initialGuessCode) override
	{
		if (initialGuessCode.compare("smooth") == 0)
			Utils::Warning("Smooth initial guess unmanaged in AMG.");

		/*CondensedLevel* fine = dynamic_cast<CondensedLevel*>(this->_fineLevel);
		auto T = fine->A_T_T->rows();
		auto F = fine->A_T_F->cols();
		Vector condensedB = b.tail(F) - fine->A_T_F->transpose() * (*fine->inv_A_T_T) * b.head(T);
		return Multigrid::Solve(condensedB, initialGuessCode);*/
		return Multigrid::Solve(b, initialGuessCode);
	}
	
protected:
	Level* CreateCoarseLevel(Level* fineLevel) override
	{
		CondensedLevel* coarseLevel = new CondensedLevel(fineLevel->Number + 1, this->_cellBlockSize, this->_faceBlockSize);
		return coarseLevel;
	}
};