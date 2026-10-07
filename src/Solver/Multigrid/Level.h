#pragma once
#include "Smoother.h"
#include "SmootherFactory.h"
#include "../Krylov/FlexibleConjugateGradient.h"
#include "../../Utils/ExportModule.h"
using namespace std;

class Level
{
public:
	int Number;
	const SparseMatrix* OperatorMatrix = nullptr;
	bool UseGalerkinOperator = false;
private:
	SparseMatrix _galerkinOperator;

public:
	Smoother* PreSmoother = nullptr;
	Smoother* PostSmoother = nullptr;

	Level* FinerLevel = nullptr;
	Level* CoarserLevel = nullptr;
	CoarseningType ComesFrom = CoarseningType::H;

	bool ExportComponents = false;
	BigNumber SetupComputationalWork = 0;

	//FlexibleConjugateGradient* FCG = nullptr; // used in K-cycle

protected:
	ExportModule Out;
	SparseMatrix R;
	SparseMatrix P;

public:
	Level(int number)
	{
		this->Number = number;
	}

	bool IsFinestLevel()
	{
		return this->FinerLevel == nullptr;
	}

	bool IsCoarsestLevel()
	{
		return this->CoarserLevel == nullptr;
	}

	void SetExportModule(const ExportModule& out)
	{
		this->Out = ExportModule(out);
		this->Out.AddFilePrefix("level" + to_string(this->Number));
	}

	void Setup()
	{
		Utils::Log() << "\tSetup level " << this->Number << endl;

		OnStartSetup();

		if (this->IsFinestLevel())
		{
			Utils::Log() << "\t\tFine grid operator  : "; Utils::Log().flush();
			Utils::Log() << Utils::MatrixInfo(*this->OperatorMatrix, "A") << endl;
		}
		else
		{
			if (this->UseGalerkinOperator)
			{
				Utils::Log() << "\t\tGalerkin operator   : "; Utils::Log().flush();
				if (!this->OperatorMatrix)
					ComputeGalerkinOperator();
			}
			else
			{
				Utils::Log() << "\t\tDiscretized operator: "; Utils::Log().flush();
				SetupDiscretizedOperator();
			}
			Utils::Log() << Utils::MatrixInfo(*this->OperatorMatrix, "A") << endl;
		}

		const SparseMatrix &A = *this->OperatorMatrix;

		if (!this->IsCoarsestLevel())
		{
			Utils::Log() << "\t\tProlongation        : "; Utils::Log().flush();
			SetupProlongation();
			if (this->P.rows() == 0)
				Utils::Log() << "matrix free" << endl;
			else
				Utils::Log() << Utils::MatrixInfo(this->P, "P") << endl;

			Utils::Log() << "\t\tRestriction         : "; Utils::Log().flush();
			SetupRestriction();
			if (this->R.rows() == 0)
				Utils::Log() << "matrix free" << endl;
			else
				Utils::Log() << Utils::MatrixInfo(this->R, "R") << endl;

			Utils::Log() << "\t\tPreSmoothing        : "; Utils::Log().flush();
			PreSmoother->Setup(A);
			if (PreSmoother->Iterations() == 0)
				Utils::Log() << "none" << endl;
			else
				Utils::Log() << PreSmoother->Iterations() << " iteration" << (PreSmoother->Iterations() > 1 ? "s" : "") << " of " << (*PreSmoother) << endl;

			Utils::Log() << "\t\tPostSmoothing       : "; Utils::Log().flush();
			PostSmoother->Setup(A);
			if (PostSmoother->Iterations() == 0)
				Utils::Log() << "none" << endl;
			else
				Utils::Log() << PostSmoother->Iterations() << " iteration" << (PostSmoother->Iterations() > 1 ? "s" : "") << " of " << (*PostSmoother) << endl;
		}

		if (ExportComponents)
		{
			this->ExportMatrix(A, "A");
			
			if (!this->IsCoarsestLevel())
			{
				this->ExportMatrix(P, "P");
				this->ExportMatrix(R, "R");
			}
		}

		//if (FCG)
			//FCG->Setup(A);

		OnEndSetup();
	}

	void ComputeGalerkinOperator()
	{
		this->_galerkinOperator = FinerLevel->R * *(FinerLevel->OperatorMatrix) * FinerLevel->P;
		this->OperatorMatrix = &this->_galerkinOperator;
	}

	virtual Vector Restrict(Vector& vectorOnThisLevel)
	{
		Vector coarseVector = R * vectorOnThisLevel;
		return coarseVector;
	}

	virtual Flops RestrictCost()
	{
		return 2 * R.nonZeros();
	}

	virtual Vector Prolong(Vector& vectorOnTheCoarserLevel)
	{
		Vector vectorOnThisLevel = P * vectorOnTheCoarserLevel;
		return vectorOnThisLevel;
	}

	virtual Flops ProlongCost()
	{
		return 2 * P.nonZeros();
	}

	// Used for full Neumann boundary conditions. Fixes the constant in order to have a well-posed problem.
	virtual void ApplyZeroMeanCondition(Vector& x) {}
	virtual Flops ApplyZeroMeanConditionCost(Vector& x) { return 0; }

	// Used for full Neumann boundary conditions. Ensures solvability (existence of a solution)
	virtual void EnforceCompatibilityCondition(Vector& b) {}
	virtual Flops EnforceCompatibilityConditionCost(Vector& b) { return 0; }

	virtual BigNumber NUnknowns()
	{
		if (this->OperatorMatrix)
			return this->OperatorMatrix->rows();
		else if (!this->IsFinestLevel())
			return this->FinerLevel->P.cols();
		Utils::FatalError("Level::NUnknowns() must be overridden to use this level");
		return 0;
	}

	virtual int PolynomialDegree() = 0;

	virtual int BlockSizeForBlockSmoothers()
	{
		Utils::FatalError("Level::BlockSizeForBlockSmoothers() must be overridden to use this level in a p-multigrid");
		return 0;
	}

	virtual void CoarsenMesh(H_CoarsStgy coarseningStgy, FaceCoarseningStrategy faceCoarseningStgy, FaceCollapsing bdryFaceCollapsing, double coarseningFactor, bool& noCoarserMeshProvided, bool& coarsestPossibleMeshReached)
	{
		Utils::FatalError("Level::CoarsenMesh() must be overridden to use this level in an h-multigrid");
	}

	void ExportVector(const Vector& v, string suffix)
	{
		Out.ExportVector(v, suffix);
	}
	void ExportMatrix(const SparseMatrix& M, string suffix)
	{
		Out.ExportMatrix(M, suffix);
	}

	virtual Smoother* CreateSmoother(string smootherCode, int nSmootherIterations, int blockSize, double omega)
	{
		return SmootherFactory::Create(smootherCode, nSmootherIterations, blockSize, omega);
	}

	virtual ~Level()
	{
		delete PreSmoother;
		delete PostSmoother;
		if (CoarserLevel)
			delete CoarserLevel;
		//if (FCG)
			//delete FCG;
	}

protected:
	virtual void OnStartSetup() {}
	virtual void OnEndSetup() {}
	virtual void SetupDiscretizedOperator() {}
	virtual void SetupRestriction() {}
	virtual void SetupProlongation() {}
};