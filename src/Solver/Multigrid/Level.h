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
	// On several processes (MPI): the operator distributed by rows (OperatorMatrix is then its block Diag); null on one
	// process
	const DistributedMatrix* DistributedOperator = nullptr;
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
	// On several processes: the prolongation, its rows the owned unknowns of this level, its columns those of the coarser
	// level; the restriction is its transpose (MultiplyTransposed())
	DistributedMatrix DistributedP;

public:
	Level(int number)
	{
		this->Number = number;
	}

	bool IsFinestLevel()
	{
		return this->FinerLevel == nullptr;
	}

	bool IsDistributed() const
	{
		return this->DistributedOperator != nullptr;
	}

	// b - A x, A x with the operator (distributed or not)
	Vector ResidualOf(const Vector& b, const Vector& x) const
	{
		if (this->DistributedOperator)
			return b - *this->DistributedOperator * x;
		return b - *this->OperatorMatrix * x;
	}

	Vector Apply(const Vector& x) const
	{
		if (this->DistributedOperator)
			return *this->DistributedOperator * x;
		return *this->OperatorMatrix * x;
	}

	// The number of non-zeros of the operator (of all the processes)
	BigNumber OperatorNonZeros() const
	{
		if (this->DistributedOperator)
			return this->DistributedOperator->GlobalNonZeros();
		return this->OperatorMatrix->nonZeros();
	}

	string OperatorInfo() const
	{
		if (!this->DistributedOperator)
			return Utils::MatrixInfo(*this->OperatorMatrix, "A");
		const DistributedMatrix& A = *this->DistributedOperator;
		BigNumber nnz = A.GlobalNonZeros();
		return "size(A)=" + to_string(A.GlobalRows()) + "x" + to_string(A.GlobalCols()) + ", \tnnz(A)=" + to_string(nnz) + " on " + to_string(A.Comm().Size()) + " processes";
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
			Utils::Log() << OperatorInfo() << endl;
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
			Utils::Log() << OperatorInfo() << endl;
		}

		const SparseMatrix &A = *this->OperatorMatrix;

		if (!this->IsCoarsestLevel())
		{
			Utils::Log() << "\t\tProlongation        : "; Utils::Log().flush();
			SetupProlongation();
			if (this->IsDistributed())
				Utils::Log() << (this->DistributedP.Diag.rows() == 0 && this->P.rows() == 0 ? "matrix free" : "distributed") << endl;
			else if (this->P.rows() == 0)
				Utils::Log() << "matrix free" << endl;
			else
				Utils::Log() << Utils::MatrixInfo(this->P, "P") << endl;

			Utils::Log() << "\t\tRestriction         : "; Utils::Log().flush();
			SetupRestriction();
			if (this->IsDistributed())
				Utils::Log() << (this->DistributedP.Diag.rows() == 0 && this->R.rows() == 0 ? "matrix free" : "transpose of the prolongation") << endl;
			else if (this->R.rows() == 0)
				Utils::Log() << "matrix free" << endl;
			else
				Utils::Log() << Utils::MatrixInfo(this->R, "R") << endl;

			Utils::Log() << "\t\tPreSmoothing        : "; Utils::Log().flush();
			if (this->IsDistributed())
				PreSmoother->Setup(*this->DistributedOperator);
			else
				PreSmoother->Setup(A);
			if (PreSmoother->Iterations() == 0)
				Utils::Log() << "none" << endl;
			else
				Utils::Log() << PreSmoother->Iterations() << " iteration" << (PreSmoother->Iterations() > 1 ? "s" : "") << " of " << (*PreSmoother) << endl;

			Utils::Log() << "\t\tPostSmoothing       : "; Utils::Log().flush();
			if (this->IsDistributed())
				PostSmoother->Setup(*this->DistributedOperator, *PreSmoother);
			else
				PostSmoother->Setup(A, *PreSmoother); // reuses the setup of the pre-smoother if they are of the same kind
			if (PostSmoother->Iterations() == 0)
				Utils::Log() << "none" << endl;
			else
				Utils::Log() << PostSmoother->Iterations() << " iteration" << (PostSmoother->Iterations() > 1 ? "s" : "") << " of " << (*PostSmoother) << endl;
		}

		if (ExportComponents && this->IsDistributed())
			Utils::FatalError("The export of the multigrid's matrices is not available on several processes.");
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
		if (this->IsDistributed())
			return this->DistributedP.MultiplyTransposed(vectorOnThisLevel);
		Vector coarseVector = R * vectorOnThisLevel;
		return coarseVector;
	}

	virtual Flops RestrictCost()
	{
		if (this->IsDistributed())
			return 2 * this->DistributedP.LocalNonZeros();
		return 2 * R.nonZeros();
	}

	virtual Vector Prolong(Vector& vectorOnTheCoarserLevel)
	{
		if (this->IsDistributed())
			return this->DistributedP * vectorOnTheCoarserLevel;
		Vector vectorOnThisLevel = P * vectorOnTheCoarserLevel;
		return vectorOnThisLevel;
	}

	virtual Flops ProlongCost()
	{
		if (this->IsDistributed())
			return 2 * this->DistributedP.LocalNonZeros();
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

	virtual Smoother* CreateSmoother(SmootherType smoother, int nSmootherIterations, int blockSize, double omega)
	{
		return SmootherFactory::Create(smoother, nSmootherIterations, blockSize, omega);
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