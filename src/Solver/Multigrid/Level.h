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
	// Mixed precision: the operator of a level in single precision (U-AMG, Multigrid::MixedPrecision), whose smoothers
	// and transfers are in single precision; OperatorMatrix is then null, but on the finest level (the matrix of the
	// solver, in double) and on the coarsest one (a copy in double, for the coarse solver). Null otherwise.
	const FloatSparseMatrix* SingleOperatorMatrix = nullptr;
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
	FloatSparseMatrix SingleR, SingleP; // a level in single precision (SingleOperatorMatrix): its transfers
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

	bool InSinglePrecision() const
	{
		return this->SingleOperatorMatrix != nullptr;
	}

	// The operator, the prolongation and the restriction in precision T (float: a level in single precision)
	template <typename T>
	const SparseMatrixOf<T>*& OperatorOf()
	{
		if constexpr (std::is_same_v<T, float>)
			return this->SingleOperatorMatrix;
		else
			return this->OperatorMatrix;
	}

	template <typename T>
	SparseMatrixOf<T>& ProlongationOf()
	{
		if constexpr (std::is_same_v<T, float>)
			return this->SingleP;
		else
			return this->P;
	}

	template <typename T>
	SparseMatrixOf<T>& RestrictionOf()
	{
		if constexpr (std::is_same_v<T, float>)
			return this->SingleR;
		else
			return this->R;
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

	// A level in single precision
	FloatVector ResidualOf(const FloatVector& b, const FloatVector& x) const
	{
		return b - *this->SingleOperatorMatrix * x;
	}

	FloatVector Apply(const FloatVector& x) const
	{
		return *this->SingleOperatorMatrix * x;
	}

	// The number of non-zeros of the operator (of all the processes)
	BigNumber OperatorNonZeros() const
	{
		if (this->DistributedOperator)
			return this->DistributedOperator->GlobalNonZeros();
		if (!this->OperatorMatrix)
			return this->SingleOperatorMatrix->nonZeros();
		return this->OperatorMatrix->nonZeros();
	}

	string OperatorInfo() const
	{
		if (!this->OperatorMatrix && this->SingleOperatorMatrix)
			return Utils::MatrixInfo(*this->SingleOperatorMatrix, "A") + " (single precision)";
		if (!this->DistributedOperator)
			return Utils::MatrixInfo(*this->OperatorMatrix, "A") + (this->SingleOperatorMatrix ? " (and in single precision)" : "");
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
				if (!this->OperatorMatrix && !this->SingleOperatorMatrix)
					ComputeGalerkinOperator();
			}
			else
			{
				Utils::Log() << "\t\tDiscretized operator: "; Utils::Log().flush();
				SetupDiscretizedOperator();
			}
			Utils::Log() << OperatorInfo() << endl;
		}

		if (!this->IsCoarsestLevel())
		{
			Utils::Log() << "\t\tProlongation        : "; Utils::Log().flush();
			SetupProlongation();
			Utils::Log() << ProlongationInfo() << endl;

			Utils::Log() << "\t\tRestriction         : "; Utils::Log().flush();
			SetupRestriction();
			Utils::Log() << RestrictionInfo() << endl;

			Utils::Log() << "\t\tPreSmoothing        : "; Utils::Log().flush();
			if (this->IsDistributed())
				PreSmoother->Setup(*this->DistributedOperator);
			else if (this->InSinglePrecision())
				PreSmoother->Setup(*this->SingleOperatorMatrix);
			else
				PreSmoother->Setup(*this->OperatorMatrix);
			if (PreSmoother->Iterations() == 0)
				Utils::Log() << "none" << endl;
			else
				Utils::Log() << PreSmoother->Iterations() << " iteration" << (PreSmoother->Iterations() > 1 ? "s" : "") << " of " << (*PreSmoother) << endl;

			Utils::Log() << "\t\tPostSmoothing       : "; Utils::Log().flush();
			if (this->IsDistributed())
				PostSmoother->Setup(*this->DistributedOperator, *PreSmoother);
			else if (this->InSinglePrecision())
				PostSmoother->Setup(*this->SingleOperatorMatrix, *PreSmoother);
			else
				PostSmoother->Setup(*this->OperatorMatrix, *PreSmoother); // reuses the setup of the pre-smoother if they are of the same kind
			if (PostSmoother->Iterations() == 0)
				Utils::Log() << "none" << endl;
			else
				Utils::Log() << PostSmoother->Iterations() << " iteration" << (PostSmoother->Iterations() > 1 ? "s" : "") << " of " << (*PostSmoother) << endl;
		}

		if (ExportComponents && (this->IsDistributed() || this->InSinglePrecision()))
			Utils::FatalError("The export of the multigrid's matrices is not available on several processes nor in single precision.");
		if (ExportComponents)
		{
			this->ExportMatrix(*this->OperatorMatrix, "A");
			
			if (!this->IsCoarsestLevel())
				ExportTransfers();
		}

		//if (FCG)
			//FCG->Setup(A);

		OnEndSetup();
	}

	// The printouts of the setup: the prolongation and the restriction
	virtual string ProlongationInfo()
	{
		if (this->IsDistributed())
			return this->DistributedP.Diag.rows() == 0 && this->P.rows() == 0 ? "matrix free" : "distributed";
		if (this->SingleP.rows() > 0)
			return Utils::MatrixInfo(this->SingleP, "P");
		return this->P.rows() == 0 ? "matrix free" : Utils::MatrixInfo(this->P, "P");
	}

	virtual string RestrictionInfo()
	{
		if (this->IsDistributed())
			return this->DistributedP.Diag.rows() == 0 && this->R.rows() == 0 ? "matrix free" : "transpose of the prolongation";
		if (this->SingleR.rows() > 0)
			return Utils::MatrixInfo(this->SingleR, "R");
		return this->R.rows() == 0 ? "matrix free" : Utils::MatrixInfo(this->R, "R");
	}

	// The export of the prolongation and the restriction (-export mg)
	virtual void ExportTransfers()
	{
		this->ExportMatrix(P, "P");
		this->ExportMatrix(R, "R");
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

	virtual FloatVector Restrict(FloatVector& vectorOnThisLevel)
	{
		return SingleR * vectorOnThisLevel;
	}

	virtual Flops RestrictCost()
	{
		if (this->IsDistributed())
			return 2 * this->DistributedP.LocalNonZeros();
		return 2 * (R.nonZeros() + SingleR.nonZeros()); // one of them is empty
	}

	virtual Vector Prolong(Vector& vectorOnTheCoarserLevel)
	{
		if (this->IsDistributed())
			return this->DistributedP * vectorOnTheCoarserLevel;
		Vector vectorOnThisLevel = P * vectorOnTheCoarserLevel;
		return vectorOnThisLevel;
	}

	virtual FloatVector Prolong(FloatVector& vectorOnTheCoarserLevel)
	{
		return SingleP * vectorOnTheCoarserLevel;
	}

	virtual Flops ProlongCost()
	{
		if (this->IsDistributed())
			return 2 * this->DistributedP.LocalNonZeros();
		return 2 * (P.nonZeros() + SingleP.nonZeros()); // one of them is empty
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
		else if (this->SingleOperatorMatrix)
			return this->SingleOperatorMatrix->rows();
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