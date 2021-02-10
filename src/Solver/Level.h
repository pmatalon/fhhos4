#pragma once
#include "Smoother.h"
#include "FlexibleConjugateGradient.h"
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

	bool ExportComponents = false;
	BigNumber SetupComputationalWork = 0;

	//FlexibleConjugateGradient* FCG = nullptr; // used in K-cycle

protected:
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

	void Setup()
	{
		cout << "\tSetup level " << this->Number << endl;

		OnStartSetup();

		if (this->IsFinestLevel())
		{
			cout << "\t\tFine grid operator  : "; cout.flush();
			cout << Utils::MatrixInfo(*this->OperatorMatrix, "A") << endl;
		}
		else
		{
			if (this->UseGalerkinOperator)
			{
				cout << "\t\tGalerkin operator   : "; cout.flush();
				if (!this->OperatorMatrix)
					ComputeGalerkinOperator();
			}
			else
			{
				cout << "\t\tDiscretized operator: "; cout.flush();
				SetupDiscretizedOperator();
			}
			cout << Utils::MatrixInfo(*this->OperatorMatrix, "A") << endl;
		}

		const SparseMatrix &A = *this->OperatorMatrix;

		if (!this->IsCoarsestLevel())
		{
			cout << "\t\tProlongation        : "; cout.flush();
			SetupProlongation();
			cout << Utils::MatrixInfo(this->P, "P") << endl;

			cout << "\t\tRestriction         : "; cout.flush();
			SetupRestriction();
			cout << Utils::MatrixInfo(this->R, "R") << endl;

			cout << "\t\tPreSmoothing        : "; cout.flush();
			PreSmoother->Setup(A);
			cout << PreSmoother->Iterations() << " iteration" << (PreSmoother->Iterations() > 1 ? "s" : "") <<  endl;

			cout << "\t\tPostSmoothing       : "; cout.flush();
			PostSmoother->Setup(A);
			cout << PostSmoother->Iterations() << " iteration" << (PostSmoother->Iterations() > 1 ? "s" : "") << endl;
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

	virtual double RestrictCost()
	{
		return 2 * R.nonZeros();
	}

	virtual Vector Prolong(Vector& vectorOnTheCoarserLevel)
	{
		Vector vectorOnThisLevel = P * vectorOnTheCoarserLevel;
		return vectorOnThisLevel;
	}

	virtual double ProlongCost()
	{
		return 2 * P.nonZeros();
	}

	virtual BigNumber NUnknowns()
	{
		if (this->OperatorMatrix)
			return this->OperatorMatrix->rows();
		else if (!this->IsFinestLevel())
			return this->FinerLevel->P.cols();
		assert(false && "NUnknowns() should be overriden");
	}

	virtual void CoarsenMesh(CoarseningStrategy coarseningStgy, int coarseningFactor, bool& noCoarserMeshProvided, bool& coarsestPossibleMeshReached) = 0;

	virtual void ExportVector(const Vector& v, string suffix, int levelNumber) = 0;
	virtual void ExportMatrix(const SparseMatrix& M, string suffix, int levelNumber) = 0;

	void ExportVector(const Vector& v, string suffix)
	{
		this->ExportVector(v, suffix, this->Number);
	}
	void ExportMatrix(const SparseMatrix& M, string suffix)
	{
		this->ExportMatrix(M, suffix, this->Number);
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