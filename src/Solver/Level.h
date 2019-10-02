#pragma once
#include <Eigen/Sparse>
#include "Smoother.h"
using namespace std;

class Level
{
public:
	int Number;
	SparseMatrix OperatorMatrix;
	bool UseGalerkinOperator = true;

	Smoother* PreSmoother = nullptr;
	Smoother* PostSmoother = nullptr;

	Level* FinerLevel = nullptr;
	Level* CoarserLevel = nullptr;

	BigNumber SetupComputationalWork = 0;

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
			cout << "\t\tFine grid operator: "; cout.flush();
			cout << Utils::MatrixInfo(this->OperatorMatrix, "A") << endl;
		}
		else
		{
			if (this->UseGalerkinOperator)
			{
				cout << "\t\tGalerkin operator : "; cout.flush();
				this->OperatorMatrix = FinerLevel->R * FinerLevel->OperatorMatrix * FinerLevel->P;
			}
			else
			{
				cout << "\t\tDiscretized operator: "; cout.flush();
				SetupDiscretizedOperator();
			}
			cout << Utils::MatrixInfo(this->OperatorMatrix, "A") << endl;
		}

		if (!this->IsCoarsestLevel())
		{
			cout << "\t\tProlongation      : "; cout.flush();
			SetupProlongation();
			cout << Utils::MatrixInfo(this->P, "P") << endl;

			cout << "\t\tRestriction       : "; cout.flush();
			SetupRestriction();
			cout << Utils::MatrixInfo(this->R, "R") << endl;
		}

		if (!this->IsCoarsestLevel())
		{
			cout << "\t\tSmoothers..." << endl;
			SetupSmoothers();
		}
	}

	Vector Restrict(Vector& vectorOnThisLevel)
	{
		Vector coarseVector = R * vectorOnThisLevel;
		return coarseVector;
	}

	double RestrictCost()
	{
		return 2 * R.nonZeros();
	}

	Vector Prolong(Vector& vectorOnTheCoarserLevel)
	{
		Vector vectorOnThisLevel = P * vectorOnTheCoarserLevel;
		return vectorOnThisLevel;
	}

	double ProlongCost()
	{
		return 2 * P.nonZeros();
	}

	virtual void ExportVector(Vector& v, string suffix) = 0;

	virtual ~Level()
	{
		delete PreSmoother;
		delete PostSmoother;
		if (CoarserLevel)
			delete CoarserLevel;
	}

protected:
	virtual void OnStartSetup() {}
	virtual void SetupDiscretizedOperator() {}
	virtual void SetupRestriction() = 0;
	virtual void SetupProlongation() = 0;

	void SetupSmoothers()
	{
		this->PreSmoother->Setup(this->OperatorMatrix);
		this->PostSmoother->Setup(this->OperatorMatrix);
	}
};