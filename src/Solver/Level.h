#pragma once
#include "Smoother.h"
using namespace std;

class Level
{
public:
	int Number;
	SparseMatrix OperatorMatrix;
	bool UseGalerkinOperator = false;

	Smoother* PreSmoother = nullptr;
	Smoother* PostSmoother = nullptr;

	Level* FinerLevel = nullptr;
	Level* CoarserLevel = nullptr;

	bool ExportMatrices = false;
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
			cout << "\t\tFine grid operator  : "; cout.flush();
			cout << Utils::MatrixInfo(this->OperatorMatrix, "A") << endl;
		}
		else
		{
			if (this->UseGalerkinOperator)
			{
				cout << "\t\tGalerkin operator   : "; cout.flush();
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
			cout << "\t\tProlongation        : "; cout.flush();
			SetupProlongation();
			cout << Utils::MatrixInfo(this->P, "P") << endl;

			cout << "\t\tRestriction         : "; cout.flush();
			SetupRestriction();
			cout << Utils::MatrixInfo(this->R, "R") << endl;

			cout << "\t\tPreSmoothing        : "; cout.flush();
			PreSmoother->Setup(this->OperatorMatrix);
			cout << PreSmoother->Iterations() << " iteration" << (PreSmoother->Iterations() > 1 ? "s" : "") <<  endl;

			cout << "\t\tPostSmoothing       : "; cout.flush();
			PostSmoother->Setup(this->OperatorMatrix);
			cout << PostSmoother->Iterations() << " iteration" << (PostSmoother->Iterations() > 1 ? "s" : "") << endl;
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

	virtual BigNumber NUnknowns() = 0;
	virtual void CoarsenMesh(CoarseningStrategy coarseningStgy, bool& noCoarserMeshProvided, bool& coarsestPossibleMeshReached) = 0;
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
};