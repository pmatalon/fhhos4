#pragma once
#include <Eigen/Sparse>
#include "Smoother.h"
using namespace std;

class Level
{
public:
	int Number;

	Eigen::SparseMatrix<double> OperatorMatrix;

	Smoother PreSmoother;
	Smoother PostSmoother;

	Level* FinerLevel = NULL;
	Level* CoarserLevel = NULL;

	Level(int number, Smoother preSmoother, Smoother postSmoother) : PreSmoother(preSmoother), PostSmoother(postSmoother)
	{
		this->Number = number;
	}

	bool IsFinestLevel()
	{
		return this->FinerLevel == NULL;
	}

	bool IsCoarsestLevel()
	{
		return this->CoarserLevel == NULL;
	}

	virtual void Setup(const Eigen::SparseMatrix<double>& A)
	{
		SetupRestriction();
		SetupProlongation();

		SetupOperator(A);

		//---- Can be done async ----//
		if (!this->IsCoarsestLevel())
			this->CoarserLevel->Setup(this->OperatorMatrix);
		//---------------------------//

		SetupSmoothers();
	}

	virtual Eigen::VectorXd Restrict(Eigen::VectorXd& vectorOnThisLevel) = 0;
	virtual Eigen::VectorXd Prolong(Eigen::VectorXd& vectorOnTheCoarserLevel) = 0;

private:
	virtual void SetupRestriction() = 0;
	virtual void SetupProlongation() = 0;
	virtual void SetupOperator(const Eigen::SparseMatrix<double>& A) = 0;

	void SetupSmoothers()
	{
		this->PreSmoother.Setup(this->OperatorMatrix);
		this->PostSmoother.Setup(this->OperatorMatrix);
	}
};