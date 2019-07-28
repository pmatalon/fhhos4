#pragma once
#include <Eigen/Sparse>
#include "IterativeSolver.h"
#include "Level.h"
#include "EigenSparseLU.h"
using namespace std;

class Multigrid : public IterativeSolver
{
protected:
	Level* _fineLevel = NULL;
	Solver* _coarseSolver = NULL;

public:
	Multigrid() : IterativeSolver() { }

protected:
	void SetupCoarseSolver()
	{
		cout << "\t\tSetup coarse solver..." << endl;

		if (this->_coarseSolver == NULL)
			this->_coarseSolver = new EigenSparseLU();

		Level* level = this->_fineLevel;
		while (level->CoarserLevel != NULL)
			level = level->CoarserLevel;

		this->_coarseSolver->Setup(level->OperatorMatrix);
	}

private:

	Eigen::VectorXd ExecuteOneIteration(const Eigen::VectorXd& b, Eigen::VectorXd& x) override
	{
		x = MultigridCycle(this->_fineLevel, b, x);
		return x;
	}

	Eigen::VectorXd MultigridCycle(Level* level, const Eigen::VectorXd& b, Eigen::VectorXd& initialGuess)
	{
		SparseMatrix A = level->OperatorMatrix;
		Eigen::VectorXd x;

		if (level->IsCoarsestLevel())
			x = _coarseSolver->Solve(b);
		else
		{
			x = initialGuess;

			//level->ExportVector(x, "mg_initialGuess");

			// Pre-smoothing //
			x = level->PreSmoother->Smooth(x, b);

			//level->ExportVector(x, "mg_afterPreSmoothing");

			// Residual computation //
			Eigen::VectorXd r = b - A * x;

			//cout << "res = " << (r.norm() / b.norm()) << endl;
			
			// Restriction of the residual on the coarse grid //
			Eigen::VectorXd rc = level->Restrict(r);

			// Residual equation Ae=r solved on the coarse grid //
			Eigen::VectorXd zero = Eigen::VectorXd::Zero(rc.rows());
			Eigen::VectorXd ec = MultigridCycle(level->CoarserLevel, rc, zero);

			// Coarse-grid correction //
			x = x + level->Prolong(ec);

			// Post-smoothing //
			x = level->PostSmoother->Smooth(x, b);
			//level->ExportVector(x, "mg_afterPostSmoothing");
		}

		return x;
	}

public:
	virtual ~Multigrid()
	{
		delete _fineLevel;
		delete _coarseSolver;
	}
};