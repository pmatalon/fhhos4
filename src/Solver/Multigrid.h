#pragma once
#include "IterativeSolver.h"
#include "Level.h"
using namespace std;

class Multigrid : public IterativeSolver
{
protected:
	Level* _fineLevel = NULL;
	Solver* _coarseSolver = NULL;
private:
	NonZeroCoefficients _cycleSchema;
public:
	int WLoops = 1; // 1 --> V-cycle, >1 --> W-cycle 

	Multigrid() : IterativeSolver() { }

	int NumberOfLevels()
	{
		int nLevels = 1;
		Level* level = this->_fineLevel;
		while (level = level->CoarserLevel)
			nLevels++;
		return nLevels;
	}

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

	IterationResult ExecuteOneIteration(const Eigen::VectorXd& b, Eigen::VectorXd& x, const IterationResult& oldResult) override
	{
		IterationResult result(oldResult);
		x = MultigridCycle(this->_fineLevel, b, x, result);
		result.SetX(x);
		return result;
	}

	Eigen::VectorXd MultigridCycle(Level* level, const Eigen::VectorXd& b, Eigen::VectorXd& initialGuess, IterationResult& result)
	{
		SparseMatrix A = level->OperatorMatrix;
		Eigen::VectorXd x;

		if (level->IsCoarsestLevel())
		{
			x = _coarseSolver->Solve(b);
			result.AddCost(_coarseSolver->SolvingComputationalWork);
		}
		else
		{
			x = initialGuess;

			//level->ExportVector(x, "mg_initialGuess");

			// Pre-smoothing //
			x = level->PreSmoother->Smooth(x, b);
			result.AddCost(level->PreSmoother->SolvingComputationalWork());

			//level->ExportVector(x, "mg_afterPreSmoothing");

			// Residual computation //
			Eigen::VectorXd r = b - A * x;
			result.AddCost(2 * A.nonZeros());

			//cout << "res = " << (r.norm() / b.norm()) << endl;
			
			// Restriction of the residual on the coarse grid //
			Eigen::VectorXd rc = level->Restrict(r);
			result.AddCost(level->RestrictCost());

			// Residual equation Ae=r solved on the coarse grid //
			Eigen::VectorXd ec = Eigen::VectorXd::Zero(rc.rows());
			for (int i = 0; i < this->WLoops; ++i)
			{
				ec = MultigridCycle(level->CoarserLevel, rc, ec, result);
				if (level->CoarserLevel->IsCoarsestLevel())
					break;
			}

			// Coarse-grid correction //
			x = x + level->Prolong(ec);
			result.AddCost(level->ProlongCost());

			// Post-smoothing //
			x = level->PostSmoother->Smooth(x, b);
			result.AddCost(level->PostSmoother->SolvingComputationalWork());
			//level->ExportVector(x, "mg_afterPostSmoothing");
		}

		return x;
	}

public:
	void PrintCycleSchema()
	{
		StoreCycleSchema(this->_fineLevel);
		int nLevels = this->NumberOfLevels();
		SparseMatrix schema2(nLevels, 50);
		_cycleSchema.Fill(schema2);
		Eigen::MatrixXd schema = schema2;
		for (int row = 0; row < schema.rows(); row++)
		{
			for (int col = 0; col < schema.cols(); col++)
			{
				if (schema(row, col) == 0)
					cout << "  ";
				else
					cout << " o";
			}
			cout << endl;
		}
		cout << endl;
	}

private:
	inline void PrintLevel(Level* level)
	{
		for (int i = 0; i < level->Number; i++)
			cout << "   ";
		cout << "o" << endl;
	}
	inline void StoreLevelInSchema(Level* level)
	{
		static int col = 0;
		_cycleSchema.Add(level->Number, col++, 1);
	}
	void StoreCycleSchema(Level* level)
	{
		static Level* currentLevel = nullptr;
		if (currentLevel != level)
		{
			StoreLevelInSchema(level);
			currentLevel = level;
		}

		if (!level->IsCoarsestLevel())
		{
			for (int i = 0; i < this->WLoops; ++i)
			{
				StoreCycleSchema(level->CoarserLevel);
				if (level->CoarserLevel->IsCoarsestLevel())
					break;
			}

			if (currentLevel != level)
			{
				StoreLevelInSchema(level);
				currentLevel = level;
			}
		}
	}

public:
	virtual ~Multigrid()
	{
		delete _fineLevel;
		delete _coarseSolver;
	}
};