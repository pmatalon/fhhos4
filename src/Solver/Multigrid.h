#pragma once
#include "IterativeSolver.h"
#include "SmootherFactory.h"
#include "Level.h"
using namespace std;

class Multigrid : public IterativeSolver
{
protected:
	Level* _fineLevel = nullptr;
	Solver* _coarseSolver = nullptr;
private:
	bool _automaticNumberOfLevels;
	int _nLevels;
	NonZeroCoefficients _cycleSchema;
public:
	int WLoops = 1; // 1 --> V-cycle, >= 2 --> W-cycle 
	int MatrixMaxSizeForCoarsestLevel;
	bool UseGalerkinOperator = false;
	string PreSmootherCode = "bgs";
	string PostSmootherCode = "rbgs";
	int PreSmoothingIterations = 1;
	int PostSmoothingIterations = 1;
	int CoarseLevelChangeSmoothingCoeff = 0;
	char CoarseLevelChangeSmoothingOperator = '+';
	int BlockSizeForBlockSmoothers = -1;
	double RelaxationParameter = 1;
	CoarseningStrategy CoarseningStgy = CoarseningStrategy::StandardCoarsening;
	bool ExportMatrices = false;

	Multigrid(int nLevels) : IterativeSolver()
	{
		this->_automaticNumberOfLevels = nLevels <= 0;

		if (this->_automaticNumberOfLevels)
		{
			this->_nLevels = 0;
			this->MatrixMaxSizeForCoarsestLevel = 1000;
		}
		else
		{
			this->_nLevels = nLevels;
			this->MatrixMaxSizeForCoarsestLevel = 0;
		}
	}

	int NumberOfLevels()
	{
		return _nLevels;
	}

	void Setup(const SparseMatrix& A) override
	{
		IterativeSolver::Setup(A);

		cout << "Setup..." << endl;

		this->_fineLevel->OperatorMatrix = A;
		this->_fineLevel->PreSmoother = this->_fineLevel->CreateSmoother(PreSmootherCode, PreSmoothingIterations, BlockSizeForBlockSmoothers, RelaxationParameter);
		this->_fineLevel->PostSmoother = this->_fineLevel->CreateSmoother(PostSmootherCode, PostSmoothingIterations, BlockSizeForBlockSmoothers, RelaxationParameter);
		this->_fineLevel->UseGalerkinOperator = false;
		this->_fineLevel->ExportMatrices = ExportMatrices;

		Level* currentLevel = this->_fineLevel;
		bool noCoarserMeshProvided = false;
		bool coarsestPossibleMeshReached = false;
		int levelNumber = 0;
		while ((_automaticNumberOfLevels && currentLevel->NUnknowns() > MatrixMaxSizeForCoarsestLevel) || (levelNumber < _nLevels - 1))
		{
			// Can we coarsen the mesh?
			currentLevel->CoarsenMesh(this->CoarseningStgy, noCoarserMeshProvided, coarsestPossibleMeshReached);
			if (noCoarserMeshProvided || coarsestPossibleMeshReached)
				break;

			// Build coarse level
			levelNumber++;
			Level* coarseLevel = CreateCoarseLevel(currentLevel);

			int preSmoothingIterations = currentLevel->PreSmoother->Iterations();
			int postSmoothingIterations = currentLevel->PostSmoother->Iterations();
			if (CoarseLevelChangeSmoothingOperator == '+')
			{
				preSmoothingIterations += CoarseLevelChangeSmoothingCoeff;
				postSmoothingIterations += CoarseLevelChangeSmoothingCoeff;
			}
			else if (CoarseLevelChangeSmoothingOperator == '-')
			{
				preSmoothingIterations -= CoarseLevelChangeSmoothingCoeff;
				postSmoothingIterations -= CoarseLevelChangeSmoothingCoeff;
			}
			else if (CoarseLevelChangeSmoothingOperator == '*')
			{
				preSmoothingIterations *= CoarseLevelChangeSmoothingCoeff;
				postSmoothingIterations *= CoarseLevelChangeSmoothingCoeff;
			}
			else
				assert(false);

			coarseLevel->PreSmoother = coarseLevel->CreateSmoother(PreSmootherCode, preSmoothingIterations, BlockSizeForBlockSmoothers, RelaxationParameter);
			coarseLevel->PostSmoother = coarseLevel->CreateSmoother(PostSmootherCode, postSmoothingIterations, BlockSizeForBlockSmoothers, RelaxationParameter);

			coarseLevel->UseGalerkinOperator = UseGalerkinOperator;
			coarseLevel->ExportMatrices = ExportMatrices;

			// Link between levels
			currentLevel->CoarserLevel = coarseLevel;
			coarseLevel->FinerLevel = currentLevel;

			// Setup fine level
			currentLevel->Setup();

			currentLevel = coarseLevel;
		}

		_nLevels = levelNumber + 1;
		currentLevel->Setup();

		if (coarsestPossibleMeshReached)
			cout << Utils::BeginYellow << "Warning: impossible to coarsen the mesh any more." << Utils::EndColor;
		if (noCoarserMeshProvided)
			cout << Utils::BeginYellow << "Warning: cannot build coarser level because no coarser mesh has been provided." << Utils::EndColor;

		this->SetupCoarseSolver();
		cout << "\t--> " << _nLevels << " levels built." << endl;

		if (this->WLoops > 1)
			PrintCycleSchema();
	}

protected:
	virtual Level* CreateCoarseLevel(Level* fineLevel) = 0;

	void SetupCoarseSolver()
	{
		cout << "\t\tSetup coarse solver..." << endl;

		if (this->_coarseSolver == nullptr)
			this->_coarseSolver = new EigenSparseLU();

		Level* level = this->_fineLevel;
		while (level->CoarserLevel != nullptr)
			level = level->CoarserLevel;

		this->_coarseSolver->Setup(level->OperatorMatrix);
	}

private:

	IterationResult ExecuteOneIteration(const Vector& b, Vector& x, const IterationResult& oldResult) override
	{
		IterationResult result(oldResult);
		x = MultigridCycle(this->_fineLevel, b, x, result);
		result.SetX(x);
		return result;
	}

	Vector MultigridCycle(Level* level, const Vector& b, Vector& initialGuess, IterationResult& result)
	{
		SparseMatrix A = level->OperatorMatrix;
		Vector x;

		if (level->IsCoarsestLevel())
		{
			x = _coarseSolver->Solve(b);
			//result.AddCost(_coarseSolver->SolvingComputationalWork);
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
			Vector r = b - A * x;
			result.AddCost(2 * A.nonZeros());

			//cout << "res = " << (r.norm() / b.norm()) << endl;
			
			// Restriction of the residual on the coarse grid //
			Vector rc = level->Restrict(r);
			result.AddCost(level->RestrictCost());

			// Residual equation Ae=r solved on the coarse grid //
			Vector ec = Vector::Zero(rc.rows());
			for (int i = 0; i < this->WLoops; ++i)
			{
				ec = MultigridCycle(level->CoarserLevel, rc, ec, result);
				if (level->CoarserLevel->IsCoarsestLevel())
					break;
			}

			// Coarse-grid correction //
			/*if (this->ExportMatrices && this->IterationCount == 1)
			{
				level->ExportVector(ec, "ce");
				level->ExportVector(x, "sol");
				auto cgc = level->Prolong(ec);
				level->ExportVector(cgc, "cgc");
			}*/

			x = x + level->Prolong(ec);

			/*if (this->ExportMatrices && this->IterationCount == 1)
				level->ExportVector(x, "sol_cgc");*/

			result.AddCost(level->ProlongCost());

			// Post-smoothing //
			//if (this->ExportMatrices && !level->IsCoarsestLevel() && level->CoarserLevel->IsCoarsestLevel())
			x = level->PostSmoother->Smooth(x, b);
			result.AddCost(level->PostSmoother->SolvingComputationalWork());
			//if (this->ExportMatrices && this->IterationCount == 1)
				//level->ExportVector(x, "sol_afterPostSmoothing");
		}

		return x;
	}

protected:
	virtual void BeginSerialize(ostream& os) const
	{
		os << "Multigrid" << endl;
	}

public:
	virtual void Serialize(ostream& os) const override
	{
		BeginSerialize(cout);

		os << "\t" << "Cycle              : ";
		if (this->WLoops == 1)
			os << "V";
		else if (this->WLoops == 2)
			os << "W";
		else
			os << "W(" << this->WLoops << " loops)";
		os << "(" << PreSmoothingIterations << "," << PostSmoothingIterations;
		if ((CoarseLevelChangeSmoothingOperator == '+' && CoarseLevelChangeSmoothingCoeff > 0) ||
			(CoarseLevelChangeSmoothingOperator == '-' && CoarseLevelChangeSmoothingCoeff > 0) ||
			(CoarseLevelChangeSmoothingOperator == '*' && CoarseLevelChangeSmoothingCoeff > 1))
			os << "," << CoarseLevelChangeSmoothingOperator << CoarseLevelChangeSmoothingCoeff;
		os << ")" << endl;

		os << "\t" << "Levels             : ";
		if (_automaticNumberOfLevels && _nLevels == 0)
			os << "automatic coarsening until matrix size <= " << MatrixMaxSizeForCoarsestLevel << endl;
		else if (_automaticNumberOfLevels && _nLevels > 0)
			os << _nLevels << " (automatic)" << endl;
		else
			os << _nLevels << endl;

		os << "\t" << "Coarsening strategy: ";
		if (CoarseningStgy == CoarseningStrategy::StandardCoarsening)
			os << "standard" << endl;
		else if (CoarseningStgy == CoarseningStrategy::AgglomerationCoarsening)
			os << "agglomeration" << endl;
		else if (CoarseningStgy == CoarseningStrategy::SplittingRefinement)
			os << "refinement by splitting from coarse mesh" << endl;
		else if (CoarseningStgy == CoarseningStrategy::BeyRefinement)
			os << "Bey's refinement from coarse mesh" << endl;
		else if (CoarseningStgy == CoarseningStrategy::FaceCoarsening)
			os << "face coarsening" << endl;
		else
			os << "unknown" << endl;

		Smoother* preSmoother = this->_fineLevel->CreateSmoother(PreSmootherCode, PreSmoothingIterations, BlockSizeForBlockSmoothers, RelaxationParameter);
		Smoother* postSmoother = this->_fineLevel->CreateSmoother(PostSmootherCode, PostSmoothingIterations, BlockSizeForBlockSmoothers, RelaxationParameter);
		os << "\t" << "Pre-smoothing      : " << *preSmoother << endl;
		os << "\t" << "Post-smoothing     : " << *postSmoother;
		delete preSmoother;
		delete postSmoother;
	}

	void PrintCycleSchema()
	{
		StoreCycleSchema(this->_fineLevel);
		int nLevels = this->NumberOfLevels();
		SparseMatrix schema2(nLevels, 50);
		_cycleSchema.Fill(schema2);
		DenseMatrix schema = schema2;
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