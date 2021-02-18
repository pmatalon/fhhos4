#pragma once
#include "IterativeSolver.h"
#include "SmootherFactory.h"
#include "Level.h"
using namespace std;

class Multigrid : public IterativeSolver
{
protected:
	Level* _fineLevel = nullptr;
private:
	bool _automaticNumberOfLevels;
	int _nLevels;
	NonZeroCoefficients _cycleSchema; // used only for drawing the cycle in the console
public:
	Solver* CoarseSolver = nullptr;
	char Cycle = 'V';
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
	FaceCoarseningStrategy FaceCoarseningStgy = FaceCoarseningStrategy::InterfaceCollapsing;
	int CoarseningFactor = 2;
	bool ExportComponents = false;
	bool DoNotCreateLevels = false;

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

protected:
	Multigrid() {};

public:
	int NumberOfLevels()
	{
		return _nLevels;
	}

	void Setup(const SparseMatrix& A) override
	{
		IterativeSolver::Setup(A);

		if (DoNotCreateLevels)
			return;

		cout << "Setup..." << endl;

		this->_fineLevel->OperatorMatrix = &A;
		this->_fineLevel->PreSmoother = this->_fineLevel->CreateSmoother(PreSmootherCode, PreSmoothingIterations, BlockSizeForBlockSmoothers, RelaxationParameter);
		this->_fineLevel->PostSmoother = this->_fineLevel->CreateSmoother(PostSmootherCode, PostSmoothingIterations, BlockSizeForBlockSmoothers, RelaxationParameter);
		this->_fineLevel->ExportComponents = ExportComponents;

		Level* currentLevel = this->_fineLevel;
		bool noCoarserMeshProvided = false;
		bool coarsestPossibleMeshReached = false;
		int levelNumber = 0;
		double operatorComplexity = 0;
		double fineNNZ = (double)this->_fineLevel->OperatorMatrix->nonZeros();
		double gridComplexity = 0;
		double fineNUnknowns = (double)this->_fineLevel->NUnknowns();
		while ((_automaticNumberOfLevels && currentLevel->NUnknowns() > MatrixMaxSizeForCoarsestLevel) || (levelNumber < _nLevels - 1))
		{
			// Can we coarsen the mesh?
			currentLevel->CoarsenMesh(this->CoarseningStgy, this->FaceCoarseningStgy, this->CoarseningFactor, noCoarserMeshProvided, coarsestPossibleMeshReached);
			if (noCoarserMeshProvided || coarsestPossibleMeshReached)
				break;

			// Build coarse level
			levelNumber++;
			Level* coarseLevel = CreateCoarseLevel(currentLevel);

			// Smoothers on coarse level
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
			coarseLevel->ExportComponents = ExportComponents;

			/*// FCG if K-cycle
			if (this->Cycle == 'K')
			{
				FlexibleConjugateGradient* fcg = new FlexibleConjugateGradient();
				coarseLevel->FCG = fcg;

				Multigrid* mg = new Multigrid();
				mg->_fineLevel = coarseLevel;
				mg->DoNotCreateLevels = true;
				mg->Cycle = this->Cycle;
				mg->WLoops = this->WLoops;
				fcg->Precond = Preconditioner(mg);

				fcg->ComputeExactSolution = false;
				fcg->PrintIterationResults = false;
				fcg->MaxIterations = 2;
				fcg->StoppingCrit = StoppingCriteria::NormalizedResidual;
				fcg->Tolerance = 0.25;
			}*/

			// Link between levels
			currentLevel->CoarserLevel = coarseLevel;
			coarseLevel->FinerLevel = currentLevel;

			// Setup fine level
			currentLevel->Setup();

			operatorComplexity += currentLevel->OperatorMatrix->nonZeros() / fineNNZ;
			gridComplexity += currentLevel->NUnknowns() / fineNUnknowns;

			currentLevel = coarseLevel;
		}

		_nLevels = levelNumber + 1;
		currentLevel->Setup();

		operatorComplexity += currentLevel->OperatorMatrix->nonZeros() / fineNNZ;
		gridComplexity += currentLevel->NUnknowns() / fineNUnknowns;

		if (coarsestPossibleMeshReached)
			Utils::Warning("impossible to coarsen the mesh any more.");
		if (noCoarserMeshProvided)
			Utils::Warning("cannot define a lower level because no coarser mesh has been provided.");

		this->SetupCoarseSolver();

		/*if (this->Cycle == 'K')
		{
			Level* level = this->_fineLevel->CoarserLevel;
			while (level)
			{
				Multigrid* mg = static_cast<Multigrid*>(level->FCG->Precond.GetSolver());
				mg->CoarseSolver = this->CoarseSolver;
				level = level->CoarserLevel;
			}
		}*/

		cout << endl;
		cout << "\tLevels             : " << _nLevels << endl;
		cout << "\tOperator complexity: " << operatorComplexity << endl;
		cout << "\tGrid complexity    : " << gridComplexity << endl;
		cout << endl;

		if (this->WLoops > 1)
			PrintCycleSchema();
	}

protected:
	virtual Level* CreateCoarseLevel(Level* fineLevel)
	{
		assert(false && "Not implemented. This method must be implemented in the subclass.");
	}

	void SetupCoarseSolver()
	{
		cout << "\t\tSetup coarse solver..." << endl;

		if (this->CoarseSolver == nullptr)
			this->CoarseSolver = new EigenSparseLU();

		Level* level = this->_fineLevel;
		while (level->CoarserLevel != nullptr)
			level = level->CoarserLevel;

		this->CoarseSolver->Setup(*level->OperatorMatrix);
	}

private:

	IterationResult ExecuteOneIteration(const Vector& b, Vector& x, const IterationResult& oldResult) override
	{
		IterationResult result(oldResult);
		MultigridCycle(this->_fineLevel, b, x, result);
		result.SetX(x);
		return result;
	}

	void MultigridCycle(Level* level, const Vector& b, Vector& x, IterationResult& result)
	{
		const SparseMatrix& A = *level->OperatorMatrix;

		if (Utils::ProgramArgs.Actions.ExportMultigridIterationVectors)
			level->ExportVector(b, "it" + to_string(this->IterationCount) + "_b");

		if (level->IsCoarsestLevel())
		{
			x = CoarseSolver->Solve(b);
			result.AddCost(CoarseSolver->SolvingComputationalWork);
		}
		else
		{
			if (Utils::ProgramArgs.Actions.ExportMultigridIterationVectors)
				level->ExportVector(x, "it" + to_string(this->IterationCount) + "_sol_beforePreSmoothing");

			//---------------//
			// Pre-smoothing //
			//---------------//

			//level->PreSmoother->Smooth(x, b);                                        result.AddCost(level->PreSmoother->SolvingComputationalWork());
			Vector r = level->PreSmoother->SmoothAndComputeResidual(x, b);           result.AddCost(level->PreSmoother->SolvingComputationalWork());

			if (Utils::ProgramArgs.Actions.ExportMultigridIterationVectors)
				level->ExportVector(x, "it" + to_string(this->IterationCount) + "_sol_afterPreSmoothing");

			//----------------------//
			// Residual computation //
			//----------------------//

			//Vector r = b - A * x;                                                    result.AddCost(Cost::DAXPY(A));
			
			//------------------------------------------------//
			// Restriction of the residual on the coarse grid //
			//------------------------------------------------//

			Vector rc = level->Restrict(r);                                          result.AddCost(level->RestrictCost());

			//--------------------------------------------------//
			// Residual equation Ae=r solved on the coarse grid //
			//--------------------------------------------------//

			Vector ec = Vector::Zero(rc.rows());
			bool ecEqualZero = true;
			if (this->Cycle == 'V' || this->Cycle == 'W')
			{
				for (int i = 0; i < this->WLoops; ++i)
				{
					MultigridCycle(level->CoarserLevel, rc, ec, result);
					if (level->CoarserLevel->IsCoarsestLevel())
						break;
				}
			}
			else if (this->Cycle == 'K')
			{
				if (level->CoarserLevel->IsCoarsestLevel())
					MultigridCycle(level->CoarserLevel, rc, ec, result); // exact solution
				else
				{
					//level->CoarserLevel->FCG->Solve(rc, ecEqualZero, ec);                    result.AddCost(level->CoarserLevel->FCG->SolvingComputationalWork);
					FCGForKCycle(level->CoarserLevel, rc, ec, result);
				}
			}

			//------------------------//
			// Coarse-grid correction //
			//------------------------//

			if (Utils::ProgramArgs.Actions.ExportMultigridIterationVectors)
			{
				level->ExportVector(ec, "it" + to_string(this->IterationCount) + "_ce");
				level->ExportVector(x, "it" + to_string(this->IterationCount) + "_sol");
				Vector cgc = level->Prolong(ec);
				level->ExportVector(cgc, "it" + to_string(this->IterationCount) + "_cgc");
			}

			x += level->Prolong(ec);                                              result.AddCost(level->ProlongCost() + Cost::AddVec(x));

			if (Utils::ProgramArgs.Actions.ExportMultigridIterationVectors)
				level->ExportVector(x, "it" + to_string(this->IterationCount) + "_sol_cgc");

			//----------------//
			// Post-smoothing //
			//----------------//

			level->PostSmoother->Smooth(x, b);                                   result.AddCost(level->PostSmoother->SolvingComputationalWork());

			if (Utils::ProgramArgs.Actions.ExportMultigridIterationVectors)
				level->ExportVector(x, "it" + to_string(this->IterationCount) + "_sol_afterPostSmoothing");
		}
	}

private:
	void FCGForKCycle(Level* level, Vector& r, Vector& x, IterationResult& result)
	{
		const SparseMatrix& A = *level->OperatorMatrix;
		double t = 0.25;

		Vector c = Vector::Zero(x.rows());
		MultigridCycle(level, r, c, result);
		Vector v = A * c;                                                    result.AddCost(Cost::MatVec(A));
		double rho1 = c.dot(v);                                              result.AddCost(Cost::Dot(c));
		double alpha1 = c.dot(r);                                            result.AddCost(Cost::Dot(c));

		double r_norm_old = r.norm();
		r -= alpha1 / rho1 * v;
		if (r.norm() <= t * r_norm_old)
			x = alpha1 / rho1 * c;
		else
		{
			Vector d = Vector::Zero(x.rows());
			MultigridCycle(level, r, d, result);
			Vector w = A * d;                                                result.AddCost(Cost::MatVec(A));
			double gamma = d.dot(v);                                         result.AddCost(Cost::Dot(d));
			double beta = d.dot(w);                                          result.AddCost(Cost::Dot(d));
			double alpha2 = d.dot(r);                                        result.AddCost(Cost::Dot(d));
			double rho2 = beta - gamma * gamma / rho1;
			x = (alpha1/rho1 - gamma*alpha2/(rho1*rho2))*c + alpha2 / rho2 * d;
		}
	}

protected:
	virtual void BeginSerialize(ostream& os) const
	{
		os << "Multigrid" << endl;
	}
	virtual void EndSerialize(ostream& os) const {}

public:
	virtual void Serialize(ostream& os) const override
	{
		BeginSerialize(cout);

		os << "\t" << "Cycle                   : ";
		if (this->Cycle == 'K')
			os << "K";
		else if (this->WLoops == 1)
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

		os << "\t" << "Levels                  : ";
		if (_automaticNumberOfLevels && _nLevels == 0)
			os << "automatic coarsening until matrix size <= " << MatrixMaxSizeForCoarsestLevel << endl;
		else if (_automaticNumberOfLevels && _nLevels > 0)
			os << _nLevels << " (automatic)" << endl;
		else
			os << _nLevels << endl;

		os << "\t" << "Coarsening strategy     : ";
		if (CoarseningStgy == CoarseningStrategy::StandardCoarsening)
			os << "standard [-cs s]" << endl;
		else if (CoarseningStgy == CoarseningStrategy::GMSHSplittingRefinement)
			os << "GMSH refinement by splitting from coarse mesh [-cs r]" << endl;
		else if (CoarseningStgy == CoarseningStrategy::BeyRefinement)
			os << "Bey's refinement from coarse mesh [-cs b]" << endl;
		else if (CoarseningStgy == CoarseningStrategy::IndependentRemeshing)
			os << "independant remeshing [-cs m]" << endl;
		else if (CoarseningStgy == CoarseningStrategy::DoublePairwiseAggregation)
			os << "double pairwise aggregation [-cs dpa]" << endl;
		//else if (CoarseningStgy == CoarseningStrategy::AgglomerationCoarsening)
			//os << "agglomeration [-cs a]" << endl;
		else if (CoarseningStgy == CoarseningStrategy::AgglomerationCoarseningByFaceNeighbours)
			os << "agglomeration by face neighbours [-cs n]" << endl;
		else if (CoarseningStgy == CoarseningStrategy::AgglomerationCoarseningByMostCoplanarFaces)
			os << "agglomeration by most coplanar faces [-cs mcf]" << endl;
		else if (CoarseningStgy == CoarseningStrategy::AgglomerationCoarseningByClosestCenter)
			os << "agglomeration by closest element center [-cs cc]" << endl;
		else if (CoarseningStgy == CoarseningStrategy::AgglomerationCoarseningByClosestFace)
			os << "agglomeration by closest face center [-cs clf]" << endl;
		else if (CoarseningStgy == CoarseningStrategy::AgglomerationCoarseningByLargestInterface)
			os << "agglomeration by largest interface [-cs li]" << endl;
		else if (CoarseningStgy == CoarseningStrategy::FaceCoarsening)
			os << "face coarsening [-cs f]" << endl;
		else
			os << "unknown" << endl;

		os << "\t" << "Face coarsening strategy: ";
		if (Utils::IsRefinementStrategy(CoarseningStgy) || CoarseningStgy == CoarseningStrategy::IndependentRemeshing)
			os << "NA" << endl;
		else if (FaceCoarseningStgy == FaceCoarseningStrategy::None)
			os << "none" << endl;
		else if (FaceCoarseningStgy == FaceCoarseningStrategy::InterfaceCollapsing)
			os << "interface collapsing [-fcs c]" << endl;
		else if (FaceCoarseningStgy == FaceCoarseningStrategy::InterfaceCollapsingAndTryAggregInteriorToInterfaces)
			os << "interface collapsing and try aggregate interior faces too [-fcs i]" << endl;
		else
			os << "unknown" << endl;

		Smoother* preSmoother = this->_fineLevel->CreateSmoother(PreSmootherCode, PreSmoothingIterations, BlockSizeForBlockSmoothers, RelaxationParameter);
		Smoother* postSmoother = this->_fineLevel->CreateSmoother(PostSmootherCode, PostSmoothingIterations, BlockSizeForBlockSmoothers, RelaxationParameter);
		os << "\t" << "Pre-smoothing           : " << *preSmoother << endl;
		os << "\t" << "Post-smoothing          : " << *postSmoother;
		delete preSmoother;
		delete postSmoother;
		if (PreSmootherCode.compare("bj") == 0 || PostSmootherCode.compare("bj") == 0)
		{
			os << endl;
			Utils::Warning("Note that without relaxation parameter, the (block) Jacobi iteration is not a smoother. You can use 'bj23' (under-relaxation parameter = 2/3) instead.");
		}
		if (this->CoarseSolver)
		{
			os << endl;
			os << "\t" << "Coarse solver           : " << (*this->CoarseSolver);
		}

		EndSerialize(cout);
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
		delete CoarseSolver;
	}
};