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
	MGType MultigridType = MGType::h_Multigrid;
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
	double CoarseningFactor = 2;
	int CoarsePolyDegree = 1; // used in p-Multigrid
	bool ExportComponents = false;
	bool DoNotCreateLevels = false;

	Timer SmoothingAndResTimer;
	Timer IntergridTransferTimer;
	Timer CoarseSolverTimer;

	Multigrid(MGType multigridType, int nLevels) : IterativeSolver()
	{
		this->MultigridType = multigridType;
		if (multigridType == MGType::h_Multigrid)
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
		else
		{
			this->_automaticNumberOfLevels = true;
			this->_nLevels = 0;
			this->MatrixMaxSizeForCoarsestLevel = 0;
		}
	}

	bool CanOptimizeResidualComputation() override
	{
		return this->_nLevels > 1 && this->_fineLevel->PostSmoother->CanOptimizeResidualComputation();
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

		if (!this->_fineLevel)
			this->_fineLevel = this->CreateFineLevel();

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

		while (CoarseLevelNeeded(currentLevel, levelNumber))
		{
			if (MultigridType == MGType::h_Multigrid)
			{
				// Can we coarsen the mesh?
				currentLevel->CoarsenMesh(this->CoarseningStgy, this->FaceCoarseningStgy, this->CoarseningFactor, noCoarserMeshProvided, coarsestPossibleMeshReached);
				if (noCoarserMeshProvided || coarsestPossibleMeshReached)
					break;
			}

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


			int blockSize = MultigridType == MGType::h_Multigrid ? BlockSizeForBlockSmoothers : coarseLevel->BlockSizeForBlockSmoothers();
			coarseLevel->PreSmoother  = coarseLevel->CreateSmoother( PreSmootherCode,  preSmoothingIterations, blockSize, RelaxationParameter);
			coarseLevel->PostSmoother = coarseLevel->CreateSmoother(PostSmootherCode, postSmoothingIterations, blockSize, RelaxationParameter);

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

		this->InitializeCoarseSolver();
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
	virtual Level* CreateFineLevel() const
	{
		assert(false && "Not implemented. This method must be implemented in the subclass.");
	}
	virtual Level* CreateCoarseLevel(Level* fineLevel)
	{
		assert(false && "Not implemented. This method must be implemented in the subclass.");
	}

	virtual void InitializeCoarseSolver()
	{
		if (!this->CoarseSolver)
			this->CoarseSolver = new EigenSparseLU();
	}

	void SetupCoarseSolver()
	{
		cout << "\t\tSetup coarse solver..." << endl;

		Level* level = this->_fineLevel;
		while (level->CoarserLevel != nullptr)
			level = level->CoarserLevel;

		this->CoarseSolver->Setup(*level->OperatorMatrix);
	}

private:

	bool CoarseLevelNeeded(Level* currentLevel, int levelNumber)
	{
		if (MultigridType == MGType::h_Multigrid)
			return (_automaticNumberOfLevels && currentLevel->NUnknowns() > MatrixMaxSizeForCoarsestLevel) || (levelNumber < _nLevels - 1);
		else if (MultigridType == MGType::p_Multigrid)
			return currentLevel->PolynomialDegree() > CoarsePolyDegree;
		assert(false);
	}

	IterationResult ExecuteOneIteration(const Vector& b, Vector& x, bool& xEquals0, bool computeResidual, bool computeAx, const IterationResult& oldResult) override
	{
		IterationResult result(oldResult);

		MultigridCycle(this->_fineLevel, b, x, xEquals0, computeResidual, computeAx, result);
		result.SetX(x);
		return result;
	}

	void MultigridCycle(Level* level, const Vector& b, Vector& x, bool& xEquals0, bool computeResidual, bool computeAx, IterationResult& result)
	{
		const SparseMatrix& A = *level->OperatorMatrix;

		if (Utils::ProgramArgs.Actions.ExportMultigridIterationVectors)
			level->ExportVector(b, "it" + to_string(this->IterationCount) + "_b");

		if (level->IsCoarsestLevel())
		{
			                                                                         CoarseSolverTimer.Start();
			x = CoarseSolver->Solve(b);                                              result.AddCost(CoarseSolver->SolvingComputationalWork);
			                                                                         CoarseSolverTimer.Pause();
			xEquals0 = false;
		}
		else
		{
			if (Utils::ProgramArgs.Actions.ExportMultigridIterationVectors)
				level->ExportVector(x, "it" + to_string(this->IterationCount) + "_sol_beforePreSmoothing");

			//---------------//
			// Pre-smoothing //
			//---------------//

			auto flopBeforePreSmooth = result._iterationComputationalWork;
			SmoothingAndResTimer.Start();

			Vector r;
			if (level->PreSmoother->CanOptimizeResidualComputation())
			{
				level->PreSmoother->SmoothAndComputeResidualOrAx(x, b, xEquals0, true, false);    result.AddCost(level->PreSmoother->SolvingComputationalWork());
				r = level->PreSmoother->Residual();
			}
			else
			{
				level->PreSmoother->Smooth(x, b, xEquals0);                               result.AddCost(level->PreSmoother->SolvingComputationalWork());
				// Residual computation
				r = b - A * x;                                                            result.AddCost(Cost::DAXPY(A));
			}

			SmoothingAndResTimer.Stop();

			if (Utils::ProgramArgs.Actions.ExportMultigridIterationVectors)
				level->ExportVector(x, "it" + to_string(this->IterationCount) + "_sol_afterPreSmoothing");

			auto flopPreSmoothAndRes = result._iterationComputationalWork - flopBeforePreSmooth;
			//cout << "- flopPreSmoothAndRes = " << flopPreSmoothAndRes << endl;

			//------------------------------------------------//
			// Restriction of the residual on the coarse grid //
			//------------------------------------------------//

			auto flopBeforeRestrict = result._iterationComputationalWork;
			IntergridTransferTimer.Start();

			Vector rc = level->Restrict(r);                                               result.AddCost(level->RestrictCost());

			IntergridTransferTimer.Pause();
			auto flopRestrict = result._iterationComputationalWork - flopBeforeRestrict;
			//cout << "- flopRestrict = " << flopRestrict << endl;

			//--------------------------------------------------//
			// Residual equation Ae=r solved on the coarse grid //
			//--------------------------------------------------//

			Vector ec;
			if (this->Cycle == 'V' || this->Cycle == 'W')
			{
				ec = Vector::Zero(rc.rows());
				bool ecEquals0 = true;
				for (int i = 0; i < this->WLoops; ++i)
				{
					IterationResult coarseResult;
					MultigridCycle(level->CoarserLevel, rc, ec, ecEquals0, false, false, coarseResult);   result.AddCost(coarseResult.SolvingComputationalWork());
					if (level->CoarserLevel->IsCoarsestLevel())
						break;
				}
			}
			else if (this->Cycle == 'K')
			{
				IterationResult coarseResult;
				if (level->CoarserLevel->IsCoarsestLevel())
				{
					bool ecEquals0 = true;
					MultigridCycle(level->CoarserLevel, rc, ec, ecEquals0, false, false, coarseResult); // exact solution
					                                                                               result.AddCost(coarseResult.SolvingComputationalWork());
				}
				else
				{
					//level->CoarserLevel->FCG->Solve(rc, ecEqualZero, ec);                    result.AddCost(level->CoarserLevel->FCG->SolvingComputationalWork);
					ec = FCGForKCycle(level->CoarserLevel, rc, coarseResult);                      result.AddCost(coarseResult.SolvingComputationalWork());
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


			IntergridTransferTimer.Start();

			x += level->Prolong(ec);                                                  result.AddCost(level->ProlongCost() + Cost::AddVec(x));

			IntergridTransferTimer.Pause();

			if (Utils::ProgramArgs.Actions.ExportMultigridIterationVectors)
				level->ExportVector(x, "it" + to_string(this->IterationCount) + "_sol_cgc");

			//----------------//
			// Post-smoothing //
			//----------------//

			assert(!xEquals0);

			auto flopBeforePostSmooth = result._iterationComputationalWork;
			SmoothingAndResTimer.Start();

			if (level->PostSmoother->CanOptimizeResidualComputation() && (computeResidual || computeAx))
			{
				level->PostSmoother->SmoothAndComputeResidualOrAx(x, b, xEquals0, computeResidual, computeAx);    result.AddCost(level->PostSmoother->SolvingComputationalWork());
				if (computeResidual)
					result.Residual = level->PostSmoother->Residual();
				if (computeAx)
					result.Ax = level->PostSmoother->Ax();
			}
			else
			{
				level->PostSmoother->Smooth(x, b, xEquals0);                          result.AddCost(level->PostSmoother->SolvingComputationalWork());
				if (computeAx)
				{
					result.Ax = A * x;                                                result.AddCost(Cost::MatVec(A));
					if (computeResidual)
					{
						result.Residual = b - result.Ax;                              result.AddCost(Cost::AddVec(b));
					}
				}
				else if (computeResidual)
				{
					result.Residual = b - A*x;                                        result.AddCost(Cost::DAXPY(A));
				}
			}

			SmoothingAndResTimer.Pause();

			if (Utils::ProgramArgs.Actions.ExportMultigridIterationVectors)
				level->ExportVector(x, "it" + to_string(this->IterationCount) + "_sol_afterPostSmoothing");

			auto flopPostSmooth = result._iterationComputationalWork - flopBeforePostSmooth;
			//cout << "- flopPostSmooth = " << flopPostSmooth << endl;
		}
	}

private:
	Vector FCGForKCycle(Level* level, Vector& r, IterationResult& result)
	{
		const SparseMatrix& A = *level->OperatorMatrix;
		double t = 0.25;
		Vector x;

		// Apply multigrid preconditioner: c = Prec(r)
		Vector c = Vector::Zero(r.rows());
		bool cEquals0 = true;
		bool computeAc = true;
		IterationResult precResult1;
		MultigridCycle(level, r, c, cEquals0, false, computeAc, precResult1);     result.AddCost(precResult1.SolvingComputationalWork());

		// 1st iteration of FCG
		//Vector v = A * c;                                                         result.AddCost(Cost::MatVec(A));
		Vector& v = precResult1.Ax; // v = A*c
		double rho1 = c.dot(v);                                                   result.AddCost(Cost::Dot(c));
		double alpha1 = c.dot(r);                                                 result.AddCost(Cost::Dot(c));

		double r_norm_old = r.norm();                                             result.AddCost(Cost::Norm(r));
		r -= (alpha1 / rho1) * v;                                                 result.AddCost(Cost::VectorDAXPY(r));
		double r_norm_new = r.norm();                                             result.AddCost(Cost::Norm(r));
		if (r_norm_new <= t * r_norm_old)
		{
			x = (alpha1 / rho1) * c;                                              result.AddCost(Cost::ConstantVec(x));
		}
		else
		{
			// Apply multigrid preconditioner: d = Prec(r)
			Vector d = Vector::Zero(r.rows());
			bool dEquals0 = true;
			bool computeAd = true;
			IterationResult precResult2;
			MultigridCycle(level, r, d, dEquals0, false, computeAd, precResult2);  result.AddCost(precResult2.SolvingComputationalWork());

			// 2nd iteration of FCG
			//Vector w = A * d;                                                     result.AddCost(Cost::MatVec(A));
			Vector& w = precResult2.Ax; // w = A*d
			double gamma = d.dot(v);                                              result.AddCost(Cost::Dot(d));
			double beta = d.dot(w);                                               result.AddCost(Cost::Dot(d));
			double alpha2 = d.dot(r);                                             result.AddCost(Cost::Dot(d));
			double rho2 = beta - gamma * gamma / rho1;
			x = (alpha1/rho1 - gamma*alpha2/(rho1*rho2))*c + (alpha2 / rho2) * d; result.AddCost(2*Cost::ConstantVec(x) + Cost::AddVec(x));
		}
		return x;
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

		if (MultigridType == MGType::h_Multigrid)
		{
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
			else if (CoarseningStgy == CoarseningStrategy::MultiplePairwiseAggregation)
				os << "multiple pairwise aggregation [-cs mpa -coarsening-factor " << Utils::ProgramArgs.Solver.MG.CoarseningFactor << "]" << endl;
			else if (CoarseningStgy == CoarseningStrategy::AgglomerationCoarseningByFaceNeighbours)
				os << "agglomeration by face neighbours [-cs n]" << endl;
			else if (CoarseningStgy == CoarseningStrategy::MultipleAgglomerationCoarseningByFaceNeighbours)
				os << "multiple agglomeration by face neighbours [-cs mn -coarsening-factor " << Utils::ProgramArgs.Solver.MG.CoarseningFactor << "]" << endl;
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
		}

		//Level* level = this->CreateFineLevel();
		Smoother* preSmoother = SmootherFactory::Create(PreSmootherCode, PreSmoothingIterations, BlockSizeForBlockSmoothers, RelaxationParameter);
		Smoother* postSmoother = SmootherFactory::Create(PostSmootherCode, PostSmoothingIterations, BlockSizeForBlockSmoothers, RelaxationParameter);
		os << "\t" << "Pre-smoothing           : " << *preSmoother << endl;
		os << "\t" << "Post-smoothing          : " << *postSmoother;
		os.flush();
		//delete level;
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
			os << "\t" << "Coarse solver           : ";
			Multigrid* coarseMG = dynamic_cast<Multigrid*>(this->CoarseSolver);
			if (coarseMG)
				os << endl << "\t------------------------" << endl << "\t";
			os << (*this->CoarseSolver);
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