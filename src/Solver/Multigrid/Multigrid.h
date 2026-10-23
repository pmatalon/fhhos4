#pragma once
#include "../IterativeSolver.h"
#include "../../Utils/SparseMatrixOps.h"
#include "SmootherFactory.h"
#include "Level.h"
using namespace std;

class Multigrid : public IterativeSolver
{
protected:
	Level* _fineLevel = nullptr;
	bool _mixedPrecision = false; // the levels in single precision (MixedPrecision, if the subclass supports it)
private:
	bool _automaticNumberOfLevels;
	int _nLevels;
	SparseMatrix _coarseOperator;  // levels in single precision: the coarsest operator in double, for the coarse solver
public:
	ExportModule Out;
	Solver* CoarseSolver = nullptr;
	CycleType Cycle = CycleType::V;
	int WLoops = 1; // 1 --> V-cycle, >= 2 --> W-cycle 
	KrylovMethod KCycleKrylov = KrylovMethod::FCG; // the iterations of the K-cycle on the coarse levels
	// Mixed precision (U-AMG only; otherwise ignored): the matrices of the levels stored and the cycles run in single
	// precision. The operators, transfers and smoothers of the levels in float (the setup in double, its results rounded
	// when stored), the vectors of the cycle in float, its sums and dot products in double, the coarsest solve in double
	// (on a copy of the coarsest operator). Out of the cycle, in double: the matrix of the solver (A x for the outer Krylov method), the residuals
	// (each cycle corrects the residual of the current x: iterative refinement), x. With the options that the single
	// precision doesn't support (see UncondensedAMG::MixedPrecisionApplies()), the levels are in double.
	bool MixedPrecision = false;
	int MatrixMaxSizeForCoarsestLevel;
	bool UseGalerkinOperator = false;
	SmootherType PreSmoother = SmootherType::BlockGaussSeidel;
	SmootherType PostSmoother = SmootherType::BackwardBlockGaussSeidel;
	int PreSmoothingIterations = 1;
	int PostSmoothingIterations = 1;
	int CoarseLevelChangeSmoothingCoeff = 0;
	SmoothingChangeOperator CoarseLevelChangeSmoothingOperator = SmoothingChangeOperator::Add;
	int BlockSizeForBlockSmoothers = -1;
	double RelaxationParameter = 1;
	HP_CoarsStgy HP_CS = HP_CoarsStgy::H_only;
	H_CoarsStgy H_CS = H_CoarsStgy::StandardCoarsening;
	P_CoarsStgy P_CS = P_CoarsStgy::Minus2;
	FaceCoarseningStrategy FaceCoarseningStgy = FaceCoarseningStrategy::InterfaceCollapsing;
	FaceCollapsing BdryFaceCollapsing = FaceCollapsing::OnlyCollinear;
	double CoarseningFactor = 2;
	int CoarsePolyDegree = 1; // used in p-Multigrid
	int NumberOfMeshes = 0;
	bool ApplyZeroMeanCondition = false;
	bool EnforceCompatibilityCondition = false;
	bool ExportComponents = false;
	bool ExportIterationVectors = false; // the vectors of each step of the cycle, at each iteration (-export mgit)
	bool DoNotCreateLevels = false;

	Timer IntergridTransferTimer;
	Timer CoarseSolverTimer;
	MFlops IntergridTransferCost = 0;
	MFlops CoarseSolverCost = 0;

	Multigrid(int nLevels) : IterativeSolver()
	{
		//this->MultigridType = multigridType;
		//if (multigridType == MGType::h_Multigrid)
		//{
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
		/*}
		else
		{
			this->_automaticNumberOfLevels = true;
			this->_nLevels = 0;
			this->MatrixMaxSizeForCoarsestLevel = 0;
		}*/
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

		SetupLevels(A, nullptr);
	}

	// On several processes: the levels distributed by rows (see Level::DistributedOperator)
	void Setup(const DistributedMatrix& A) override
	{
		SetDistributedMatrix(A);
		SetupLevels(A.Diag, &A);
	}

private:
	// The levels of the fine operator A, distributed (distributedA, A being its block Diag) or not (null)
	void SetupLevels(const SparseMatrix& A, const DistributedMatrix* distributedA)
	{
		Utils::Log() << "Setup..." << endl;

		if (!this->_fineLevel)
			this->_fineLevel = this->CreateFineLevel();

		this->_fineLevel->OperatorMatrix = &A;
		this->_fineLevel->DistributedOperator = distributedA;
		this->_fineLevel->PreSmoother = this->_fineLevel->CreateSmoother(PreSmoother, PreSmoothingIterations, BlockSizeForBlockSmoothers, RelaxationParameter);
		this->_fineLevel->PostSmoother = this->_fineLevel->CreateSmoother(PostSmoother, PostSmoothingIterations, BlockSizeForBlockSmoothers, RelaxationParameter);
		this->_fineLevel->ExportComponents = ExportComponents;
		this->_fineLevel->SetExportModule(this->Out);

		Level* currentLevel = this->_fineLevel;
		bool noCoarserMeshProvided = false;
		bool coarsestPossibleMeshReached = false;
		int levelNumber = 0;
		int numberOfMeshes = 1;
		double operatorComplexity = 0;
		double fineNNZ = (double)this->_fineLevel->OperatorNonZeros();
		double gridComplexity = 0;
		double fineNUnknowns = (double)this->_fineLevel->NUnknowns();

		// Type of coarsening that must be performed
		CoarseningType coarseningType = FirstCoarseningType();


		// Creation of the coarse levels
		while (CoarseLevelNeeded(currentLevel, levelNumber))
		{
			coarseningType = ChangeCoarseningType(currentLevel, coarseningType);

			// Mesh coarsening if needed
			if (coarseningType == CoarseningType::H || coarseningType == CoarseningType::HP)
			{
				if (this->NumberOfMeshes == 0 || numberOfMeshes < this->NumberOfMeshes)
				{
					// Mesh coarsening
					currentLevel->CoarsenMesh(this->H_CS, this->FaceCoarseningStgy, this->BdryFaceCollapsing, this->CoarseningFactor, noCoarserMeshProvided, coarsestPossibleMeshReached);

					if (noCoarserMeshProvided || coarsestPossibleMeshReached)
					{
						// If no coarser mesh, then p-coarsening if needed
						if ((this->HP_CS == HP_CoarsStgy::H_then_P || this->HP_CS == HP_CoarsStgy::HP_then_P) && currentLevel->PolynomialDegree() > this->CoarsePolyDegree)
						{
							if (currentLevel->PolynomialDegree() > this->CoarsePolyDegree)
								coarseningType = CoarseningType::P;
							else
							{
								if (coarsestPossibleMeshReached)
									Utils::Warning("Cannot define a lower level because the mesh cannot be coarsened any more and the degree cannot be lowered.");
								if (noCoarserMeshProvided)
									Utils::Warning("Cannot define a lower level because no coarser mesh has been provided and the degree cannot be lowered.");
								break; // stop building levels
							}
						}
						else
						{
							if (coarsestPossibleMeshReached)
								Utils::Warning("Cannot define a lower level because the mesh cannot be coarsened any more.");
							if (noCoarserMeshProvided)
								Utils::Warning("Cannot define a lower level because no coarser mesh has been provided.");
							break; // stop building levels
						}
					}
					else
						numberOfMeshes++;
				}
				else if (this->HP_CS == HP_CoarsStgy::H_then_P || this->HP_CS == HP_CoarsStgy::HP_then_P)
				{
					if (currentLevel->PolynomialDegree() > this->CoarsePolyDegree)
						coarseningType = CoarseningType::P;
					else
					{
						Utils::Warning("Cannot define a lower level because the requested number of meshes is reached and the degree cannot be lowered.");
						break; // stop building levels
					}
				}
				else
				{
					Utils::Warning("Cannot define a lower level because the requested number of meshes is reached.");
					break; // stop building levels
				}
			}

			if (coarseningType == CoarseningType::P && currentLevel->PolynomialDegree() == this->CoarsePolyDegree)
			{
				Utils::Warning("Cannot define a lower level because the degree cannot be lowered.");
				break; // stop building levels
			}

			int coarseDegree = -1;
			int currentDegree = currentLevel->PolynomialDegree();
			if (coarseningType == CoarseningType::P || coarseningType == CoarseningType::HP)
			{
				if (coarseningType == CoarseningType::P)
				{
					if (this->P_CS == P_CoarsStgy::Minus1)
						coarseDegree = currentDegree - 1;
					else if (this->P_CS == P_CoarsStgy::Minus2)
						coarseDegree = currentDegree - 2;
					else if (this->P_CS == P_CoarsStgy::DivideBy2)
						coarseDegree = currentDegree / 2;
					else if (this->P_CS == P_CoarsStgy::DirectToLow)
						coarseDegree = this->CoarsePolyDegree;
					else
						assert(false);
				}
				else
					coarseDegree = currentDegree - 1;

				if (this->HP_CS == HP_CoarsStgy::P_then_HP)
				{
					int degreeWhenHPCoarseningMustStart = (this->NumberOfMeshes - 1) + this->CoarsePolyDegree;
					coarseDegree = max(coarseDegree, degreeWhenHPCoarseningMustStart);
				}
				else
					coarseDegree = max(coarseDegree, this->CoarsePolyDegree);
			}
			else
				coarseDegree = currentDegree;

			// Build coarse level
			levelNumber++;
			if (this->HP_CS != HP_CoarsStgy::H_only)
				Utils::Log() << "\tCreation of level " << levelNumber << " (" << (coarseningType == CoarseningType::H ? "h" : (coarseningType == CoarseningType::P ? "p" : "hp")) << "-coarsening)" << endl;
			Level* coarseLevel = CreateCoarseLevel(currentLevel, coarseningType, coarseDegree);
			coarseLevel->ComesFrom = coarseningType;
			coarseLevel->SetExportModule(this->Out);

			// Smoothers on coarse level
			int preSmoothingIterations = currentLevel->PreSmoother->Iterations();
			int postSmoothingIterations = currentLevel->PostSmoother->Iterations();
			if (CoarseLevelChangeSmoothingOperator == SmoothingChangeOperator::Add)
			{
				preSmoothingIterations += CoarseLevelChangeSmoothingCoeff;
				postSmoothingIterations += CoarseLevelChangeSmoothingCoeff;
			}
			else if (CoarseLevelChangeSmoothingOperator == SmoothingChangeOperator::Subtract)
			{
				preSmoothingIterations -= CoarseLevelChangeSmoothingCoeff;
				postSmoothingIterations -= CoarseLevelChangeSmoothingCoeff;
			}
			else if (CoarseLevelChangeSmoothingOperator == SmoothingChangeOperator::Multiply)
			{
				preSmoothingIterations *= CoarseLevelChangeSmoothingCoeff;
				postSmoothingIterations *= CoarseLevelChangeSmoothingCoeff;
			}
			else
				assert(false);


			int blockSize = this->HP_CS == HP_CoarsStgy::H_only ? BlockSizeForBlockSmoothers : coarseLevel->BlockSizeForBlockSmoothers();
			coarseLevel->PreSmoother  = coarseLevel->CreateSmoother( PreSmoother,  preSmoothingIterations, blockSize, RelaxationParameter);
			coarseLevel->PostSmoother = coarseLevel->CreateSmoother(PostSmoother, postSmoothingIterations, blockSize, RelaxationParameter);

			coarseLevel->UseGalerkinOperator = UseGalerkinOperator;
			coarseLevel->ExportComponents = ExportComponents;

			/*// FCG if K-cycle
			if (this->Cycle == CycleType::K)
			{
				FlexibleConjugateGradient* fcg = new FlexibleConjugateGradient();
				coarseLevel->FCG = fcg;

				Multigrid* mg = new Multigrid();
				mg->_fineLevel = coarseLevel;
				mg->DoNotCreateLevels = true;
				mg->Cycle = this->Cycle;
				mg->WLoops = this->WLoops;
				fcg->Precond = SolverPreconditioner(mg);

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

			operatorComplexity += currentLevel->OperatorNonZeros() / fineNNZ;
			gridComplexity += currentLevel->NUnknowns() / fineNUnknowns;

			currentLevel = coarseLevel;
		}

		_nLevels = levelNumber + 1;
		currentLevel->Setup();

		operatorComplexity += currentLevel->OperatorNonZeros() / fineNNZ;
		gridComplexity += currentLevel->NUnknowns() / fineNUnknowns;

		this->InitializeCoarseSolver();
		this->SetupCoarseSolver();

		/*if (this->Cycle == CycleType::K)
		{
			Level* level = this->_fineLevel->CoarserLevel;
			while (level)
			{
				Multigrid* mg = static_cast<Multigrid*>(level->FCG->Precond.GetSolver());
				mg->CoarseSolver = this->CoarseSolver;
				level = level->CoarserLevel;
			}
		}*/

		Utils::Log() << endl;
		Utils::Log() << "\tLevels             : " << _nLevels << endl;
		Utils::Log() << "\tOperator complexity: " << operatorComplexity << endl;
		Utils::Log() << "\tGrid complexity    : " << gridComplexity << endl;
		Utils::Log() << endl;

		if (this->WLoops > 1)
			PrintCycleSchema();
	}

protected:
	void DoBeforeSolving() override
	{
		this->IntergridTransferCost = 0;
		this->CoarseSolverCost = 0;
	}

	virtual Level* CreateFineLevel() const
	{
		Utils::FatalError("The method CreateFineLevel() must be implemented in the Level subclass.");
		return nullptr;
	}
	virtual Level* CreateCoarseLevel(Level* fineLevel, CoarseningType coarseningType, int coarseDegree)
	{
		Utils::FatalError("The method CreateCoarseLevel() must be implemented in the Level subclass.");
		return nullptr;
	}

	virtual void InitializeCoarseSolver()
	{
		if (!this->CoarseSolver)
			this->CoarseSolver = new EigenSparseLU();
	}

	void SetupCoarseSolver()
	{
		Utils::Log() << "\t\tSetup coarse solver..." << endl;

		Level* level = this->_fineLevel;
		while (level->CoarserLevel != nullptr)
			level = level->CoarserLevel;

		if (level->IsDistributed())
			this->CoarseSolver->Setup(*level->DistributedOperator);
		else
		{
			if (!level->OperatorMatrix) // a level in single precision (not the finest): its operator in double
			{
				_coarseOperator = SparseMatrixOps::Cast<double>(*level->SingleOperatorMatrix);
				level->OperatorMatrix = &_coarseOperator;
			}
			this->CoarseSolver->Setup(*level->OperatorMatrix);
		}
	}

private:
	CoarseningType FirstCoarseningType()
	{
		if (this->HP_CS == HP_CoarsStgy::H_only || this->HP_CS == HP_CoarsStgy::H_then_P)
			return CoarseningType::H;
		else if (this->HP_CS == HP_CoarsStgy::P_only)
		{
			if (_fineLevel->PolynomialDegree() > this->CoarsePolyDegree)
				return CoarseningType::P;
			else
				Utils::FatalError("p-coarsening cannot be performed because the polynomial degree of the fine level = " + to_string(_fineLevel->PolynomialDegree()));
		}
		else if (this->HP_CS == HP_CoarsStgy::P_then_H)
			return _fineLevel->PolynomialDegree() > this->CoarsePolyDegree ? CoarseningType::P : CoarseningType::H;
		else if (this->HP_CS == HP_CoarsStgy::P_then_HP)
		{
			if (_fineLevel->PolynomialDegree() > this->CoarsePolyDegree)
			{
				if (this->NumberOfMeshes == 0)
					Utils::FatalError("With -hp-cs p_hp, the number of meshes must be fixed using -num-meshes.");
				return CoarseningType::P;
			}
			else
				Utils::FatalError("Neither p- or hp-coarsening can be performed because the polynomial degree of the fine level = " + to_string(_fineLevel->PolynomialDegree()));
		}
		else if (this->HP_CS == HP_CoarsStgy::HP_then_H)
			return _fineLevel->PolynomialDegree() > this->CoarsePolyDegree ? CoarseningType::HP : CoarseningType::H;
		else if (this->HP_CS == HP_CoarsStgy::HP_then_P)
		{
			if (_fineLevel->PolynomialDegree() > this->CoarsePolyDegree)
				return CoarseningType::HP;
			else
				Utils::FatalError("Neither hp- or p-coarsening can be performed because the polynomial degree of the fine level = " + to_string(_fineLevel->PolynomialDegree()));
		}
		else if (this->HP_CS == HP_CoarsStgy::Alternate)
			return _fineLevel->PolynomialDegree() > this->CoarsePolyDegree ? CoarseningType::H : CoarseningType::P; // return H if we want to start with P, because ChangeCoarseningType() is applied after
		else
			Utils::FatalError("This hp-cs is not implemented.");

		return CoarseningType::H; // just to avoid warning
	}

	CoarseningType ChangeCoarseningType(Level* currentLevel, CoarseningType currentCoarseningType)
	{
		if (this->HP_CS == HP_CoarsStgy::P_then_HP && currentCoarseningType == CoarseningType::P)
		{
			if (this->NumberOfMeshes == 0)
				Utils::FatalError("With -hp-cs p_hp, the number of meshes must be fixed using -num-meshes.");

			int degreeWhenHPCoarseningMustStart = (this->NumberOfMeshes - 1) + this->CoarsePolyDegree;
			return currentLevel->PolynomialDegree() > degreeWhenHPCoarseningMustStart ? CoarseningType::P : CoarseningType::HP;
		}
		else if (this->HP_CS == HP_CoarsStgy::P_then_H && currentCoarseningType == CoarseningType::P && currentLevel->PolynomialDegree() == this->CoarsePolyDegree)
			return CoarseningType::H;
		else if (currentCoarseningType == CoarseningType::HP && currentLevel->PolynomialDegree() == this->CoarsePolyDegree)
		{
			if (this->HP_CS == HP_CoarsStgy::HP_then_H)
				return CoarseningType::H;
		}
		else if (this->HP_CS == HP_CoarsStgy::Alternate)
		{
			if (currentCoarseningType == CoarseningType::H && currentLevel->PolynomialDegree() > this->CoarsePolyDegree)
				return CoarseningType::P;
			else
				return CoarseningType::H;
		}
		return currentCoarseningType;
	}

	bool CoarseLevelNeeded(Level* currentLevel, int levelNumber)
	{
		bool stillNeedToCoarsen = (_automaticNumberOfLevels && currentLevel->NUnknowns() > MatrixMaxSizeForCoarsestLevel) || (levelNumber < _nLevels - 1);
		if (!stillNeedToCoarsen)
			return false;
		if (this->HP_CS == HP_CoarsStgy::P_only && currentLevel->PolynomialDegree() == CoarsePolyDegree)
			return false; // cannot p-coarsen anymore
		return true;
	}

	IterationResult ExecuteOneIteration(const Vector& b, Vector& x, bool& xEquals0, bool computeResidual, bool computeAx, const IterationResult& oldResult) override
	{
		IterationResult result(oldResult);

		if (_mixedPrecision)
			MixedPrecisionCycle(b, x, xEquals0, computeResidual, computeAx, result);
		else
			MultigridCycle(this->_fineLevel, b, x, xEquals0, computeResidual, computeAx, result);
		result.SetX(x);
		return result;
	}

	void MultigridCycle(Level* level, const Vector& b, Vector& x, bool& xEquals0, bool computeResidual, bool computeAx, IterationResult& result)
	{
		CycleIn<double>(level, b, x, xEquals0, computeResidual ? &result.Residual : nullptr, computeAx ? &result.Ax : nullptr, result);
	}

private:
	// Mixed precision: one cycle in single precision on the residual r = b - A x (b if x = 0), in double with the matrix of
	// the solver, its correction e added to x (x = e if x = 0); the residual and A x in double (A full: a parallel product)
	void MixedPrecisionCycle(const Vector& b, Vector& x, bool& xEquals0, bool computeResidual, bool computeAx, IterationResult& result)
	{
		const SparseMatrix& A = *this->_fineLevel->OperatorMatrix;
		BigNumber n = b.rows();
		Vector r;
		if (!xEquals0)
		{
			r = b - A * x;                                                         result.AddWorkInFlops(Cost::DAXPY(A));
		}
		const Vector& rhs = xEquals0 ? b : r;
		FloatVector rSingle(n), eSingle;
		#pragma omp parallel for if (n > 20000)
		for (BigNumber i = 0; i < n; i++)
			rSingle[i] = (float)rhs[i];
		bool eSingleEquals0 = true;
		CycleIn<float>(this->_fineLevel, rSingle, eSingle, eSingleEquals0, nullptr, nullptr, result);
		if (xEquals0)
			x.resize(n);
		#pragma omp parallel for if (n > 20000)
		for (BigNumber i = 0; i < n; i++)
			x[i] = xEquals0 ? (double)eSingle[i] : x[i] + eSingle[i];
		xEquals0 = false;
		if (computeAx || computeResidual)
		{
			result.Ax = A * x;                                                     result.AddWorkInFlops(Cost::MatVec(A));
			if (computeResidual)
			{
				result.Residual = b - result.Ax;                                   result.AddWorkInFlops(Cost::AddVec(b));
			}
		}
	}

	// The cycle on the level in precision T: double, or float (the levels set up in single precision, the coarsest solve
	// in double); the residual and A x after the post-smoothing in residual and Ax if not null
	template <typename T>
	void CycleIn(Level* level, const VectorOf<T>& b, VectorOf<T>& x, bool& xEquals0, VectorOf<T>* residual, VectorOf<T>* Ax, IterationResult& result)
	{
		constexpr bool single = std::is_same_v<T, float>;

		if constexpr (!single)
		{
			if (this->ExportIterationVectors)
				level->ExportVector(b, "it" + to_string(this->IterationCount) + "_b");
		}

		if (level->IsCoarsestLevel())
		{
			                                                                         CoarseSolverTimer.Start();
			if constexpr (single)
				x = CoarseSolver->Solve(Vector(b.template cast<double>())).template cast<float>();
			else
				x = CoarseSolver->Solve(b);
			                                                                         result.AddWorkInMFlops(CoarseSolver->SolvingComputationalWork);
			                                                                         CoarseSolverTimer.Pause();
																					 CoarseSolverCost += CoarseSolver->SolvingComputationalWork;
			xEquals0 = false;
		}
		else
		{
			if constexpr (!single)
			{
				if (this->ExportIterationVectors)
					level->ExportVector(x, "it" + to_string(this->IterationCount) + "_sol_beforePreSmoothing");
			}

			//---------------//
			// Pre-smoothing //
			//---------------//

			//auto flopBeforePreSmooth = result.IterationComputationalWork();

			VectorOf<T> r;
			if (level->PreSmoother->CanOptimizeResidualComputation())
			{
				// Smooth and compute residual
				level->PreSmoother->SmoothAndComputeResidualOrAx(x, b, xEquals0, true, false);    result.AddWorkInMFlops(level->PreSmoother->SolvingComputationalWork());
				
				// Move residual vector
				r = level->PreSmoother->template TakeResidual<T>();
			}
			else
			{
				// Smooth
				level->PreSmoother->Smooth(x, b, xEquals0);                               result.AddWorkInMFlops(level->PreSmoother->SolvingComputationalWork());
				// Residual computation
				r = level->ResidualOf(b, x);                                              result.AddWorkInFlops(Cost::DAXPY(level->template LocalOperatorNonZeros<T>()));
			}
			
			if constexpr (!single)
			{
				if (this->ExportIterationVectors)
					level->ExportVector(x, "it" + to_string(this->IterationCount) + "_sol_afterPreSmoothing");
			}

			//auto flopPreSmoothAndRes = result.IterationComputationalWork() - flopBeforePreSmooth;
			//cout << "- flopPreSmoothAndRes = " << flopPreSmoothAndRes << endl;

			//------------------------------------------------//
			// Restriction of the residual on the coarse grid //
			//------------------------------------------------//

			auto flopBeforeRestrict = result.IterationComputationalWork();
			IntergridTransferTimer.Start();

			VectorOf<T> rc = level->Restrict(r);                                          result.AddWorkInFlops(level->RestrictCost());
			                                                                              this->IntergridTransferCost += level->RestrictCost()*1e-6;
			IntergridTransferTimer.Pause();
			auto flopRestrict = result.IterationComputationalWork() - flopBeforeRestrict;
			assert(abs(round(flopRestrict*1e6 - level->RestrictCost())) == 0);
			//cout << "- flopRestrict = " << flopRestrict << endl;

			//--------------------------------------------------//
			// Residual equation Ae=r solved on the coarse grid //
			//--------------------------------------------------//

			if constexpr (!single)
			{
				if (this->EnforceCompatibilityCondition)
				{
					level->CoarserLevel->EnforceCompatibilityCondition(rc);               result.AddWorkInFlops(level->EnforceCompatibilityConditionCost(rc));
				}
			}

			VectorOf<T> ec;
			if (this->Cycle == CycleType::V || this->Cycle == CycleType::W)
			{
				ec = VectorOf<T>::Zero(rc.rows());
				bool ecEquals0 = true;
				for (int i = 0; i < this->WLoops; ++i)
				{
					IterationResult coarseResult;
					CycleIn<T>(level->CoarserLevel, rc, ec, ecEquals0, nullptr, nullptr, coarseResult);   result.AddWorkInMFlops(coarseResult.SolvingComputationalWork());
					if (level->CoarserLevel->IsCoarsestLevel())
						break;
				}
			}
			else if (this->Cycle == CycleType::K)
			{
				IterationResult coarseResult;
				if (level->CoarserLevel->IsCoarsestLevel())
				{
					bool ecEquals0 = true;
					CycleIn<T>(level->CoarserLevel, rc, ec, ecEquals0, nullptr, nullptr, coarseResult); // exact solution
					                                                                               result.AddWorkInMFlops(coarseResult.SolvingComputationalWork());
				}
				else
				{
					//level->CoarserLevel->FCG->Solve(rc, ecEqualZero, ec);                    result.AddWorkInMFlops(level->CoarserLevel->FCG->SolvingComputationalWork);
					ec = KrylovForKCycle<T>(level->CoarserLevel, rc, coarseResult);                result.AddWorkInMFlops(coarseResult.SolvingComputationalWork());
				}
			}

			//------------------------//
			// Coarse-grid correction //
			//------------------------//

			if constexpr (!single)
			{
				if (this->ExportIterationVectors)
				{
					level->ExportVector(ec, "it" + to_string(this->IterationCount) + "_ce");
					level->ExportVector(x, "it" + to_string(this->IterationCount) + "_sol");
					Vector cgc = level->Prolong(ec);
					level->ExportVector(cgc, "it" + to_string(this->IterationCount) + "_cgc");
				}
			}

			// Prolongation
			IntergridTransferTimer.Start();

			x += level->Prolong(ec);                                                  result.AddWorkInFlops(level->ProlongCost() + Cost::AddVec(x));
			                                                                          this->IntergridTransferCost += (level->ProlongCost() + Cost::AddVec(x))*1e-6;
			IntergridTransferTimer.Pause();

			if constexpr (!single)
			{
				if (this->ExportIterationVectors)
					level->ExportVector(x, "it" + to_string(this->IterationCount) + "_sol_cgc");
			}

			//----------------//
			// Post-smoothing //
			//----------------//

			assert(!xEquals0);

			//auto flopBeforePostSmooth = result.IterationComputationalWork();

			if (level->PostSmoother->CanOptimizeResidualComputation() && (residual || Ax))
			{
				// Smooth and compute residual
				level->PostSmoother->SmoothAndComputeResidualOrAx(x, b, xEquals0, residual != nullptr, Ax != nullptr);    result.AddWorkInMFlops(level->PostSmoother->SolvingComputationalWork());
				if (residual)
					*residual = level->PostSmoother->template TakeResidual<T>();
				if (Ax)
					*Ax = level->PostSmoother->template TakeAx<T>();
			}
			else
			{
				// Smooth
				level->PostSmoother->Smooth(x, b, xEquals0);                          result.AddWorkInMFlops(level->PostSmoother->SolvingComputationalWork());
				if (Ax)
				{
					// Residual computation
					*Ax = level->Apply(x);                                            result.AddWorkInFlops(Cost::SpMatVec(level->template LocalOperatorNonZeros<T>()));
					if (residual)
					{
						*residual = b - *Ax;                                          result.AddWorkInFlops(Cost::AddVec(b));
					}
				}
				else if (residual)
				{
					*residual = level->ResidualOf(b, x);                              result.AddWorkInFlops(Cost::DAXPY(level->template LocalOperatorNonZeros<T>()));
				}
			}

			if constexpr (!single)
			{
				if (this->ExportIterationVectors)
					level->ExportVector(x, "it" + to_string(this->IterationCount) + "_sol_afterPostSmoothing");
			}

			//auto flopPostSmooth = result.IterationComputationalWork() - flopBeforePostSmooth;
			//cout << "- flopPostSmooth = " << flopPostSmooth << endl;

			if constexpr (!single)
			{
				if (this->ApplyZeroMeanCondition)
				{
					level->ApplyZeroMeanCondition(x);                                 result.AddWorkInFlops(level->ApplyZeroMeanConditionCost(x));
				}
			}
		}
	}

	// 1 or 2 iterations of FCG or GCR (KCycleKrylov) preconditioned by the cycle on this level: Notay, "An
	// aggregation-based algebraic multigrid method", ETNA 37 (2010), Algorithm 3.2, step 4. The two methods differ by
	// their inner products only: FCG minimizes the A-norm of the error, GCR the residual norm. In precision T (the dot
	// products in double).
	template <typename T>
	VectorOf<T> KrylovForKCycle(Level* level, VectorOf<T>& r, IterationResult& result)
	{
		//const SparseMatrix& A = *level->OperatorMatrix;
		double t = 0.25;
		bool gcr = this->KCycleKrylov == KrylovMethod::GCR;
		VectorOf<T> x;

		// Apply multigrid preconditioner: c = Prec(r)
		VectorOf<T> c = VectorOf<T>::Zero(r.rows());
		bool cEquals0 = true;
		IterationResult precResult1;
		VectorOf<T> v; // v = A*c
		CycleIn<T>(level, r, c, cEquals0, nullptr, &v, precResult1);             result.AddWorkInMFlops(precResult1.SolvingComputationalWork());

		// 1st iteration
		//Vector v = A * c;                                                         result.AddWorkInFlops(Cost::MatVec(A));
		double rho1 = gcr ? Parallel::SquaredNorm(v) : Parallel::Dot(c, v);       result.AddWorkInFlops(Cost::Dot(c));
		double alpha1 = gcr ? Parallel::Dot(v, r) : Parallel::Dot(c, r);          result.AddWorkInFlops(Cost::Dot(c));

		double r_norm_old = Parallel::Norm(r);                                    result.AddWorkInFlops(Cost::Norm(r));
		r -= T(alpha1 / rho1) * v;                                                result.AddWorkInFlops(Cost::VectorDAXPY(r));
		double r_norm_new = Parallel::Norm(r);                                    result.AddWorkInFlops(Cost::Norm(r));
		if (r_norm_new <= t * r_norm_old)
		{
			x = T(alpha1 / rho1) * c;                                             result.AddWorkInFlops(Cost::ConstantVec(x));
		}
		else
		{
			// Apply multigrid preconditioner: d = Prec(r)
			VectorOf<T> d = VectorOf<T>::Zero(r.rows());
			bool dEquals0 = true;
			IterationResult precResult2;
			VectorOf<T> w; // w = A*d
			CycleIn<T>(level, r, d, dEquals0, nullptr, &w, precResult2);         result.AddWorkInMFlops(precResult2.SolvingComputationalWork());

			// 2nd iteration
			//Vector w = A * d;                                                     result.AddWorkInFlops(Cost::MatVec(A));
			double gamma = gcr ? Parallel::Dot(w, v) : Parallel::Dot(d, v);       result.AddWorkInFlops(Cost::Dot(d));
			double beta = gcr ? Parallel::SquaredNorm(w) : Parallel::Dot(d, w);   result.AddWorkInFlops(Cost::Dot(d));
			double alpha2 = gcr ? Parallel::Dot(w, r) : Parallel::Dot(d, r);      result.AddWorkInFlops(Cost::Dot(d));
			double rho2 = beta - gamma * gamma / rho1;
			x = T(alpha1/rho1 - gamma*alpha2/(rho1*rho2))*c + T(alpha2 / rho2) * d; result.AddWorkInFlops(2*Cost::ConstantVec(x) + Cost::AddVec(x));
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
		BeginSerialize(os);

		os << "\t" << "Cycle                   : ";
		if (this->Cycle == CycleType::K)
			os << "K";
		else if (this->WLoops == 1)
			os << "V";
		else if (this->WLoops == 2)
			os << "W";
		else
			os << "W(" << this->WLoops << " loops)";
		os << "(" << PreSmoothingIterations << "," << PostSmoothingIterations;
		switch (CoarseLevelChangeSmoothingOperator)
		{
			case SmoothingChangeOperator::Add:
				if (CoarseLevelChangeSmoothingCoeff > 0)
					os << ",+" << CoarseLevelChangeSmoothingCoeff;
				break;
			case SmoothingChangeOperator::Subtract:
				if (CoarseLevelChangeSmoothingCoeff > 0)
					os << ",-" << CoarseLevelChangeSmoothingCoeff;
				break;
			case SmoothingChangeOperator::Multiply:
				if (CoarseLevelChangeSmoothingCoeff > 1)
					os << ",*" << CoarseLevelChangeSmoothingCoeff;
				break;
		}
		os << ")";
		if (this->Cycle == CycleType::K)
			os << ", " << (this->KCycleKrylov == KrylovMethod::GCR ? "GCR" : "FCG") << " on the coarse levels";
		os << endl;

		if (this->HP_CS != HP_CoarsStgy::P_only)
		{
			os << "\t" << "Levels                  : ";
			if (_automaticNumberOfLevels && _nLevels == 0)
				os << "automatic coarsening until matrix size <= " << MatrixMaxSizeForCoarsestLevel << endl;
			else if (_automaticNumberOfLevels && _nLevels > 0)
				os << _nLevels << " (automatic)" << endl;
			else
				os << _nLevels << endl;

			os << "\t" << "Coarsening strategy     : ";
			if (H_CS == H_CoarsStgy::StandardCoarsening)
				os << "standard [-cs s]" << endl;
			else if (H_CS == H_CoarsStgy::GMSHSplittingRefinement)
				os << "GMSH refinement by splitting from coarse mesh [-cs r]" << endl;
			else if (H_CS == H_CoarsStgy::BeyRefinement)
				os << "Bey's refinement from coarse mesh [-cs b]" << endl;
			else if (H_CS == H_CoarsStgy::IndependentRemeshing)
				os << "independant remeshing [-cs m]" << endl;
			else if (H_CS == H_CoarsStgy::DoublePairwiseAggregation)
				os << "double pairwise aggregation [-cs dpa]" << endl;
			else if (H_CS == H_CoarsStgy::MultiplePairwiseAggregation)
				os << "multiple pairwise aggregation [-cs mpa -coarsening-factor " << CoarseningFactor << "]" << endl;
			else if (H_CS == H_CoarsStgy::AgglomerationCoarseningByFaceNeighbours)
				os << "agglomeration by face neighbours [-cs n]" << endl;
			else if (H_CS == H_CoarsStgy::MultipleAgglomerationCoarseningByFaceNeighbours)
				os << "multiple agglomeration by face neighbours [-cs mn -coarsening-factor " << CoarseningFactor << "]" << endl;
			else if (H_CS == H_CoarsStgy::AgglomerationCoarseningByMostCoplanarFaces)
				os << "agglomeration by most coplanar faces [-cs mcf]" << endl;
			else if (H_CS == H_CoarsStgy::AgglomerationCoarseningByClosestCenter)
				os << "agglomeration by closest element center [-cs cc]" << endl;
			else if (H_CS == H_CoarsStgy::AgglomerationCoarseningByClosestFace)
				os << "agglomeration by closest face center [-cs clf]" << endl;
			else if (H_CS == H_CoarsStgy::AgglomerationCoarseningByLargestInterface)
				os << "agglomeration by largest interface [-cs li]" << endl;
			else if (H_CS == H_CoarsStgy::FaceCoarsening)
				os << "face coarsening [-cs f]" << endl;
			else
				os << "unknown" << endl;

			if (!Utils::IsRefinementStrategy(H_CS) && H_CS != H_CoarsStgy::IndependentRemeshing)
			{
				os << "\t" << "Face coarsening strategy: ";
				if (FaceCoarseningStgy == FaceCoarseningStrategy::None)
					os << "none [-fcs n]" << endl;
				else if (FaceCoarseningStgy == FaceCoarseningStrategy::InterfaceCollapsing)
					os << "interface collapsing [-fcs c]" << endl;
				else if (FaceCoarseningStgy == FaceCoarseningStrategy::InterfaceCollapsingAndTryAggregInteriorToInterfaces)
					os << "interface collapsing and try aggregate interior faces too [-fcs i]" << endl;
				else
					os << "unknown" << endl;

				os << "\t" << "Boundary face collapsing: ";
				if (BdryFaceCollapsing == FaceCollapsing::Disabled)
					os << "disabled [-bfc d]" << endl;
				else if (BdryFaceCollapsing == FaceCollapsing::OnlyCollinear)
					os << "collinear only [-bfc c]" << endl;
				else if (BdryFaceCollapsing == FaceCollapsing::ByPairs)
					os << "by pairs [-bfc p]" << endl;
				else if (BdryFaceCollapsing == FaceCollapsing::Max)
					os << "maximum [-bfc m]" << endl;
				else
					os << "unknown" << endl;
			}
		}

		//Level* level = this->CreateFineLevel();
		Smoother* preSmoother = SmootherFactory::Create(PreSmoother, PreSmoothingIterations, BlockSizeForBlockSmoothers, RelaxationParameter);
		Smoother* postSmoother = SmootherFactory::Create(PostSmoother, PostSmoothingIterations, BlockSizeForBlockSmoothers, RelaxationParameter);
		os << "\t" << "Pre-smoothing           : " << *preSmoother << endl;
		os << "\t" << "Post-smoothing          : " << *postSmoother;
		os.flush();
		//delete level;
		delete preSmoother;
		delete postSmoother;
		if (PreSmoother == SmootherType::BlockJacobi || PostSmoother == SmootherType::BlockJacobi)
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

		EndSerialize(os);
	}

	void PrintCycleSchema()
	{
		// Level number of each successive step of the cycle
		vector<int> schema;
		Level* lastStoredLevel = nullptr;
		StoreCycleSchema(this->_fineLevel, schema, lastStoredLevel);

		for (int levelNumber = 0; levelNumber < this->NumberOfLevels(); levelNumber++)
		{
			for (int stepLevel : schema)
				Utils::Log() << (stepLevel == levelNumber ? " o" : "  ");
			Utils::Log() << endl;
		}
		Utils::Log() << endl;
	}

private:
	void StoreCycleSchema(Level* level, vector<int>& schema, Level*& lastStoredLevel)
	{
		if (lastStoredLevel != level)
		{
			schema.push_back(level->Number);
			lastStoredLevel = level;
		}

		if (!level->IsCoarsestLevel())
		{
			for (int i = 0; i < this->WLoops; ++i)
			{
				StoreCycleSchema(level->CoarserLevel, schema, lastStoredLevel);
				if (level->CoarserLevel->IsCoarsestLevel())
					break;
			}

			if (lastStoredLevel != level)
			{
				schema.push_back(level->Number);
				lastStoredLevel = level;
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