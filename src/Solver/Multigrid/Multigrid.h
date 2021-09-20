#pragma once
#include "../IterativeSolver.h"
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
	HP_CoarsStgy HP_CS = HP_CoarsStgy::H_only;
	H_CoarsStgy H_CS = H_CoarsStgy::StandardCoarsening;
	P_CoarsStgy P_CS = P_CoarsStgy::Minus2;
	FaceCoarseningStrategy FaceCoarseningStgy = FaceCoarseningStrategy::InterfaceCollapsing;
	double CoarseningFactor = 2;
	int CoarsePolyDegree = 1; // used in p-Multigrid
	int NumberOfMeshes = 0;
	bool ExportComponents = false;
	bool DoNotCreateLevels = false;

	Timer IntergridTransferTimer;
	Timer CoarseSolverTimer;
	size_t IntergridTransferCost = 0;
	size_t CoarseSolverCost = 0;

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
		int numberOfMeshes = 1;
		double operatorComplexity = 0;
		double fineNNZ = (double)this->_fineLevel->OperatorMatrix->nonZeros();
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
					currentLevel->CoarsenMesh(this->H_CS, this->FaceCoarseningStgy, this->CoarseningFactor, noCoarserMeshProvided, coarsestPossibleMeshReached);

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
				cout << "\tCreation of level " << levelNumber << " (" << (coarseningType == CoarseningType::H ? "h" : (coarseningType == CoarseningType::P ? "p" : "hp")) << "-coarsening)" << endl;
			Level* coarseLevel = CreateCoarseLevel(currentLevel, coarseningType, coarseDegree);
			coarseLevel->ComesFrom = coarseningType;

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


			int blockSize = this->HP_CS == HP_CoarsStgy::H_only ? BlockSizeForBlockSmoothers : coarseLevel->BlockSizeForBlockSmoothers();
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
		cout << "\t\tSetup coarse solver..." << endl;

		Level* level = this->_fineLevel;
		while (level->CoarserLevel != nullptr)
			level = level->CoarserLevel;

		this->CoarseSolver->Setup(*level->OperatorMatrix);
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
				Utils::FatalError("p-coarsening cannot be performed because the polynomial degree of the fine level = " + to_string(this->CoarsePolyDegree));
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
				Utils::FatalError("Neither p- or hp-coarsening can be performed because the polynomial degree of the fine level = " + to_string(this->CoarsePolyDegree));
		}
		else if (this->HP_CS == HP_CoarsStgy::HP_then_H)
			return _fineLevel->PolynomialDegree() > this->CoarsePolyDegree ? CoarseningType::HP : CoarseningType::H;
		else if (this->HP_CS == HP_CoarsStgy::HP_then_P)
		{
			if (_fineLevel->PolynomialDegree() > this->CoarsePolyDegree)
				return CoarseningType::HP;
			else
				Utils::FatalError("Neither hp- or p-coarsening can be performed because the polynomial degree of the fine level = " + to_string(this->CoarsePolyDegree));
		}
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
																					 CoarseSolverCost += CoarseSolver->SolvingComputationalWork;
			xEquals0 = false;
		}
		else
		{
			if (Utils::ProgramArgs.Actions.ExportMultigridIterationVectors)
				level->ExportVector(x, "it" + to_string(this->IterationCount) + "_sol_beforePreSmoothing");

			//---------------//
			// Pre-smoothing //
			//---------------//

			//auto flopBeforePreSmooth = result.IterationComputationalWork();

			Vector r;
			if (level->PreSmoother->CanOptimizeResidualComputation())
			{
				// Smooth and compute residual
				level->PreSmoother->SmoothAndComputeResidualOrAx(x, b, xEquals0, true, false);    result.AddCost(level->PreSmoother->SolvingComputationalWork());
				
				// Move residual vector
				r = level->PreSmoother->Residual();
			}
			else
			{
				// Smooth
				level->PreSmoother->Smooth(x, b, xEquals0);                               result.AddCost(level->PreSmoother->SolvingComputationalWork());
				// Residual computation
				r = b - A * x;                                                            result.AddCost(Cost::DAXPY(A));
			}
			
			if (Utils::ProgramArgs.Actions.ExportMultigridIterationVectors)
				level->ExportVector(x, "it" + to_string(this->IterationCount) + "_sol_afterPreSmoothing");

			//auto flopPreSmoothAndRes = result.IterationComputationalWork() - flopBeforePreSmooth;
			//cout << "- flopPreSmoothAndRes = " << flopPreSmoothAndRes << endl;

			//------------------------------------------------//
			// Restriction of the residual on the coarse grid //
			//------------------------------------------------//

			auto flopBeforeRestrict = result.IterationComputationalWork();
			IntergridTransferTimer.Start();

			Vector rc = level->Restrict(r);                                               result.AddCost(level->RestrictCost());
			                                                                              this->IntergridTransferCost += level->RestrictCost();
			IntergridTransferTimer.Pause();
			auto flopRestrict = result.IterationComputationalWork() - flopBeforeRestrict;
			assert(flopRestrict == level->RestrictCost());
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

			// Prolongation
			IntergridTransferTimer.Start();

			x += level->Prolong(ec);                                                  result.AddCost(level->ProlongCost() + Cost::AddVec(x));
			                                                                          this->IntergridTransferCost += level->ProlongCost() + Cost::AddVec(x);
			IntergridTransferTimer.Pause();

			if (Utils::ProgramArgs.Actions.ExportMultigridIterationVectors)
				level->ExportVector(x, "it" + to_string(this->IterationCount) + "_sol_cgc");

			//----------------//
			// Post-smoothing //
			//----------------//

			assert(!xEquals0);

			//auto flopBeforePostSmooth = result.IterationComputationalWork();

			if (level->PostSmoother->CanOptimizeResidualComputation() && (computeResidual || computeAx))
			{
				// Smooth and compute residual
				level->PostSmoother->SmoothAndComputeResidualOrAx(x, b, xEquals0, computeResidual, computeAx);    result.AddCost(level->PostSmoother->SolvingComputationalWork());
				if (computeResidual)
					result.Residual = level->PostSmoother->Residual();
				if (computeAx)
					result.Ax = level->PostSmoother->Ax();
			}
			else
			{
				// Smooth
				level->PostSmoother->Smooth(x, b, xEquals0);                          result.AddCost(level->PostSmoother->SolvingComputationalWork());
				if (computeAx)
				{
					// Residual computation
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

			if (Utils::ProgramArgs.Actions.ExportMultigridIterationVectors)
				level->ExportVector(x, "it" + to_string(this->IterationCount) + "_sol_afterPostSmoothing");

			//auto flopPostSmooth = result.IterationComputationalWork() - flopBeforePostSmooth;
			//cout << "- flopPostSmooth = " << flopPostSmooth << endl;
		}
	}

private:
	Vector FCGForKCycle(Level* level, Vector& r, IterationResult& result)
	{
		//const SparseMatrix& A = *level->OperatorMatrix;
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
				os << "multiple pairwise aggregation [-cs mpa -coarsening-factor " << Utils::ProgramArgs.Solver.MG.CoarseningFactor << "]" << endl;
			else if (H_CS == H_CoarsStgy::AgglomerationCoarseningByFaceNeighbours)
				os << "agglomeration by face neighbours [-cs n]" << endl;
			else if (H_CS == H_CoarsStgy::MultipleAgglomerationCoarseningByFaceNeighbours)
				os << "multiple agglomeration by face neighbours [-cs mn -coarsening-factor " << Utils::ProgramArgs.Solver.MG.CoarseningFactor << "]" << endl;
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

			os << "\t" << "Face coarsening strategy: ";
			if (Utils::IsRefinementStrategy(H_CS) || H_CS == H_CoarsStgy::IndependentRemeshing)
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