#pragma once
#include "../Utils/Timer.h"

#include "Direct/EigenLU.h"
#include "Krylov/ConjugateGradient.h"
#include "Krylov/FlexibleConjugateGradient.h"
#include "Multigrid/MultigridForHHO/MultigridForHHO.h"
#include "Multigrid/MultigridForHHO/P_MultigridForHHO.h"
#include "Multigrid/UncondensedAMG/UncondensedAMG.h"
#include "Multigrid/AggregAMG/AggregAMG.h"
#include "FixedPoint/BlockJacobi.h"
#include "Krylov/EigenCG.h"
#include "Multigrid/AggregAMG/HighOrderAggregAMG.h"
#include "Multigrid/AGMG.h"

template <int Dim>
class SolverFactory
{
public:
	static Solver* CreateSolver(const ProgramArguments& args, int blockSize, const ExportModule& out)
	{
		Solver* solver = nullptr;
		string GaussSeidelRelaxation1 = "The relaxation parameter of Gauss-Seidel is 1. Delete -relax arg to remove this warning, or use -s sor.";

		if (args.Solver.SolverCode.compare("lu") == 0)
			solver = new EigenSparseLU();
		else if (args.Solver.SolverCode.compare("ch") == 0)
			solver = new EigenSparseCholesky();
		else if (args.Solver.SolverCode.compare("j") == 0)
			solver = new BlockJacobi(1, args.Solver.RelaxationParameter);
		else if (args.Solver.SolverCode.compare("gs") == 0)
		{
			if (args.Solver.RelaxationParameter != 1)
				Utils::Warning(GaussSeidelRelaxation1);
			solver = new GaussSeidel(Direction::Forward);
		}
		else if (args.Solver.SolverCode.compare("rgs") == 0)
		{
			if (args.Solver.RelaxationParameter != 1)
				Utils::Warning(GaussSeidelRelaxation1);
			solver = new GaussSeidel(Direction::Backward);
		}
		else if (args.Solver.SolverCode.compare("sgs") == 0)
		{
			if (args.Solver.RelaxationParameter != 1)
				Utils::Warning(GaussSeidelRelaxation1);
			solver = new GaussSeidel(Direction::Symmetric);
		}
		else if (args.Solver.SolverCode.compare("sor") == 0)
			solver = new BlockSOR(1, args.Solver.RelaxationParameter, Direction::Forward);
		else if (args.Solver.SolverCode.compare("rsor") == 0)
			solver = new BlockSOR(1, args.Solver.RelaxationParameter, Direction::Backward);
		else if (args.Solver.SolverCode.compare("ssor") == 0)
			solver = new BlockSOR(1, args.Solver.RelaxationParameter, Direction::Symmetric);
		else if (args.Solver.SolverCode.compare("bj") == 0)
			solver = new BlockJacobi(blockSize, args.Solver.RelaxationParameter);
		else if (args.Solver.SolverCode.compare("bj23") == 0)
			solver = new BlockJacobi(blockSize, 2.0 / 3.0);
		else if (args.Solver.SolverCode.compare("bgs") == 0)
		{
			if (args.Solver.RelaxationParameter != 1)
				Utils::Warning(GaussSeidelRelaxation1);
			solver = new BlockSOR(blockSize, 1, Direction::Forward);
		}
		else if (args.Solver.SolverCode.compare("rbgs") == 0)
		{
			if (args.Solver.RelaxationParameter != 1)
				Utils::Warning(GaussSeidelRelaxation1);
			solver = new BlockSOR(blockSize, 1, Direction::Backward);
		}
		else if (args.Solver.SolverCode.compare("sbgs") == 0)
		{
			if (args.Solver.RelaxationParameter != 1)
				Utils::Warning(GaussSeidelRelaxation1);
			solver = new BlockSOR(blockSize, 1, Direction::Symmetric);
		}
		else if (args.Solver.SolverCode.compare("bsor") == 0)
			solver = new BlockSOR(blockSize, args.Solver.RelaxationParameter, Direction::Forward);
		else if (args.Solver.SolverCode.compare("rbsor") == 0)
			solver = new BlockSOR(blockSize, args.Solver.RelaxationParameter, Direction::Backward);
		else if (args.Solver.SolverCode.compare("sbsor") == 0)
			solver = new BlockSOR(blockSize, args.Solver.RelaxationParameter, Direction::Symmetric); 
		else if (args.Solver.SolverCode.compare("cg") == 0)
			solver = new ConjugateGradient();
		else if (args.Solver.SolverCode.compare("eigencg") == 0)
			solver = new EigenCG();
		else if (args.Solver.SolverCode.compare("agmg") == 0)
			solver = new AGMG();
		else if (args.Solver.SolverCode.compare("aggregamg") == 0)
		{
			AggregAMG* mg = new AggregAMG(blockSize, 0.25, args.Solver.MG.Levels);
			SetMultigridParameters(mg, args, blockSize, out);
			mg->UseGalerkinOperator = 1;
			solver = mg;
		}
		else if (args.Solver.SolverCode.compare("hoaggregamg") == 0)
		{
			HighOrderAggregAMG* mg = new HighOrderAggregAMG(blockSize, 0.25, 2.0 / 3.0, 2);
			SetMultigridParameters(mg, args, blockSize, out);
			mg->UseGalerkinOperator = 1;
			solver = mg;
		}
		else
			Utils::FatalError("Unknown solver '" + args.Solver.SolverCode + "' or not applicable.");

		IterativeSolver* iterativeSolver = dynamic_cast<IterativeSolver*>(solver);
		if (iterativeSolver)
		{
			iterativeSolver->StoppingCrit = args.Solver.StoppingCrit;
			iterativeSolver->Tolerance = args.Solver.Tolerance;
			iterativeSolver->StagnationConvRate = args.Solver.StagnationConvRate;
			iterativeSolver->MaxIterations = args.Solver.MaxIterations;
			iterativeSolver->PrintIterationResults = args.Solver.PrintIterationResults;
		}

		return solver;
	}




	static Solver* CreateSolver(const ProgramArguments& args, Diffusion_HHO<Dim>* problem, int blockSize, const ExportModule& out)
	{
		Solver* solver = nullptr;
		if (args.Solver.SolverCode.compare("cg") == 0)
		{
			ConjugateGradient* cg = new ConjugateGradient();
			if (!args.Solver.PreconditionerCode.empty() && args.Solver.PreconditionerCode.compare("default") != 0)
			{
				ProgramArguments precondArgs = args;
				precondArgs.Solver.SolverCode = args.Solver.PreconditionerCode;
				precondArgs.Solver.PreconditionerCode = "";
				Solver* precondSolver = CreateSolver(precondArgs, problem, blockSize, out);
				cg->Precond = SolverPreconditioner(precondSolver);
			}
			solver = cg;
		}
		else if (args.Solver.SolverCode.rfind("cg", 0) == 0) // if SolverCode starts with "cg"
		{
			ConjugateGradient* cg = new ConjugateGradient();
			string preconditionerCode = args.Solver.SolverCode.substr(2, args.Solver.SolverCode.length() - 2);
			if (!preconditionerCode.empty())
			{
				ProgramArguments precondArgs = args;
				precondArgs.Solver.SolverCode = preconditionerCode;
				precondArgs.Solver.PreconditionerCode = "";
				Solver* precondSolver = CreateSolver(precondArgs, problem, blockSize, out);
				cg->Precond = SolverPreconditioner(precondSolver);
			}
			solver = cg;
		}
		else if (args.Solver.SolverCode.compare("fcg") == 0)
		{
			FlexibleConjugateGradient* fcg = new FlexibleConjugateGradient(1);
			if (!args.Solver.PreconditionerCode.empty() && args.Solver.PreconditionerCode.compare("default") != 0)
			{
				ProgramArguments precondArgs = args;
				precondArgs.Solver.SolverCode = args.Solver.PreconditionerCode;
				precondArgs.Solver.PreconditionerCode = "";
				Solver* precondSolver = CreateSolver(precondArgs, problem, blockSize, out);
				fcg->Precond = SolverPreconditioner(precondSolver);
			}
			solver = fcg;
		}
		else if (args.Solver.SolverCode.rfind("fcg", 0) == 0) // if SolverCode starts with "fcg"
		{
			FlexibleConjugateGradient* fcg = new FlexibleConjugateGradient(1);
			string preconditionerCode = args.Solver.SolverCode.substr(3, args.Solver.SolverCode.length() - 3);
			if (!preconditionerCode.empty())
			{
				ProgramArguments precondArgs = args;
				precondArgs.Solver.SolverCode = preconditionerCode;
				precondArgs.Solver.PreconditionerCode = "";
				Solver* precondSolver = CreateSolver(precondArgs, problem, blockSize, out);
				fcg->Precond = SolverPreconditioner(precondSolver);
			}
			solver = fcg;
		}
		else if (args.Solver.SolverCode.compare("mg") == 0)
		{
			if (args.Discretization.StaticCondensation)
			{
				MultigridForHHO<Dim>* mg = new MultigridForHHO<Dim>(args.Solver.MG.Levels);
				mg->UseHigherOrderReconstruction = args.Solver.MG.UseHigherOrderReconstruction;
				mg->H_Prolongation = args.Solver.MG.GMG_H_Prolong;
				mg->P_Prolongation = args.Solver.MG.GMG_P_Prolong;
				mg->P_Restriction = args.Solver.MG.GMG_P_Restrict;
				mg->UseHeterogeneousWeighting = args.Solver.MG.UseHeterogeneousWeighting;
				if (problem)
					mg->InitializeWithProblem(problem);
				SetMultigridParameters(mg, args, blockSize, out);
				solver = mg;
			}
			else
				Utils::FatalError("The Multigrid for HHO (-s mg) only applicable on HHO discretization with static condensation.");
		}
		else if (args.Solver.SolverCode.compare("p_mg") == 0)
		{
			if (args.Discretization.StaticCondensation)
			{
				P_MultigridForHHO<Dim>* mg = new P_MultigridForHHO<Dim>(problem);
				SetMultigridParameters(mg, args, blockSize, out);
				mg->HP_CS = HP_CoarsStgy::P_only;
				solver = mg;
			}
			else
				Utils::FatalError("The Multigrid for HHO only applicable on HHO discretization with static condensation.");
		}
		else if (args.Solver.SolverCode.compare("uamg") == 0)
		{
			UncondensedAMG* mg = new UncondensedAMG(Dim, problem->HHO->FaceBasis->GetDegree(), problem->HHO->nCellUnknowns, problem->HHO->nFaceUnknowns, 0.25, args.Solver.MG.UAMGFaceProlong, args.Solver.MG.UAMGCoarseningProlong, args.Solver.MG.UAMGMultigridProlong, args.Solver.MG.Levels);
			SetMultigridParameters(mg, args, blockSize, out);
			mg->CoarsePolyDegree = 0;
			solver = mg;
		}
		else
			solver = CreateSolver(args, blockSize, out);


		IterativeSolver* iterativeSolver = dynamic_cast<IterativeSolver*>(solver);
		if (iterativeSolver)
		{
			iterativeSolver->Tolerance = args.Solver.Tolerance;
			iterativeSolver->MaxIterations = args.Solver.MaxIterations;
			iterativeSolver->PrintIterationResults = args.Solver.PrintIterationResults;
		}

		return solver;
	}

private:
	static void SetMultigridParameters(Multigrid* mg, const ProgramArguments& args, int blockSize, const ExportModule& out)
	{
		mg->Out = ExportModule(out);
		mg->MatrixMaxSizeForCoarsestLevel = args.Solver.MG.MatrixMaxSizeForCoarsestLevel;
		mg->Cycle = args.Solver.MG.CycleLetter;
		mg->WLoops = args.Solver.MG.WLoops;
		mg->UseGalerkinOperator = args.Solver.MG.UseGalerkinOperator;
		mg->PreSmootherCode = args.Solver.MG.PreSmootherCode;
		mg->PostSmootherCode = args.Solver.MG.PostSmootherCode;
		mg->PreSmoothingIterations = args.Solver.MG.PreSmoothingIterations;
		mg->PostSmoothingIterations = args.Solver.MG.PostSmoothingIterations;
		mg->RelaxationParameter = args.Solver.RelaxationParameter;
		mg->BlockSizeForBlockSmoothers = blockSize;
		mg->CoarseLevelChangeSmoothingCoeff = args.Solver.MG.CoarseLevelChangeSmoothingCoeff;
		mg->CoarseLevelChangeSmoothingOperator = args.Solver.MG.CoarseLevelChangeSmoothingOperator;
		mg->HP_CS = args.Solver.MG.HP_CS;
		mg->H_CS = args.Solver.MG.H_CS;
		mg->P_CS = args.Solver.MG.P_CS;
		mg->FaceCoarseningStgy = args.Solver.MG.FaceCoarseningStgy;
		mg->NumberOfMeshes = args.Solver.MG.NumberOfMeshes;
		mg->CoarseningFactor = args.Solver.MG.CoarseningFactor;
		mg->ExportComponents = args.Actions.Export.MultigridComponents;

		// Coarse solver
		ProgramArguments argsCoarseSolver;
		argsCoarseSolver.Solver.SolverCode = args.Solver.MG.CoarseSolverCode;
		argsCoarseSolver.Solver.Tolerance = args.Solver.Tolerance;
		argsCoarseSolver.Solver.PrintIterationResults = false;
		if (Utils::EndsWith(args.Solver.MG.CoarseSolverCode, "aggregamg"))
		{
			argsCoarseSolver.Solver.MG.H_CS = H_CoarsStgy::AgglomerationCoarseningByFaceNeighbours;
			argsCoarseSolver.Solver.MG.CycleLetter = 'K';
		}
		else if (args.Solver.MG.CoarseSolverCode.compare("mg") == 0 || args.Solver.MG.CoarseSolverCode.compare("fcgmg") == 0)
		{
			argsCoarseSolver.Solver.MaxIterations = 1;
			argsCoarseSolver.Solver.MG.GMG_H_Prolong = args.Solver.MG.GMG_H_Prolong;
			argsCoarseSolver.Solver.MG.FaceProlongationCode = args.Solver.MG.FaceProlongationCode;
			argsCoarseSolver.Solver.MG.H_CS = args.Solver.MG.H_CS;
			argsCoarseSolver.Solver.MG.FaceCoarseningStgy = args.Solver.MG.FaceCoarseningStgy;
			argsCoarseSolver.Solver.MG.PreSmoothingIterations = 0;
			argsCoarseSolver.Solver.MG.PostSmoothingIterations = Dim == 2 ? 3 : 6;
		}
		mg->CoarseSolver = CreateSolver(argsCoarseSolver, nullptr, 1, out);
	}

public:
	static void PrintStats(Solver* solver, const Timer& setupTimer, const Timer& solvingTimer, const Timer& totalTimer)
	{
		IterativeSolver* iterativeSolver = dynamic_cast<IterativeSolver*>(solver);

		int sizeTime = 12;
		int sizeWork = 8;
		int sizeMatVec = 8;

		MFlops oneFineMatVec = Cost::MatVec(*solver->Matrix) * 1e-6;

		cout << "        |   CPU time   | Elapsed time ";
		if (iterativeSolver != nullptr)
			cout << "|  MFlops  |  MatVec  ";
		cout << endl;
		cout << "---------------------------------------";
		if (iterativeSolver != nullptr)
			cout << "---------------------";
		cout << endl;

		cout << "Setup   | " << setw(sizeTime) << setupTimer.CPU() << " | " << setw(sizeTime) << setupTimer.Elapsed();
		if (iterativeSolver != nullptr)
			cout << " | " << setw(sizeWork) << (int)round(iterativeSolver->SetupComputationalWork) << " | " << setw(sizeMatVec) << (int)round(iterativeSolver->SetupComputationalWork / oneFineMatVec);
		cout << endl;
		cout << "        | " << setw(sizeTime - 2) << setupTimer.CPU().InSeconds() << " s | " << setw(sizeTime - 2) << setupTimer.Elapsed().InSeconds() << " s ";
		if (iterativeSolver != nullptr)
			cout << "| " << setw(sizeWork) << " " << " | " << setw(sizeMatVec);
		cout << endl;
		cout << "---------------------------------------";
		if (iterativeSolver != nullptr)
			cout << "---------------------";
		cout << endl;

		cout << "Solving | " << setw(sizeTime) << solvingTimer.CPU() << " | " << setw(sizeTime) << solvingTimer.Elapsed();
		if (iterativeSolver != nullptr)
			cout << " | " << setw(sizeWork) << (int)round(iterativeSolver->SolvingComputationalWork) << " | " << setw(sizeMatVec) << (int)round(iterativeSolver->SolvingComputationalWork / oneFineMatVec);
		cout << endl;
		cout << "        | " << setw(sizeTime - 2) << solvingTimer.CPU().InSeconds() << " s | " << setw(sizeTime - 2) << solvingTimer.Elapsed().InSeconds() << " s ";
		if (iterativeSolver != nullptr)
			cout << "| " << setw(sizeWork) << " " << " | " << setw(sizeMatVec);
		cout << endl;
		cout << "---------------------------------------";
		if (iterativeSolver != nullptr)
			cout << "---------------------";
		cout << endl;

		cout << "Total   | " << setw(sizeTime) << totalTimer.CPU() << " | " << setw(sizeTime) << totalTimer.Elapsed();
		if (iterativeSolver != nullptr)
			cout << " | " << setw(sizeWork) << (int)round((iterativeSolver->SetupComputationalWork + iterativeSolver->SolvingComputationalWork)) << " | " << setw(sizeMatVec) << (int)round((iterativeSolver->SetupComputationalWork + iterativeSolver->SolvingComputationalWork) / oneFineMatVec);
		cout << endl;
		cout << "        | " << setw(sizeTime - 2) << totalTimer.CPU().InSeconds() << " s | " << setw(sizeTime - 2) << totalTimer.Elapsed().InSeconds() << " s ";
		if (iterativeSolver != nullptr)
			cout << "| " << setw(sizeWork) << " " << " | " << setw(sizeMatVec);
		cout << endl;

		cout << endl;
	}
};