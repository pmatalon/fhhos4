#pragma once
#include "../Utils/Timer.h"
#include "AlgebraicSolverFactory.h"
#include "LibraryUAMG.h"
#include "Multigrid/MultigridForHHO/MultigridForHHO.h"
#include "Multigrid/MultigridForHHO/P_MultigridForHHO.h"

template <int Dim>
class SolverFactory
{
public:
	static Solver* CreateSolver(const ProgramArguments& args, int blockSize, const ExportModule& out)
	{
		return AlgebraicSolverFactory::CreateSolver(args, Dim, blockSize, out, CoarseSolvers(out));
	}

	static Solver* CreateSolver(const ProgramArguments& args, Diffusion_HHO<Dim>* problem, int blockSize, const ExportModule& out)
	{
		Solver* solver = nullptr;
		bool libraryFCG = args.Solver.SolverCode.compare("fcglibuamg") == 0 || (args.Solver.SolverCode.compare("fcg") == 0 && args.Solver.PreconditionerCode.compare("libuamg") == 0);
		if (!libraryFCG && (args.Solver.PreconditionerCode.compare("libuamg") == 0 || (args.Solver.SolverCode.rfind("cg", 0) == 0 && args.Solver.SolverCode.find("libuamg") != string::npos)))
			Utils::FatalError("libuamg can only be used alone or with -s fcg (its FCG is inside the library): U-AMG is a non-linear preconditioner.");
		if (args.Solver.SolverCode.compare("libuamg") == 0 || libraryFCG)
		{
			if (!problem)
				Utils::FatalError("The solver libuamg needs the HHO problem: it cannot be a coarse solver.");
			solver = CreateLibraryUAMG(args, Dim, problem->HHO->FaceBasis->GetDegree(), problem->HHO->CellBasis->GetDegree(), libraryFCG);
		}
		else if (args.Solver.SolverCode.compare("cg") == 0)
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
			solver = AlgebraicSolverFactory::CreateUncondensedAMG(args, Dim, problem->HHO->FaceBasis->GetDegree(), problem->HHO->nCellUnknowns, problem->HHO->nFaceUnknowns, blockSize, out, CoarseSolvers(out));
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
		AlgebraicSolverFactory::SetMultigridParameters(mg, args, Dim, blockSize, out, CoarseSolvers(out));
	}

	// The coarse solvers of the multigrids: all the solvers, without the problem
	static AlgebraicSolverFactory::CoarseSolverFactory CoarseSolvers(const ExportModule& out)
	{
		return [&out](const ProgramArguments& args) { return CreateSolver(args, nullptr, 1, out); };
	}

public:
	static void PrintStats(Solver* solver, const Timer& setupTimer, const Timer& solvingTimer, const Timer& totalTimer)
	{
		IterativeSolver* iterativeSolver = dynamic_cast<IterativeSolver*>(solver);

		int sizeTime = 12;
		int sizeWork = 8;
		int sizeMatVec = 8;

		MFlops oneFineMatVec = 1;
		if (solver->Matrix)
			oneFineMatVec = Cost::MatVec(*solver->Matrix) * 1e-6;

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