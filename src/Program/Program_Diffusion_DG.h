#pragma once
#include <iomanip>
#include "../ProgramArguments.h"
#include "../Discretizations/DG/Diffusion_DG.h"
#include "../TestCases/Diffusion/DiffTestCaseFactory.h"
#include "../Mesher/MeshFactory.h"
#include "../FunctionalBasis/FunctionalBasisFactory.h"
#include "../Solver/SolverFactory.h"

template <int Dim>
class Program_Diffusion_DG
{
public:
	static void Execute(ProgramArguments& args)
	{
		GaussLegendre::Init();

		cout << "-----------------------------------------------------------" << endl;
		cout << "-                         Problem                         -" << endl;
		cout << "-----------------------------------------------------------" << endl;

		//-------------------------------------------------------------------------------------------//
		//   Test case defining the source function, boundary conditions and diffusion coefficient   //
		//-------------------------------------------------------------------------------------------//

		DiffusionTestCase<Dim>* testCase = DiffTestCaseFactory<Dim>::Create(args.Problem);
		testCase->PrintPhysicalProblem();

		cout << endl;

		//----------//
		//   Mesh   //
		//----------//

		cout << "-----------------------------------------------------------" << endl;
		cout << "-                   Mesh construction                     -" << endl;
		cout << "-----------------------------------------------------------" << endl;

		Mesh<Dim>* mesh = MeshFactory<Dim>::BuildMesh(args, testCase);
		if (args.Discretization.Mesher.compare("gmsh") == 0)
			GMSHMesh<Dim>::CloseGMSH();

		cout << "Mesh storage > " << Utils::MemoryString(mesh->MemoryUsage()) << endl;

		mesh->SetDiffusionField(&testCase->DiffField);

		mesh->SetBoundaryConditions(&testCase->BC);

		ExportModule out(args.OutputDirectory);

		// Export source
		if (args.Actions.Export.SourceToGMSH && args.Discretization.Mesher.compare("gmsh") == 0)
			dynamic_cast<GMSHMesh<Dim>*>(mesh)->ExportToGMSH_Elements(testCase->SourceFunction, args.OutputDirectory + "/source", "source");

		//----------------------//
		//       Assembly       //
		//----------------------//

		FunctionalBasis<Dim>* basis = FunctionalBasisFactory<Dim>::Create(args.Discretization.ElemBasisCode, args.Discretization.PolyDegree, args.Discretization.UsePolynomialSpaceQ);
		Diffusion_DG<Dim>* problem = new Diffusion_DG<Dim>(mesh, testCase, args.OutputDirectory, basis, args.Discretization.PenalizationCoefficient);

		cout << endl;
		cout << "----------------------------------------------------------" << endl;
		cout << "-                       Assembly                         -" << endl;
		cout << "----------------------------------------------------------" << endl;
		Timer assemblyTimer;
		assemblyTimer.Start();

		problem->Assemble(args.Actions, out);
		cout << "System storage: " << Utils::MemoryString(Utils::MemoryUsage(problem->A) + Utils::MemoryUsage(problem->b)) << endl;

		assemblyTimer.Stop();
		cout << endl << "Assembly time: CPU = " << assemblyTimer.CPU() << ", elapsed = " << assemblyTimer.Elapsed() << endl;

		//------------------------------------//
		//       Linear system solution       //
		//------------------------------------//

		if (args.Actions.SolveLinearSystem)
		{
			cout << endl;
			cout << "----------------------------------------------------------" << endl;
			cout << "-                 Linear system solution                 -" << endl;
			cout << "----------------------------------------------------------" << endl;

			// Solver creation
			int blockSizeForBlockSolver = args.Solver.BlockSize != -1 ? args.Solver.BlockSize : problem->Basis->Size();
			Solver* solver = SolverFactory<Dim>::CreateSolver(args, blockSizeForBlockSolver, out);

			cout << "Solver: " << *solver << endl << endl;

			Timer setupTimer;
			Timer solvingTimer;
			Timer totalTimer;

			// Setup
			IterativeSolver* iterativeSolver = dynamic_cast<IterativeSolver*>(solver);
			if (iterativeSolver)
				iterativeSolver->ComputeExactSolution = Utils::ProgramArgs.Actions.Export.ErrorToGMSH || problem->A.rows() <= 2000;

			totalTimer.Start();

			setupTimer.Start();
			solver->Setup(problem->A);
			setupTimer.Stop();

			// Solving
			cout << "Solving..." << endl;

			solvingTimer.Start();
			if (iterativeSolver)
			{
				problem->SystemSolution = iterativeSolver->Solve(problem->b, args.Solver.InitialGuessCode);
				cout << iterativeSolver->IterationCount << " iterations." << endl << endl;
			}
			else
			{
				problem->SystemSolution = solver->Solve(problem->b);
				cout << endl;
			}
			solvingTimer.Stop();
			totalTimer.Stop();

			SolverFactory<Dim>::PrintStats(solver, setupTimer, solvingTimer, totalTimer);

			// Export algebraic error
			if (args.Actions.Export.ErrorToGMSH && iterativeSolver)
				mesh->ExportToGMSH_Elements(problem->Basis, iterativeSolver->ExactSolution - problem->SystemSolution, out.GetFilePathPrefix(), "error");

			delete solver;

			//-----------------------------//
			//       Solution export       //
			//-----------------------------//

			if (args.Actions.Export.SolutionVectors)
				out.ExportVector(problem->SystemSolution, "solution");
			
			if (args.Actions.Export.SolutionToGMSH && args.Discretization.Mesher.compare("gmsh") == 0)
				mesh->ExportToGMSH_Elements(problem->Basis, problem->SystemSolution, out.GetFilePathPrefix(), "potential");

			//----------------------//
			//       L2 error       //
			//----------------------//

			if (testCase->ExactSolution)
			{
				double error = problem->L2Error(testCase->ExactSolution);
				cout << endl << "L2 Error = " << std::scientific << error << endl;
			}
		}

		//--------------------------//
		//       Deallocation       //
		//--------------------------//

		delete mesh;
		delete basis;
		delete problem;
		delete testCase;
		GaussLegendre::Free();
	}
};