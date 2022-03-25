#pragma once
#include <iomanip>
#include "../ProgramArguments.h"
#include "../FEM/Diffusion_FEM.h"
#include "../TestCases/Diffusion/DiffTestCaseFactory.h"
#include "../Mesher/MeshFactory.h"
#include "../Solver/SolverFactory.h"

template <int Dim>
class Program_Diffusion_FEM
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

		if (args.Discretization.MeshCode.compare("tri") != 0 && args.Discretization.MeshCode.compare("stri") != 0)
			Utils::FatalError("In FEM, only triangular elements are implemented.");

		Mesh<Dim>* mesh = MeshFactory<Dim>::BuildMesh(args, testCase);
		if (args.Discretization.Mesher.compare("gmsh") == 0)
			GMSHMesh<Dim>::CloseGMSH();

		cout << "Mesh storage > " << Utils::MemoryString(mesh->MemoryUsage()) << endl;

		mesh->SetDiffusionField(&testCase->DiffField);

		mesh->FillBoundaryAndInteriorVertexLists();
		mesh->SetBoundaryConditions(&testCase->BC);
		mesh->FillDirichletAndNeumannVertexLists();

		ExportModule out(args.OutputDirectory);

		// Export source
		if (args.Actions.Export.SourceToGMSH && args.Discretization.Mesher.compare("gmsh") == 0)
			dynamic_cast<GMSHMesh<Dim>*>(mesh)->ExportToGMSH_Nodes(testCase->SourceFunction, args.OutputDirectory + "/source", "source");

		// Export exact solution
		if (args.Actions.Export.ExactSolutionToGMSH && testCase->ExactSolution && args.Discretization.Mesher.compare("gmsh") == 0)
			dynamic_cast<GMSHMesh<Dim>*>(mesh)->ExportToGMSH_Nodes(testCase->ExactSolution, args.OutputDirectory + "/exsol", "exact solution");

		//----------------------//
		//       Assembly       //
		//----------------------//

		FunctionalBasis<Dim>* basis = new FunctionalBasis<Dim>(args.Discretization.ElemBasisCode, args.Discretization.PolyDegree, args.Discretization.UsePolynomialSpaceQ);
		Diffusion_FEM<Dim>* problem = new Diffusion_FEM<Dim>(mesh, testCase, basis);

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
			int blockSizeForBlockSolver = args.Solver.BlockSize != -1 ? args.Solver.BlockSize : 1;
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

			Vector systemSolution;

			if (iterativeSolver)
			{
				systemSolution = iterativeSolver->Solve(problem->b, args.Solver.InitialGuessCode);
				cout << iterativeSolver->IterationCount << " iterations." << endl << endl;
			}
			else
			{
				systemSolution = solver->Solve(problem->b);
				cout << endl;
			}
			solvingTimer.Stop();
			totalTimer.Stop();

			SolverFactory<Dim>::PrintStats(solver, setupTimer, solvingTimer, totalTimer);

			// Export algebraic error
			if (args.Actions.Export.ErrorToGMSH && iterativeSolver)
				mesh->ExportToGMSH_Nodes(iterativeSolver->ExactSolution - systemSolution, out.GetFilePathPrefix(), "error");

			delete solver;

			//-----------------------------//
			//       Solution export       //
			//-----------------------------//

			Vector solution = problem->AddDirichletValues(systemSolution);

			if (args.Actions.Export.SolutionVectors)
				out.ExportVector(solution, "solution");
			
			if (args.Actions.Export.SolutionToGMSH && args.Discretization.Mesher.compare("gmsh") == 0)
				mesh->ExportToGMSH_Nodes(solution, out.GetFilePathPrefix(), "potential");

			//----------------------//
			//       L2 error       //
			//----------------------//

			if (testCase->ExactSolution)
			{
				double error = problem->L2Error(testCase->ExactSolution, solution);
				cout << endl << "L2 Error = " << std::scientific << error << endl;
				problem->AssertSchemeConvergence(error);
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