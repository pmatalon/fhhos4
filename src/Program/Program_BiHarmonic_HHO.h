#pragma once
#include <iomanip>
#include "../ProgramArguments.h"
#include "../HHO/BiHarmonicMixedForm_HHO.h"
#include "../TestCases/BiHarmonic/BiHarTestCaseFactory.h"
#include "../Mesher/MeshFactory.h"
#include "../Solver/SolverFactory.h"
#include "../Solver/Krylov/BiHarmonicCG.h"
#include "../Solver/BiHarmonic/BiHarmonicGradientDescent.h"
#include "../Utils/ExportModule.h"

// Bi-harmonic equation in mixed form with mixed (homogeneous) Dirichel-Neumann BC

template <int Dim>
class Program_BiHarmonic_HHO
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

		BiHarmonicTestCase<Dim>* testCase = BiHarTestCaseFactory<Dim>::Create(args.Problem);

		testCase->PrintPhysicalProblem();

		cout << endl;

		ExportModule out(args.OutputDirectory);

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

		// Export source
		if (args.Actions.ExportSourceToGMSH && args.Discretization.Mesher.compare("gmsh") == 0)
			dynamic_cast<GMSHMesh<Dim>*>(mesh)->ExportToGMSH(testCase->SourceFunction, args.OutputDirectory + "/source", "source");

		// Export exact solution
		if (args.Actions.ExportExactSolutionToGMSH && testCase->ExactSolution && args.Discretization.Mesher.compare("gmsh") == 0)
			dynamic_cast<GMSHMesh<Dim>*>(mesh)->ExportToGMSH(testCase->ExactSolution, args.OutputDirectory + "/exsol", "exact solution");

		//----------------------//
		//       Assembly       //
		//----------------------//

		auto FullNeumann = BoundaryConditions::HomogeneousNeumannEverywhere();
		mesh->SetBoundaryConditions(&FullNeumann);

		FunctionalBasis<Dim>* reconstructionBasis = new FunctionalBasis<Dim>(args.Discretization.ElemBasisCode, args.Discretization.PolyDegree, args.Discretization.UsePolynomialSpaceQ);
		FunctionalBasis<Dim>* cellBasis = new FunctionalBasis<Dim>(args.Discretization.ElemBasisCode, args.Discretization.PolyDegree - 1, args.Discretization.UsePolynomialSpaceQ);
		FunctionalBasis<Dim - 1>* faceBasis = new FunctionalBasis<Dim - 1>(args.Discretization.FaceBasisCode, args.Discretization.PolyDegree - 1, args.Discretization.UsePolynomialSpaceQ);

		HHOParameters<Dim>* hho = new HHOParameters<Dim>(mesh, args.Discretization.Stabilization, reconstructionBasis, cellBasis, faceBasis, args.Discretization.OrthogonalizeElemBasesCode, args.Discretization.OrthogonalizeFaceBasesCode);

		bool saveMatrixBlocks = args.Solver.SolverCode.compare("uamg") == 0 || args.Solver.SolverCode.compare("fcguamg") == 0;
		BiHarmonicMixedForm_HHO<Dim>* biHarPb = new BiHarmonicMixedForm_HHO<Dim>(mesh, testCase, hho, args.Solver.BiHarReconstructBoundary, args.Actions.EnforceDirichletBC, saveMatrixBlocks);

		cout << endl;
		cout << "------------------------------------------------------" << endl;
		cout << "-          Assembly of diffusion problem             -" << endl;
		cout << "------------------------------------------------------" << endl;
		Timer assemblyTimer;
		assemblyTimer.Start();

		biHarPb->AssembleDiffPb();

		if (args.Actions.EnforceDirichletBC)
		{
			Mesh<Dim>* meshLast = MeshFactory<Dim>::BuildMesh(args, testCase);
			if (args.Discretization.Mesher.compare("gmsh") == 0)
				GMSHMesh<Dim>::CloseGMSH();
			biHarPb->SetupLastPb(meshLast);
		}

		assemblyTimer.Stop();
		cout << endl << "Assembly time: CPU = " << assemblyTimer.CPU() << ", elapsed = " << assemblyTimer.Elapsed() << endl;

		//------------------------------------//
		//       Linear system solution       //
		//------------------------------------//

		if (args.Actions.SolveLinearSystem)
		{
			cout << endl;
			cout << "----------------------------------------------" << endl;
			cout << "-           Setup Laplacian solver           -" << endl;
			cout << "----------------------------------------------" << endl;

			int blockSizeForBlockSolver = args.Solver.BlockSize != -1 ? args.Solver.BlockSize : faceBasis->Size();

			// Solve 1st diffusion biHarPb
			Solver* diffSolver = SolverFactory<Dim>::CreateSolver(args, &biHarPb->DiffPb(), blockSizeForBlockSolver, out);

			Timer setupTimer;
			Timer solvingTimer1;
			Timer solvingTimer2;
			Timer totalTimer;
			totalTimer.Start();

			cout << "Solver: " << *diffSolver << endl << endl;
			IterativeSolver* iterativeSolver = dynamic_cast<IterativeSolver*>(diffSolver);
			if (iterativeSolver)
			{
				setupTimer.Start();
				if (Utils::ProgramArgs.Solver.SolverCode.compare("uamg") == 0 || Utils::ProgramArgs.Solver.SolverCode.compare("fcguamg") == 0)
					iterativeSolver->Setup(biHarPb->DiffPb().A, biHarPb->DiffPb().A_T_T, biHarPb->DiffPb().A_T_ndF, biHarPb->DiffPb().A_ndF_ndF);
				else
					diffSolver->Setup(biHarPb->DiffPb().A);
				setupTimer.Stop();
			}
			else
			{
				setupTimer.Start();
				diffSolver->Setup(biHarPb->DiffPb().A);
				setupTimer.Stop();
			}

			biHarPb->SetDiffSolver(diffSolver);


			if (args.Actions.EnforceDirichletBC)
			{
				cout << "Enforce Dirichlet BC to the solution, so setup of last solver..." << endl << endl;

				Solver* lastSolver = SolverFactory<Dim>::CreateSolver(args, biHarPb->LastPb(), blockSizeForBlockSolver, out);
				lastSolver->Setup(biHarPb->LastPb()->A);
				biHarPb->SetLastPbSolver(lastSolver);
			}


			cout << "-------------------------------------" << endl;
			cout << "-     Solve bi-harmonic problem     -" << endl;
			cout << "-------------------------------------" << endl;

			IterativeSolver* biHarSolver = nullptr;
			if (args.Solver.BiHarmonicSolverCode.compare("cg") == 0)
				biHarSolver = new BiHarmonicCG<Dim>(*biHarPb);
			else if (args.Solver.BiHarmonicSolverCode.compare("gd") == 0)
				biHarSolver = new BiHarmonicGradientDescent<Dim>(*biHarPb);
			else
				Utils::FatalError("Unknown bi-harmonic solver '" + args.Solver.BiHarmonicSolverCode + "'");

			biHarSolver->Tolerance = args.Solver.Tolerance;
			biHarSolver->MaxIterations = args.Solver.MaxIterations;
			// Compute L2-error at each iteration
			if ((args.Solver.ComputeIterL2Error || args.Actions.ExportIterationL2Errors) && testCase->ExactSolution)
			{
				biHarSolver->OnNewSolution = [&biHarPb, &testCase](IterationResult& result, const Vector& theta)
				{
					Vector reconstructedSolution = biHarPb->ComputeSolution(theta);
					result.L2Error = biHarPb->DiffPb().L2Error(testCase->ExactSolution, reconstructedSolution);
				};
			}
			// Export iteration results
			if (args.Actions.ExportIterationResiduals || args.Actions.ExportIterationL2Errors)
			{
				if (args.Actions.ExportIterationResiduals)
					out.CleanFile("iteration_residuals");
				if (args.Actions.ExportIterationL2Errors)
					out.CleanFile("iteration_l2errors");

				biHarSolver->OnIterationEnd = [&args, &out](const IterationResult& result)
				{
					if (args.Actions.ExportIterationResiduals)
						out.ExportNewVectorValue(result.NormalizedResidualNorm, "iteration_residuals");
					if (args.Actions.ExportIterationL2Errors && result.L2Error != -1)
						out.ExportNewVectorValue(result.L2Error, "iteration_l2errors");
				};
			}

			Vector theta = biHarSolver->Solve();

			delete biHarSolver;

			Vector reconstructedSolution = biHarPb->ComputeSolution(theta);

			//-----------------------------//
			//       Solution export       //
			//-----------------------------//

			if (args.Actions.ExportSolutionVectors)
				out.ExportVector(reconstructedSolution, "solutionHigherOrder");

			if (args.Actions.ExportMeshToMatlab)
			{
				mesh->ExportToMatlab(args.OutputDirectory);
				mesh->ExportToMatlab2(args.OutputDirectory + "/mesh.m");
			}

			if (args.Actions.ExportSolutionToGMSH && args.Discretization.Mesher.compare("gmsh") == 0)
				biHarPb->DiffPb().ExportSolutionToGMSH(reconstructedSolution, out);

			//----------------------//
			//       L2 error       //
			//----------------------//

			if (testCase->ExactSolution)
			{
				double error = biHarPb->DiffPb().L2Error(testCase->ExactSolution, reconstructedSolution);
				cout << endl << "L2 Error = " << std::scientific << error << endl;
			}
		}

		//--------------------------//
		//       Deallocation       //
		//--------------------------//

		delete mesh;
		delete cellBasis;
		delete faceBasis;
		delete reconstructionBasis;
		delete hho;
		delete biHarPb;
		delete testCase;
		GaussLegendre::Free();
	}
};