#pragma once
#include <iomanip>
#include "../ProgramArguments.h"
#include "../Discretizations/HHO/BiHarmonicDDMixedForm_HHO.h"
#include "../TestCases/BiHarmonicDD/BiHarDDTestCaseFactory.h"
#include "../Mesher/MeshFactory.h"
#include "../Solver/SolverFactory.h"
#include "../Utils/ExportModule.h"

// Bi-harmonic equation in mixed form with Dirichlet BC enforced on both Laplacian problems (hence the suffix DD)

template <int Dim>
class Program_BiHarmonicDD_HHO
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

		BiHarmonicDDTestCase<Dim>* testCase = BiHarDDTestCaseFactory<Dim>::Create(args.Problem);

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

		mesh->SetBoundaryConditions(&testCase->DirichletBC);

		// Export source
		if (args.Actions.Export.SourceToGMSH && args.Discretization.Mesher.compare("gmsh") == 0)
			dynamic_cast<GMSHMesh<Dim>*>(mesh)->ExportToGMSH(testCase->SourceFunction, args.OutputDirectory + "/source", "source");

		// Export exact solution
		if (args.Actions.Export.ExactSolutionToGMSH && testCase->ExactSolution && args.Discretization.Mesher.compare("gmsh") == 0)
			dynamic_cast<GMSHMesh<Dim>*>(mesh)->ExportToGMSH(testCase->ExactSolution, args.OutputDirectory + "/exsol", "exact solution");

		//----------------------//
		//       Assembly       //
		//----------------------//

		FunctionalBasis<Dim>* reconstructionBasis = new FunctionalBasis<Dim>(args.Discretization.ElemBasisCode, args.Discretization.PolyDegree, args.Discretization.UsePolynomialSpaceQ);
		FunctionalBasis<Dim>* cellBasis = new FunctionalBasis<Dim>(args.Discretization.ElemBasisCode, args.Discretization.PolyDegree - 1, args.Discretization.UsePolynomialSpaceQ);
		FunctionalBasis<Dim - 1>* faceBasis = new FunctionalBasis<Dim - 1>(args.Discretization.FaceBasisCode, args.Discretization.PolyDegree - 1, args.Discretization.UsePolynomialSpaceQ);

		HHOParameters<Dim>* hho = new HHOParameters<Dim>(mesh, args.Discretization.Stabilization, reconstructionBasis, cellBasis, faceBasis, args.Discretization.OrthogonalizeElemBasesCode, args.Discretization.OrthogonalizeFaceBasesCode);

		bool saveMatrixBlocks = args.Solver.SolverCode.compare("uamg") == 0 || args.Solver.SolverCode.compare("fcguamg") == 0;
		BiHarmonicDDMixedForm_HHO<Dim>* biHarPb = new BiHarmonicDDMixedForm_HHO<Dim>(mesh, testCase, hho, saveMatrixBlocks);

		cout << endl;
		cout << "----------------------------------------------------------" << endl;
		cout << "-          Assembly of 1st diffusion problem             -" << endl;
		cout << "----------------------------------------------------------" << endl;
		Timer assemblyTimer;
		assemblyTimer.Start();

		biHarPb->AssembleDiffPb1();
		cout << "System storage: " << Utils::MemoryString(Utils::MemoryUsage(biHarPb->DiffPb1().A) + Utils::MemoryUsage(biHarPb->DiffPb1().b)) << endl;

		assemblyTimer.Stop();
		cout << endl << "Assembly time: CPU = " << assemblyTimer.CPU() << ", elapsed = " << assemblyTimer.Elapsed() << endl;

		//------------------------------------//
		//       Linear system solution       //
		//------------------------------------//

		if (args.Actions.SolveLinearSystem)
		{
			Vector faceSolution;

			cout << endl;
			cout << "----------------------------------------------------------" << endl;
			cout << "-                   Solve linear system                  -" << endl;
			cout << "----------------------------------------------------------" << endl;

			int blockSizeForBlockSolver = args.Solver.BlockSize != -1 ? args.Solver.BlockSize : faceBasis->Size();

			// Solve 1st diffusion biHarPb
			Solver* solver = SolverFactory<Dim>::CreateSolver(args, &biHarPb->DiffPb1(), blockSizeForBlockSolver, out);

			Timer setupTimer;
			Timer solvingTimer1;
			Timer solvingTimer2;
			Timer totalTimer;
			totalTimer.Start();

			cout << "Solver: " << *solver << endl << endl;
			IterativeSolver* iterativeSolver = dynamic_cast<IterativeSolver*>(solver);
			if (iterativeSolver)
			{
				setupTimer.Start();
				if (Utils::ProgramArgs.Solver.SolverCode.compare("uamg") == 0 || Utils::ProgramArgs.Solver.SolverCode.compare("fcguamg") == 0)
					iterativeSolver->Setup(biHarPb->DiffPb1().A, biHarPb->DiffPb1().A_T_T, biHarPb->DiffPb1().A_T_ndF, biHarPb->DiffPb1().A_ndF_ndF);
				else
					solver->Setup(biHarPb->DiffPb1().A);
				setupTimer.Stop();

				cout << "Solving first problem..." << endl;
				solvingTimer1.Start();
				faceSolution = iterativeSolver->Solve(biHarPb->DiffPb1().b, args.Solver.InitialGuessCode);
				solvingTimer1.Stop();
			}
			else
			{
				setupTimer.Start();
				solver->Setup(biHarPb->DiffPb1().A);
				setupTimer.Stop();

				cout << "Solving first problem..." << endl;
				solvingTimer1.Start();
				faceSolution = solver->Solve(biHarPb->DiffPb1().b);
				solvingTimer1.Stop();
				cout << endl;
			}

			Vector reconstructedSolution = biHarPb->DiffPb1().ReconstructHigherOrderApproximationFromFaceCoeffs(faceSolution);

			//----------------------//
			//       L2 error       //
			//----------------------//

			if (testCase->MinusLaplacianOfSolution)
			{
				double error = biHarPb->DiffPb1().L2Error(testCase->MinusLaplacianOfSolution, reconstructedSolution);
				cout << endl << "L2 Error = " << std::scientific << error << endl;
			}

			cout << "----------------------------------------------------------" << endl;
			cout << "-                  2nd diffusion problem                 -" << endl;
			cout << "----------------------------------------------------------" << endl;

			//bool higherDegree = true;
			//if (higherDegree)
				biHarPb->AssembleDiffPb2(reconstructedSolution);
			//else
				//biHarPb->AssembleDiffPb2(biHarPb->DiffPb1().GlobalHybridSolution.head(biHarPb->DiffPb1().HHO->nTotalCellUnknowns));

			cout << "Solving second problem..." << endl;
			solvingTimer2.Start();
			faceSolution = solver->Solve(biHarPb->DiffPb2().b);
			solvingTimer2.Stop();
			cout << endl;

			totalTimer.Stop();

			delete solver;

			//-----------------------------//
			//       Solution export       //
			//-----------------------------//

			if (args.Actions.Export.SolutionVectors || args.Actions.Export.SolutionToGMSH || testCase->ExactSolution)
			{
				cout << "----------------------------------------------------------" << endl;
				cout << "-                     Post-processing                    -" << endl;
				cout << "----------------------------------------------------------" << endl;

				reconstructedSolution = biHarPb->DiffPb2().ReconstructHigherOrderApproximationFromFaceCoeffs(faceSolution);
				if (args.Actions.Export.SolutionVectors)
					out.ExportVector(reconstructedSolution, "solution");
			}

			if (args.Actions.Export.MeshToMatlab)
			{
				mesh->ExportToMatlab(args.OutputDirectory);
				mesh->ExportToMatlab2(args.OutputDirectory + "/mesh.m");
			}

			if (args.Actions.Export.SolutionToGMSH && args.Discretization.Mesher.compare("gmsh") == 0)
				biHarPb->DiffPb2().ExportSolutionToGMSH(reconstructedSolution, out);



			//----------------------//
			//       L2 error       //
			//----------------------//

			if (testCase->ExactSolution)
			{
				double error = biHarPb->DiffPb2().L2Error(testCase->ExactSolution, reconstructedSolution);
				cout << endl << "L2 Error = " << std::scientific << error << endl;
				//biHarPb->AssertSchemeConvergence(error);
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