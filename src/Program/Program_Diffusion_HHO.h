#pragma once
#include <iomanip>
#include "../ProgramArguments.h"
#include "../HHO/Diffusion_HHO.h"
#include "../TestCases/Diffusion/DiffTestCaseFactory.h"
#include "../Mesher/MeshFactory.h"
#include "../Solver/SolverFactory.h"
#include "../Utils/ExportModule.h"

template <int Dim>
class Program_Diffusion_HHO
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

		//------------------------//
		//       Unit tests       //
		//------------------------//

		if (args.Actions.UnitTests)
		{
			// Unit tests
			Triangle::Test();
			Quadrilateral::Test();
#ifdef CGAL_ENABLED
			Polygon::Test();
#endif
			Tetrahedron::Test();
			TriangleIn3D::Test();

			cout << "Sanity check..." << endl;
			mesh->SanityCheck();
			if (args.Discretization.N <= 2)
				cout << *mesh << endl << endl;

			/*if (args.Discretization.MeshCode.compare("gmsh_tri") == 0)
			{
				GMSHMesh<Dim>* gmshMesh = dynamic_cast<GMSHMesh<Dim>*>(mesh);
				gmshMesh->RenumberLikeMe();
			}*/
		}

		//---------------------------//
		//   Coarsening unit tests   //
		//---------------------------//

		if (args.Actions.UnitTests && Dim == 2)
		{
			if (args.Solver.SolverCode.compare("mg") == 0 || args.Solver.SolverCode.compare("cgmg") == 0 || args.Solver.SolverCode.compare("fcgmg") == 0)
			{
				if (!Utils::IsRefinementStrategy(args.Solver.MG.H_CS))
				{
					/*Mesh<Dim>* m = mesh;
					while (true)
					{
						cout << "Coarsening..." << endl;
						m->CoarsenMesh(args.Solver.MG.H_CS, args.Solver.MG.CoarseningFactor);
						m = m->CoarseMesh;
						m->ExportToMatlab2();
						if (!m)
							break;
						cout << "Sanity check..." << endl;
						m->SanityCheck();
					}*/

					cout << "Export..." << endl;
					mesh->ExportFacesToMatlab(args.OutputDirectory + "/fine.dat");
					mesh->ExportElementCentersToMatlab(args.OutputDirectory + "/elem_fine.m");

					// 1st coarsening
					cout << "Coarsening..." << endl;
					mesh->CoarsenMesh(args.Solver.MG.H_CS, args.Solver.MG.FaceCoarseningStgy, args.Solver.MG.CoarseningFactor);
					/*cout << "Export..." << endl;
					mesh->CoarseMesh->ExportFacesToMatlab(args.OutputDirectory + "/coarse1.dat");
					mesh->CoarseMesh->ExportElementCentersToMatlab(args.OutputDirectory + "/elem_coarse1.m");*/
					cout << "Sanity check..." << endl;
					mesh->SanityCheck();
					// 2nd coarsening
					cout << "Coarsening..." << endl;
					mesh->CoarseMesh->CoarsenMesh(args.Solver.MG.H_CS, args.Solver.MG.FaceCoarseningStgy, args.Solver.MG.CoarseningFactor);
					/*cout << "Export..." << endl;
					mesh->CoarseMesh->CoarseMesh->ExportFacesToMatlab(args.OutputDirectory + "/coarse2.dat");
					mesh->CoarseMesh->CoarseMesh->ExportElementCentersToMatlab(args.OutputDirectory + "/elem_coarse2.m");*/
					cout << "Sanity check..." << endl;
					mesh->CoarseMesh->SanityCheck();
					// 3rd coarsening
					cout << "Coarsening..." << endl;
					mesh->CoarseMesh->CoarseMesh->CoarsenMesh(args.Solver.MG.H_CS, args.Solver.MG.FaceCoarseningStgy, args.Solver.MG.CoarseningFactor);
					/*cout << "Export..." << endl;
					mesh->CoarseMesh->CoarseMesh->CoarseMesh->ExportFacesToMatlab(args.OutputDirectory + "/coarse3.dat");
					mesh->CoarseMesh->CoarseMesh->CoarseMesh->ExportElementCentersToMatlab(args.OutputDirectory + "/elem_coarse3.m");*/
					cout << "Sanity check..." << endl;
					mesh->CoarseMesh->CoarseMesh->SanityCheck();
					//cout << *mesh << endl << endl;
					//cout << "Coarse mesh" << endl << *(mesh->CoarseMesh) << endl << endl;
				}
				else
				{
					mesh->ExportFacesToMatlab(args.OutputDirectory + "/fine.dat");
					mesh->ExportElementCentersToMatlab(args.OutputDirectory + "/elem_fine.m");

					if (mesh->CoarseMesh)
					{
						mesh->CoarseMesh->ExportFacesToMatlab(args.OutputDirectory + "/coarse1.dat");
						mesh->CoarseMesh->ExportElementCentersToMatlab(args.OutputDirectory + "/elem_coarse1.m");
						if (mesh->CoarseMesh->CoarseMesh)
						{
							mesh->CoarseMesh->CoarseMesh->ExportFacesToMatlab(args.OutputDirectory + "/coarse2.dat");
							mesh->CoarseMesh->CoarseMesh->ExportElementCentersToMatlab(args.OutputDirectory + "/elem_coarse2.m");
						}
					}
				}
			}
		}

		// Export source
		if (args.Actions.ExportSourceToGMSH && args.Discretization.Mesher.compare("gmsh") == 0)
			dynamic_cast<GMSHMesh<Dim>*>(mesh)->ExportToGMSH(testCase->SourceFunction, args.OutputDirectory + "/source", "source");

		// Export exact solution
		if (args.Actions.ExportExactSolutionToGMSH && testCase->ExactSolution && args.Discretization.Mesher.compare("gmsh") == 0)
			dynamic_cast<GMSHMesh<Dim>*>(mesh)->ExportToGMSH(testCase->ExactSolution, args.OutputDirectory + "/exsol", "exact solution");

		//----------------------//
		//       Assembly       //
		//----------------------//

		FunctionalBasis<Dim>* reconstructionBasis = new FunctionalBasis<Dim>(args.Discretization.ElemBasisCode, args.Discretization.PolyDegree, args.Discretization.UsePolynomialSpaceQ);
		FunctionalBasis<Dim>* cellBasis = new FunctionalBasis<Dim>(args.Discretization.ElemBasisCode, args.Discretization.PolyDegree - 1, args.Discretization.UsePolynomialSpaceQ);
		FunctionalBasis<Dim - 1>* faceBasis = new FunctionalBasis<Dim - 1>(args.Discretization.FaceBasisCode, args.Discretization.PolyDegree - 1, args.Discretization.UsePolynomialSpaceQ);

		HHOParameters<Dim>* hho = new HHOParameters<Dim>(mesh, args.Discretization.Stabilization, reconstructionBasis, cellBasis, faceBasis, args.Discretization.OrthogonalizeElemBasesCode, args.Discretization.OrthogonalizeFaceBasesCode);

		bool saveMatrixBlocks = args.Solver.SolverCode.compare("uamg") == 0 || args.Solver.SolverCode.compare("fcguamg") == 0;
		Diffusion_HHO<Dim>* problem = new Diffusion_HHO<Dim>(mesh, testCase, hho, args.Discretization.StaticCondensation, saveMatrixBlocks);


		// If full Neumann, check compatibility condition
		if (testCase->BC.Type == PbBoundaryConditions::FullNeumann)
		{
			// Compatibility condition: (f|1) + <neumann|1> = 0
			double integralF = problem->IntegralOverDomain(testCase->SourceFunction);
			double integralN = problem->IntegralOverBoundary(testCase->BC.NeumannFunction);
			if (abs(integralF + integralN) >= Utils::Eps)
				Utils::Error("Compatibility condition not respected: (f|1) + <neumann|1> = " + to_string(integralF + integralN));
		}

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
			int blockSizeForBlockSolver = args.Solver.BlockSize != -1 ? args.Solver.BlockSize : faceBasis->Size();

			Solver* solver = SolverFactory<Dim>::CreateSolver(args, problem, blockSizeForBlockSolver, out);

			Timer setupTimer;
			Timer solvingTimer;
			Timer totalTimer;
			totalTimer.Start();

			cout << "Solver: " << *solver << endl << endl;

			IterativeSolver* iterativeSolver = dynamic_cast<IterativeSolver*>(solver);
			if (iterativeSolver)
			{
				iterativeSolver->ComputeExactSolution = Utils::ProgramArgs.Actions.ExportErrorToGMSH || problem->A.rows() <= 2000;

				setupTimer.Start();
				if (Utils::ProgramArgs.Solver.SolverCode.compare("uamg") == 0 || Utils::ProgramArgs.Solver.SolverCode.compare("fcguamg") == 0)
					iterativeSolver->Setup(problem->A, problem->A_T_T, problem->A_T_ndF, problem->A_ndF_ndF);
				else
					solver->Setup(problem->A);
				setupTimer.Stop();

				cout << "Solving..." << endl;
				solvingTimer.Start();
				problem->SystemSolution = iterativeSolver->Solve(problem->b, args.Solver.InitialGuessCode);
				solvingTimer.Stop();
				cout << iterativeSolver->IterationCount << " iterations." << endl << endl;

				Multigrid* mg = dynamic_cast<Multigrid*>(iterativeSolver);
				if (mg)
				{
					int sizeTime = 7;
					int sizeWork = 5;

					double totalTime = solvingTimer.CPU().InMilliseconds;
					MFlops totalWork = iterativeSolver->SolvingComputationalWork;
					cout << "\t                             | CPU time |   Work " << endl;
					cout << "\t-------------------------------------------------" << endl;
					cout << "\tSmoothing and res. computing | " << setw(sizeTime) << (int)round((totalTime - mg->IntergridTransferTimer.CPU().InMilliseconds - mg->CoarseSolverTimer.CPU().InMilliseconds) / totalTime * 100) << "% | " << setw(sizeWork) << (int)round((totalWork - mg->IntergridTransferCost - mg->CoarseSolverCost) / totalWork * 100) << "%" << endl;
					cout << "\tIntergrid transfers          | " << setw(sizeTime) << (int)round(mg->IntergridTransferTimer.CPU().InMilliseconds / totalTime * 100) << "% | " << setw(sizeWork) << (int)round(mg->IntergridTransferCost / totalWork * 100) << "%" << endl;
					cout << "\tCoarse solver                | " << setw(sizeTime) << (int)round(mg->CoarseSolverTimer.CPU().InMilliseconds / totalTime * 100) << "% | " << setw(sizeWork) << (int)round(mg->CoarseSolverCost / totalWork * 100) << "%" << endl;
					cout << endl << endl;
				}
			}
			else
			{
				setupTimer.Start();
				solver->Setup(problem->A);
				setupTimer.Stop();

				cout << "Solving..." << endl;
				solvingTimer.Start();
				problem->SystemSolution = solver->Solve(problem->b);
				solvingTimer.Stop();
				cout << endl;
			}
			totalTimer.Stop();

			SolverFactory<Dim>::PrintStats(solver, setupTimer, solvingTimer, totalTimer);

			if (args.Actions.ExportErrorToGMSH && iterativeSolver)
				problem->ExportErrorToGMSH(iterativeSolver->ExactSolution - problem->SystemSolution, out);

			delete solver;


			if (args.Actions.ExportSolutionVectors || args.Actions.ExportSolutionToGMSH || testCase->ExactSolution || testCase->BC.Type == PbBoundaryConditions::FullNeumann)
			{
				cout << "----------------------------------------------------------" << endl;
				cout << "-                     Post-processing                    -" << endl;
				cout << "----------------------------------------------------------" << endl;

				problem->ReconstructHigherOrderApproximation();
			}

			//-----------------------------//
			//       Solution export       //
			//-----------------------------//

			if (args.Actions.ExportSolutionVectors)
			{
				if (args.Discretization.StaticCondensation)
					out.ExportVector(problem->SystemSolution, "solutionFaces");
				out.ExportVector(problem->GlobalHybridSolution, "solutionHybrid");
				out.ExportVector(problem->ReconstructedSolution, "solutionHigherOrder");
			}

			if (args.Actions.ExportMeshToMatlab)
			{
				mesh->ExportToMatlab(args.OutputDirectory);
				mesh->ExportToMatlab2(args.OutputDirectory + "/mesh.m");
			}

			if (args.Actions.ExportSolutionToGMSH && args.Discretization.Mesher.compare("gmsh") == 0)
				problem->ExportSolutionToGMSH(out);

			//----------------------//
			//       L2 error       //
			//----------------------//

			if (testCase->ExactSolution)
			{
				double error = problem->L2Error(testCase->ExactSolution);
				cout << endl << "L2 Error = " << std::scientific << error << endl;
				problem->AssertSchemeConvergence(error);
			}

			// Check mean value if full Neumann conditions
			if (testCase->BC.Type == PbBoundaryConditions::FullNeumann)
			{
				double meanValue = problem->MeanValueFromReconstructedCoeffs(problem->ReconstructedSolution);
				cout << "Mean value = " << meanValue << endl;
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
		delete problem;
		delete testCase;
		GaussLegendre::Free();
	}
};