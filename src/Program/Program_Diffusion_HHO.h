#pragma once
#include "../ProgramArguments.h"
#include "../Discretizations/HHO/Diffusion_HHO.h"
#include "../TestCases/Diffusion/DiffTestCaseFactory.h"
#include "../Mesher/MeshFactory.h"
#include "../FunctionalBasis/FunctionalBasisFactory.h"
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

		if (args.Problem.ComputeNormalDerivative && testCase->ExactSolution_Neumann)
			mesh->RenumberBoundaryElementsFirst();

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
#ifdef ENABLE_3D
			Tetrahedron::Test();
			TriangleIn3D::Test();
#endif

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
			if (args.Solver.SolverCode.compare("mg") == 0 || args.Solver.PreconditionerCode.compare("mg") == 0)
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
					mesh->ExportElementNumbersToMatlab(args.OutputDirectory + "/elem_fine.m");

					// 1st coarsening
					cout << "Coarsening..." << endl;
					mesh->CoarsenMesh(args.Solver.MG.H_CS, args.Solver.MG.FaceCoarseningStgy, args.Solver.MG.BoundaryFaceCollapsing, args.Solver.MG.CoarseningFactor);
					/*cout << "Export..." << endl;
					mesh->CoarseMesh->ExportFacesToMatlab(args.OutputDirectory + "/coarse1.dat");
					mesh->CoarseMesh->ExportElementNumbersToMatlab(args.OutputDirectory + "/elem_coarse1.m");*/
					cout << "Sanity check..." << endl;
					mesh->SanityCheck();
					// 2nd coarsening
					cout << "Coarsening..." << endl;
					mesh->CoarseMesh->CoarsenMesh(args.Solver.MG.H_CS, args.Solver.MG.FaceCoarseningStgy, args.Solver.MG.BoundaryFaceCollapsing, args.Solver.MG.CoarseningFactor);
					/*cout << "Export..." << endl;
					mesh->CoarseMesh->CoarseMesh->ExportFacesToMatlab(args.OutputDirectory + "/coarse2.dat");
					mesh->CoarseMesh->CoarseMesh->ExportElementNumbersToMatlab(args.OutputDirectory + "/elem_coarse2.m");*/
					cout << "Sanity check..." << endl;
					mesh->CoarseMesh->SanityCheck();
					// 3rd coarsening
					cout << "Coarsening..." << endl;
					mesh->CoarseMesh->CoarseMesh->CoarsenMesh(args.Solver.MG.H_CS, args.Solver.MG.FaceCoarseningStgy, args.Solver.MG.BoundaryFaceCollapsing, args.Solver.MG.CoarseningFactor);
					/*cout << "Export..." << endl;
					mesh->CoarseMesh->CoarseMesh->CoarseMesh->ExportFacesToMatlab(args.OutputDirectory + "/coarse3.dat");
					mesh->CoarseMesh->CoarseMesh->CoarseMesh->ExportElementNumbersToMatlab(args.OutputDirectory + "/elem_coarse3.m");*/
					cout << "Sanity check..." << endl;
					mesh->CoarseMesh->CoarseMesh->SanityCheck();
					//cout << *mesh << endl << endl;
					//cout << "Coarse mesh" << endl << *(mesh->CoarseMesh) << endl << endl;
				}
				else
				{
					mesh->ExportFacesToMatlab(args.OutputDirectory + "/fine.dat");
					mesh->ExportElementNumbersToMatlab(args.OutputDirectory + "/elem_fine.m");

					if (mesh->CoarseMesh)
					{
						mesh->CoarseMesh->ExportFacesToMatlab(args.OutputDirectory + "/coarse1.dat");
						mesh->CoarseMesh->ExportElementNumbersToMatlab(args.OutputDirectory + "/elem_coarse1.m");
						if (mesh->CoarseMesh->CoarseMesh)
						{
							mesh->CoarseMesh->CoarseMesh->ExportFacesToMatlab(args.OutputDirectory + "/coarse2.dat");
							mesh->CoarseMesh->CoarseMesh->ExportElementNumbersToMatlab(args.OutputDirectory + "/elem_coarse2.m");
						}
					}
				}
			}
		}

		// Export source
		if (args.Actions.Export.SourceToGMSH && args.Discretization.Mesher.compare("gmsh") == 0)
			dynamic_cast<GMSHMesh<Dim>*>(mesh)->ExportToGMSH_Elements(testCase->SourceFunction, args.OutputDirectory + "/source", "source");

		// Export exact solution
		if (args.Actions.Export.ExactSolutionToGMSH && testCase->ExactSolution && args.Discretization.Mesher.compare("gmsh") == 0)
			dynamic_cast<GMSHMesh<Dim>*>(mesh)->ExportToGMSH_Elements(testCase->ExactSolution, args.OutputDirectory + "/exsol", "exact solution");

		//----------------------//
		//       Assembly       //
		//----------------------//

		int k = args.Discretization.PolyDegree - 1;
		int reconstructDegree = k + 1;
		int faceDegree = k;
		int cellDegree = k + args.Discretization.RelativeCellPolyDegree;

		FunctionalBasis<Dim>* reconstructionBasis = FunctionalBasisFactory<Dim>::Create(args.Discretization.ElemBasisCode, reconstructDegree, args.Discretization.UsePolynomialSpaceQ);
		FunctionalBasis<Dim>* cellBasis = FunctionalBasisFactory<Dim>::Create(args.Discretization.ElemBasisCode, cellDegree, args.Discretization.UsePolynomialSpaceQ);
		FunctionalBasis<Dim - 1>* faceBasis = FunctionalBasisFactory<Dim-1>::Create(args.Discretization.FaceBasisCode, faceDegree, args.Discretization.UsePolynomialSpaceQ);

		HHOParameters<Dim>* hho = new HHOParameters<Dim>(mesh, args.Discretization.Stabilization, reconstructionBasis, cellBasis, faceBasis, args.Discretization.OrthogonalizeElemBasesCode, args.Discretization.OrthogonalizeFaceBasesCode);

		bool saveMatrixBlocks = false;
		if (args.Solver.SolverCode.compare("uamg") == 0 || args.Solver.PreconditionerCode.compare("uamg") == 0)
			saveMatrixBlocks = true;
		else if (args.Problem.ComputeNormalDerivative && testCase->ExactSolution_Neumann)
			saveMatrixBlocks = true;

		Diffusion_HHO<Dim>* problem = new Diffusion_HHO<Dim>(mesh, testCase, hho, args.Discretization.StaticCondensation, saveMatrixBlocks);

		cout << endl;
		cout << "----------------------------------------------------------" << endl;
		cout << "-                       Assembly                         -" << endl;
		cout << "----------------------------------------------------------" << endl;

		// If full Neumann, check compatibility condition
		if (testCase->BC.Type == PbBoundaryConditions::FullNeumann)
		{
			// Compatibility condition: (f|1) + <neumann|1> = 0
			double integralF = problem->CellSpace.Integral(testCase->SourceFunction);
			double integralN = problem->BoundarySpace.Integral(testCase->BC.NeumannFunction);
			if (abs(integralF + integralN) < Utils::NumericalZero)
				cout << "Compatibility condition: (f|1) + <neumann|1> = " << (integralF + integralN) << endl;
			else
				Utils::Error("Compatibility condition not respected: (f|1) + <neumann|1> = " + to_string(integralF + integralN));
			cout << endl;
		}

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


			if (testCase->BC.Type == PbBoundaryConditions::FullNeumann)
			{
				NumericImageEnforcer imageEnforcer(&problem->SkeletonSpace);
				imageEnforcer.Setup();
				cout << "Kernel coefficient = " << imageEnforcer.CheckKernel(problem->A) << endl;
				cout << "Orthogonality factor of rhs = " << imageEnforcer.OrthogonalityFactor(problem->b) << endl;
				imageEnforcer.ProjectOntoImage(problem->b);
				cout << "Orthogonality factor of rhs = " << imageEnforcer.OrthogonalityFactor(problem->b) << " (after projection onto Im(A)) " << endl;
				cout << endl;
			}

			Vector systemSolution;

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
				iterativeSolver->ComputeExactSolution = Utils::ProgramArgs.Actions.Export.ErrorToGMSH || problem->A.rows() <= 2000;
				if (args.Solver.ComputeIterL2Error && testCase->ExactSolution && args.Discretization.StaticCondensation)
				{
					iterativeSolver->OnNewSolution = [&problem, &testCase](IterationResult& result, const Vector& faceSolution)
					{
						Vector reconstructedSolution = problem->ReconstructHigherOrderApproximationFromFaceCoeffs(faceSolution);
						result.L2Error = problem->L2Error(testCase->ExactSolution, reconstructedSolution);

						//cout << "                                                                                           " << std::setprecision(8) << (0.5*faceSolution.dot(problem->A * faceSolution) - problem->b.dot(faceSolution)) << endl;
					};
				}

				setupTimer.Start();
				if (Utils::ProgramArgs.Solver.SolverCode.compare("uamg") == 0 || Utils::ProgramArgs.Solver.PreconditionerCode.compare("uamg") == 0)
					iterativeSolver->Setup(problem->A, problem->A_T_T, problem->A_T_ndF, problem->A_ndF_ndF);
				else
					solver->Setup(problem->A);
				setupTimer.Stop();

				cout << "Solving..." << endl;
				solvingTimer.Start();
				systemSolution = iterativeSolver->Solve(problem->b, args.Solver.InitialGuessCode);
				solvingTimer.Stop();
				cout << iterativeSolver->IterationCount << " iterations." << endl << endl;

				Multigrid* mg = dynamic_cast<Multigrid*>(iterativeSolver);
				if (mg)
				{
					int sizeTime = 7;
					int sizeWork = 5;

					double totalTime = solvingTimer.CPU().InMilliseconds();
					MFlops totalWork = iterativeSolver->SolvingComputationalWork;
					cout << "\t                             | CPU time |   Work " << endl;
					cout << "\t-------------------------------------------------" << endl;
					cout << "\tSmoothing and res. computing | " << setw(sizeTime) << (int)round((totalTime - mg->IntergridTransferTimer.CPU().InMilliseconds() - mg->CoarseSolverTimer.CPU().InMilliseconds()) / totalTime * 100) << "% | " << setw(sizeWork) << (int)round((totalWork - mg->IntergridTransferCost - mg->CoarseSolverCost) / totalWork * 100) << "%" << endl;
					cout << "\tIntergrid transfers          | " << setw(sizeTime) << (int)round(mg->IntergridTransferTimer.CPU().InMilliseconds() / totalTime * 100) << "% | " << setw(sizeWork) << (int)round(mg->IntergridTransferCost / totalWork * 100) << "%" << endl;
					cout << "\tCoarse solver                | " << setw(sizeTime) << (int)round(mg->CoarseSolverTimer.CPU().InMilliseconds() / totalTime * 100) << "% | " << setw(sizeWork) << (int)round(mg->CoarseSolverCost / totalWork * 100) << "%" << endl;
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
				systemSolution = solver->Solve(problem->b);
				solvingTimer.Stop();
				cout << endl;
			}
			totalTimer.Stop();

			SolverFactory<Dim>::PrintStats(solver, setupTimer, solvingTimer, totalTimer);

			if (args.Actions.Export.ErrorToGMSH && iterativeSolver)
				problem->ExportErrorToGMSH(iterativeSolver->ExactSolution - systemSolution, out);

			delete solver;

			Vector hybridSolution;
			Vector reconstructedSolution;
			if (args.Actions.Export.SolutionVectors || args.Actions.Export.SolutionToGMSH || testCase->ExactSolution || testCase->BC.Type == PbBoundaryConditions::FullNeumann)
			{
				cout << "----------------------------------------------------------" << endl;
				cout << "-                     Post-processing                    -" << endl;
				cout << "----------------------------------------------------------" << endl;

				if (args.Discretization.StaticCondensation)
				{
					cout << "Solving cell unknowns..." << endl;
					hybridSolution = problem->HybridCoeffsBySolvingCellUnknowns(systemSolution);
				}
				else
					hybridSolution = problem->HybridCoeffsByAddingDirichletCoeffs(systemSolution);
				
				cout << "Reconstruction of higher order approximation..." << endl;
				reconstructedSolution = problem->ReconstructHigherOrderApproximationFromHybridCoeffs(hybridSolution);
			}

			//-----------------------------//
			//       Solution export       //
			//-----------------------------//

			if (args.Actions.Export.SolutionVectors)
			{
				if (args.Discretization.StaticCondensation)
					out.ExportVector(systemSolution, "solutionFaces");
				out.ExportVector(hybridSolution, "solutionHybrid");
				out.ExportVector(reconstructedSolution, "solutionHigherOrder");
			}

			if (args.Actions.Export.MeshToMatlab)
			{
				mesh->ExportToMatlab(args.OutputDirectory);
				mesh->ExportToMatlab2(args.OutputDirectory + "/mesh.m");
			}

			if (args.Actions.Export.SolutionToGMSH && args.Discretization.Mesher.compare("gmsh") == 0)
				problem->ExportSolutionToGMSH(reconstructedSolution, out, args.Actions.Export.VisuTolerance, args.Actions.Export.VisuMaxRefinements);

			//----------------------//
			//       L2 error       //
			//----------------------//

			if (args.Actions.ComputeErrors)
			{
				cout << std::scientific << std::setprecision(2);

				// Check mean value if full Neumann conditions
				if (testCase->BC.Type == PbBoundaryConditions::FullNeumann)
				{
					if (testCase->ExactSolution)
						cout << "Mean value of exact solution = " << problem->ReconstructSpace.MeanValue(testCase->ExactSolution) << endl;
					double meanValue = problem->ReconstructSpace.MeanValue(reconstructedSolution);
					cout << "Mean value = " << meanValue << endl;

					if (abs(meanValue) > Utils::NumericalZero)
					{
						ZeroMeanEnforcer integralZeroOnDomain(&problem->ReconstructSpace);
						integralZeroOnDomain.Setup();
						integralZeroOnDomain.Enforce(reconstructedSolution);
						meanValue = problem->ReconstructSpace.MeanValue(reconstructedSolution);
						cout << "Mean value = " << meanValue << " (after correction)" << endl;
					}
				}


				if (testCase->ExactSolution)
				{
					double error = problem->L2Error(testCase->ExactSolution, reconstructedSolution);
					cout << endl << "L2 Error = " << std::scientific << error << endl;
					problem->AssertSchemeConvergence(error);
				}

				if (args.Problem.ComputeNormalDerivative && testCase->ExactSolution_Neumann)
				{
					if (args.Problem.BCCode.compare("d") == 0)
					{
						SparseMatrix Theta_T_bF_transpose = problem->Theta_T_bF_transpose();
						auto nBoundaryElemUnknowns = Theta_T_bF_transpose.cols();

						SparseMatrix S_iF_bF_transpose = Theta_T_bF_transpose * problem->A_T_ndF.topRows(nBoundaryElemUnknowns) + problem->A_ndF_dF.transpose();

						SparseMatrix tmp = problem->A_T_dF.topRows(nBoundaryElemUnknowns).transpose() * Theta_T_bF_transpose.transpose();
						SparseMatrix S_bF_bF = tmp + problem->A_dF_dF;

						Vector g_D = problem->DirichletSpace.Project(testCase->BC.DirichletFunction);
						Vector source = problem->CellSpace.InnerProdWithBasis(testCase->SourceFunction);
						Vector sourceElemBoundary = problem->ExtractElemBoundary(source);

						Vector normalDerivativeRHS = S_iF_bF_transpose * systemSolution + S_bF_bF * g_D - Theta_T_bF_transpose * sourceElemBoundary;
						Vector normalDerivative = problem->BoundarySpace.SolveMassMatrix(normalDerivativeRHS);

						Vector exact = problem->BoundarySpace.Project(testCase->ExactSolution_Neumann);
						double error = problem->BoundarySpace.RelativeL2Error(normalDerivative, exact);
						cout << endl << "L2 Error (normal derivative) (old) = " << std::scientific << error << endl;

						error = problem->L2ErrorNormalDerivative(testCase->ExactSolution_Neumann, normalDerivative);
						cout << endl << "L2 Error (normal derivative) (new) = " << std::scientific << error << endl;

					}
					else
						Utils::Warning("The normal derivative is computed only for Dirichlet problems.");
				}
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