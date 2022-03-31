#pragma once
#include <iomanip>
#include "../ProgramArguments.h"
#include "../Discretizations/FEM/BiHarmonicMixedFormGlowinski_FEM.h"
#include "../TestCases/BiHarmonic/BiHarTestCaseFactory.h"
#include "../Mesher/MeshFactory.h"
#include "../Solver/SolverFactory.h"
#include "../Solver/BiHarmonic/BiHarmonicCG.h"
#include "../Solver/BiHarmonic/BiHarmonicGradientDescent.h"
#include "../Utils/ExportModule.h"

// Biharmonic equation in mixed form with mixed (homogeneous) Dirichlet-Neumann BC

template <int Dim>
class Program_BiHarmonic_FEM
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

		ExportModule out(args.OutputDirectory, "", args.Actions.Export.ValueSeparator);

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

		// Export source
		if (args.Actions.Export.SourceToGMSH && args.Discretization.Mesher.compare("gmsh") == 0)
			dynamic_cast<GMSHMesh<Dim>*>(mesh)->ExportToGMSH_Nodes(testCase->SourceFunction, args.OutputDirectory + "/source", "source");

		// Export exact solution
		if (args.Actions.Export.ExactSolutionToGMSH && testCase->ExactSolution && args.Discretization.Mesher.compare("gmsh") == 0)
			dynamic_cast<GMSHMesh<Dim>*>(mesh)->ExportToGMSH_Nodes(testCase->ExactSolution, args.OutputDirectory + "/exsol", "exact solution");

		cout << endl;

		//----------------------//
		//       Assembly       //
		//----------------------//

		if (args.Problem.Scheme.compare("f") == 0)
			Utils::FatalError("Falk's scheme is not implemented in FEM.");
		else if (args.Problem.Scheme.compare("g") == 0)
		{
			auto Dirichlet = BoundaryConditions::HomogeneousDirichletEverywhere();
			mesh->SetBoundaryConditions(&Dirichlet);
			mesh->FillBoundaryAndInteriorVertexLists();
			mesh->FillDirichletAndNeumannVertexLists();
			cout << "Scheme: Glowinski" << endl;
		}
		else
			Utils::FatalError("Unknown scheme '" + args.Problem.Scheme + "'. Check -sch parameter. Possible values are 'f' and 'g'.");

		FunctionalBasis<Dim>* basis = new FunctionalBasis<Dim>(args.Discretization.ElemBasisCode, args.Discretization.PolyDegree, args.Discretization.UsePolynomialSpaceQ);

		BiHarmonicMixedForm_FEM<Dim>* biHarPb = new BiHarmonicMixedFormGlowinski_FEM<Dim>(mesh, testCase, basis);

		cout << endl;
		cout << "------------------------------------------------------" << endl;
		cout << "-          Assembly of diffusion problem             -" << endl;
		cout << "------------------------------------------------------" << endl;
		Timer assemblyTimer;
		assemblyTimer.Start();

		biHarPb->Setup();

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

			int blockSizeForBlockSolver = args.Solver.BlockSize != -1 ? args.Solver.BlockSize : 1;
			Solver* diffSolver = SolverFactory<Dim>::CreateSolver(args, nullptr, blockSizeForBlockSolver, out);

			Timer setupTimer;
			Timer solvingTimer1;
			Timer solvingTimer2;
			Timer totalTimer;
			totalTimer.Start();

			cout << "Solver: " << *diffSolver << endl << endl;

			setupTimer.Start();
			diffSolver->Setup(biHarPb->DiffPb().A);
			setupTimer.Stop();

			biHarPb->SetDiffSolver(diffSolver);


			cout << "-------------------------------------" << endl;
			cout << "-     Solve biharmonic problem     -" << endl;
			cout << "-------------------------------------" << endl;

			// Solve 1st problem with f as source
			Vector theta_f = biHarPb->FindCompatibleTheta();
			Vector lambda_f = biHarPb->Solve1stDiffProblemWithFSource(theta_f);
			// Solve 2nd problem and extract the boundary (or normal derivative)
			Vector u_boundary_f = biHarPb->Solve2ndDiffProblem(lambda_f, true);

			// A(theta) -> -boundary(u) is s.p.d.
			// Note that because of the minus sign, u_boundary0 = -Af(theta_f)
			// We want to solve A(theta0) = u_boundary_f, in order to have
			// theta = theta_f+theta0 with
			//         A(theta) = Af(theta_f) + A(theta0) = -u_boundary_f + u_boundary_f = 0.
			Vector b = u_boundary_f;


			DenseMatrix A; // computed explicitly only if explicit solver or export requested
			if (args.Solver.BiHarmonicSolverCode.compare("lu") == 0 || args.Solver.BiHarmonicSolverCode.compare("jcg") == 0 || args.Actions.Export.LinearSystem)
			{
				cout << "Computation of the matrix..." << endl;
				A = biHarPb->Matrix();

				if (args.Actions.Export.LinearSystem)
				{
					cout << "Export linear system..." << endl;
					out.ExportMatrix(A, "matrix");
					out.ExportVector(b, "b");
				}
			}


			//-------------------------------------//
			//            Create solver            //
			//-------------------------------------//

			Solver* biHarSolver = nullptr;
			if (Utils::EndsWith(args.Solver.BiHarmonicSolverCode, "cg"))
				biHarSolver = new BiHarmonicCG(biHarPb, args.Solver.Restart);
			else if (args.Solver.BiHarmonicSolverCode.compare("gd") == 0)
				biHarSolver = new BiHarmonicGradientDescent(biHarPb);
			else if (args.Solver.BiHarmonicSolverCode.compare("lu") == 0)
				biHarSolver = new EigenLU();
			else
				Utils::FatalError("Unknown biharmonic solver '" + args.Solver.BiHarmonicSolverCode + "'");

			cout << "Solver: " << *biHarSolver << endl << endl;

			IterativeSolver* biHarIterSolver = dynamic_cast<IterativeSolver*>(biHarSolver);
			if (biHarIterSolver)
			{
				biHarIterSolver->StoppingCrit = args.Solver.StoppingCrit;
				biHarIterSolver->Tolerance = args.Solver.Tolerance;
				biHarIterSolver->StagnationConvRate = args.Solver.StagnationConvRate;
				biHarIterSolver->MaxIterations = args.Solver.MaxIterations;
				// Compute L2-error at each iteration
				if ((args.Solver.ComputeIterL2Error || args.Actions.Export.IterationL2Errors) && testCase->ExactSolution)
				{
					biHarIterSolver->OnNewSolution = [&biHarPb, &testCase, &theta_f, &u_boundary_f, &b](IterationResult& result, const Vector& theta_0)
					{
						Vector reconstructedSolution = biHarPb->ComputeSolution(theta_f + theta_0);
						result.L2Error = biHarPb->DiffPb().L2Error(testCase->ExactSolution, reconstructedSolution);

						Vector lambda2 = biHarPb->Solve1stDiffProblemWithFSource(theta_f + theta_0);
						Vector u_boundary = biHarPb->Solve2ndDiffProblem(lambda2, true);
						/* // same thing
						Vector lambda0 = biHarPb->Solve1stDiffProblemWithZeroSource(theta);
						Vector u_boundary_0 = biHarPb->Solve2ndDiffProblem(lambda0, true);
						Vector u_boundary_b = u_boundary_f + u_boundary_0;*/
						//cout << "                                                                ||u_boundary|| = " << std::scientific << sqrt(biHarPb->L2InnerProdOnBoundary(u_boundary_b, u_boundary_b)) << endl;
						result.BoundaryL2Norm = sqrt(biHarPb->L2InnerProdOnBoundary(u_boundary, u_boundary));

						// Energy functional
						//Vector lambda = biHarPb->Solve1stDiffProblemWithZeroSource(theta);
						//Vector Atheta = -biHarPb->Solve2ndDiffProblem(lambda, true);
						//cout << "                                                                                " << std::setprecision(8) << 0.5 * theta.dot(Atheta) - b.dot(theta) << endl;
					};
				}
				// Export iteration results
				if (args.Actions.Export.Iterations || args.Actions.Export.IterationResiduals || args.Actions.Export.IterationL2Errors)
				{
					if (args.Actions.Export.Iterations)
						out.CleanFile("iterations.csv");
					if (args.Actions.Export.IterationResiduals)
						out.CleanFile("iteration_residuals.dat");
					if (args.Actions.Export.IterationL2Errors)
						out.CleanFile("iteration_l2errors.dat");

					biHarIterSolver->OnIterationEnd = [&args, &out](const IterationResult& result)
					{
						if (args.Actions.Export.Iterations)
							out.Export(result, "iterations");
						if (args.Actions.Export.IterationResiduals)
							out.ExportNewVectorValue(result.NormalizedResidualNorm, "iteration_residuals");
						if (args.Actions.Export.IterationL2Errors && result.L2Error != -1)
							out.ExportNewVectorValue(result.L2Error, "iteration_l2errors");
					};
				}

				if (Utils::EndsWith(args.Solver.BiHarmonicSolverCode, "cg"))
				{
					BiHarmonicCG* cg = static_cast<BiHarmonicCG*>(biHarIterSolver);
					if (args.Solver.BiHarmonicSolverCode.compare("jcg") == 0)
					{
						DenseBlockJacobiPreconditioner* jp = new DenseBlockJacobiPreconditioner(1);
						jp->Setup(A);
						cg->Precond = jp;
					}
					else if (args.Solver.BiHarmonicSolverCode.compare("cg") == 0)
					{
						IdentityPreconditioner* jp = new IdentityPreconditioner();
						cg->Precond = jp;
					}
					else
						Utils::FatalError("Unknown preconditioner. Check -bihar-solver " + args.Solver.BiHarmonicSolverCode + ".");
				}

				// Give the Laplacian matrix so that the Work Units are computed with respect to it
				biHarIterSolver->Setup(biHarPb->DiffPb().A);
			}
			else // direct solver
			{
				cout << "Factorization..." << endl;
				biHarSolver->Setup(A);
			}


			//--------------------------------------------------------------------------//
			/*ProblemArguments pbDiff;
			pbDiff.TestCaseCode = "fullneumann";
			pbDiff.SourceCode = "exp";
			pbDiff.BCCode = "n2";
			DiffusionTestCase<Dim>* diffTestCase = DiffTestCaseFactory<Dim>::Create(pbDiff);

			Diffusion_HHO<Dim>* diffPb = new Diffusion_HHO<Dim>(mesh, diffTestCase, hho, true, false);
			diffPb->InitHHO();

			HigherOrderBoundary<Dim> higherOrderBoundary(diffPb);
			higherOrderBoundary.Setup(false, true);

			Vector projectedSolution = diffPb->ReconstructSpace.Project(diffTestCase->ExactSolution);
			Vector projectedSolutionBoundary = higherOrderBoundary.ExtractBoundaryElements(projectedSolution);
			Vector normalDeriv = higherOrderBoundary.NormalDerivative(projectedSolutionBoundary);

			Vector projectedNeumann = higherOrderBoundary.BoundarySpace.Project(diffTestCase->BC.NeumannFunction);

			double error = (normalDeriv - projectedNeumann).norm() / normalDeriv.norm();




			DomFunction source = [](const DomPoint& p)
			{
				return p.X * p.Y;
			};
			Vector projectedSource = diffPb->ReconstructSpace.Project(source);

			Vector b_sourceFromContinuous = diffPb->AssembleSourceTerm(source);
			Vector b_sourceFromProjection = diffPb->AssembleSourceTerm(projectedSource);

			error = (b_sourceFromContinuous - b_sourceFromProjection).norm() / b_sourceFromContinuous.norm();
			cout << error << endl;*/

			if (args.Actions.UnitTests)
			{
				auto n = A.rows();
				A = (A + A.transpose()) / 2.0;

				Eigen::EigenSolver<DenseMatrix> es(A);
				double det = A.determinant();
				//auto v = A.eigenvalues();
				//cout << v << endl;
				Eigen::VectorXcd eigenvalues = es.eigenvalues();
				cout << eigenvalues << endl;
				Eigen::MatrixXcd eigenvectors = es.eigenvectors();
				Eigen::VectorXcd kernelVector = eigenvectors.col(n - 1);
				cout << "---------------------" << endl << kernelVector << endl;
				cout << "---------------------" << endl << eigenvectors.col(n - 2) << endl;
				cout << "---------------------" << endl << eigenvectors.col(n - 3) << endl;

				Vector lambda = biHarPb->Solve1stDiffProblemWithZeroSource(kernelVector.real());
				//cout << lambda.norm() << endl;

				//DiffPb().ExportReconstructedVectorToGMSH(lambda, out, "lambda");
				Vector solPb2 = biHarPb->Solve2ndDiffProblem(lambda, false);
				//biHarPb->DiffPb().ExportReconstructedVectorToGMSH(solPb2, out, "solPb2");


				Vector zero = -biHarPb->Solve2ndDiffProblem(lambda, true);
				cout << zero.norm() << endl;

				//out.ExportMatrix(A, "matrix");
			}
			//--------------------------------------------------------------------------//

			if (!biHarIterSolver)
			{
				Utils::Empty(A);
				cout << "Solve linear system..." << endl;
			}

			Vector theta_0 = biHarSolver->Solve(b);

			delete biHarSolver;

			cout << "Compute solution..." << endl;
			Vector solution = biHarPb->ComputeSolution(theta_f + theta_0);

			//-----------------------------//
			//       Solution export       //
			//-----------------------------//

			if (args.Actions.Export.SolutionVectors)
				out.ExportVector(solution, "solution");

			if (args.Actions.Export.MeshToMatlab)
			{
				mesh->ExportToMatlab(args.OutputDirectory);
				mesh->ExportToMatlab2(args.OutputDirectory + "/mesh.m");
			}

			if (args.Actions.Export.SolutionToGMSH && args.Discretization.Mesher.compare("gmsh") == 0)
				mesh->ExportToGMSH_Nodes(solution, out.GetFilePathPrefix(), "solution");

			//----------------------//
			//       L2 error       //
			//----------------------//

			if (testCase->ExactSolution)
			{
				double error = biHarPb->DiffPb().L2Error(testCase->ExactSolution, solution);
				cout << endl << "L2 Error = " << std::scientific << error << endl;
			}
		}

		//--------------------------//
		//       Deallocation       //
		//--------------------------//

		delete mesh;
		delete basis;
		delete biHarPb;
		delete testCase;
		GaussLegendre::Free();
	}
};