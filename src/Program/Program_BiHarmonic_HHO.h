#pragma once
#include <iomanip>
#include "../ProgramArguments.h"
#include "../Discretizations/HHO/BiHarmonicMixedFormFalk_HHO.h"
#include "../Discretizations/HHO/BiHarmonicMixedFormGlowinski_HHO.h"
#include "../TestCases/BiHarmonic/BiHarTestCaseFactory.h"
#include "../Mesher/MeshFactory.h"
#include "../Solver/SolverFactory.h"
#include "../Solver/BiHarmonic/BiHarmonicCG.h"
#include "../Solver/BiHarmonic/BiHarmonicGradientDescent.h"
#include "../Solver/BiHarmonic/BiharPatchPreconditioner.h"
#include "../Utils/ExportModule.h"

// Biharmonic equation in mixed form with mixed (homogeneous) Dirichlet-Neumann BC

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

		ExportModule out(args.OutputDirectory, "", args.Actions.Export.ValueSeparator);

		//----------//
		//   Mesh   //
		//----------//

		cout << "-----------------------------------------------------------" << endl;
		cout << "-                   Mesh construction                     -" << endl;
		cout << "-----------------------------------------------------------" << endl;

		Mesh<Dim>* mesh = MeshFactory<Dim>::BuildMesh(args, testCase);
		if (args.Discretization.Mesher.compare("gmsh") == 0)
			GMSHMesh<Dim>::CloseGMSH();

		mesh->RenumberBoundaryElementsFirst();

		cout << "Mesh storage > " << Utils::MemoryString(mesh->MemoryUsage()) << endl;

		// Export source
		if (args.Actions.Export.SourceToGMSH && args.Discretization.Mesher.compare("gmsh") == 0)
			dynamic_cast<GMSHMesh<Dim>*>(mesh)->ExportToGMSH_Elements(testCase->SourceFunction, args.OutputDirectory + "/source", "source");

		// Export exact solution
		if (args.Actions.Export.ExactSolutionToGMSH && testCase->ExactSolution && args.Discretization.Mesher.compare("gmsh") == 0)
			dynamic_cast<GMSHMesh<Dim>*>(mesh)->ExportToGMSH_Elements(testCase->ExactSolution, args.OutputDirectory + "/exsol", "exact solution");

		cout << endl;

		//----------------------//
		//       Assembly       //
		//----------------------//

		if (args.Problem.Scheme.compare("f") == 0)
		{
			auto FullNeumann = BoundaryConditions::HomogeneousNeumannEverywhere();
			mesh->SetBoundaryConditions(&FullNeumann);
			cout << "Scheme: Falk" << endl;
		}
		else if (args.Problem.Scheme.compare("g") == 0)
		{
			auto Dirichlet = BoundaryConditions::HomogeneousDirichletEverywhere();
			mesh->SetBoundaryConditions(&Dirichlet);
			cout << "Scheme: Glowinski" << endl;
		}
		else
			Utils::FatalError("Unknown scheme '" + args.Problem.Scheme + "'. Check -sch parameter. Possible values are 'f' and 'g'.");

		int k = args.Discretization.PolyDegree - 1;
		int reconstructDegree = k + 1;
		int faceDegree = k;
		int cellDegree = k + args.Discretization.RelativeCellPolyDegree;

		FunctionalBasis<Dim>* reconstructionBasis = new FunctionalBasis<Dim>(args.Discretization.ElemBasisCode, reconstructDegree, args.Discretization.UsePolynomialSpaceQ);
		FunctionalBasis<Dim>* cellBasis = new FunctionalBasis<Dim>(args.Discretization.ElemBasisCode, cellDegree, args.Discretization.UsePolynomialSpaceQ);
		FunctionalBasis<Dim - 1>* faceBasis = new FunctionalBasis<Dim - 1>(args.Discretization.FaceBasisCode, faceDegree, args.Discretization.UsePolynomialSpaceQ);

		HHOParameters<Dim>* hho = new HHOParameters<Dim>(mesh, args.Discretization.Stabilization, reconstructionBasis, cellBasis, faceBasis, args.Discretization.OrthogonalizeElemBasesCode, args.Discretization.OrthogonalizeFaceBasesCode);

		bool saveMatrixBlocks = args.Solver.SolverCode.compare("uamg") == 0 || args.Solver.PreconditionerCode.compare("uamg") == 0;
		BiHarmonicMixedForm_HHO<Dim>* biHarPb;
		if (args.Problem.Scheme.compare("f") == 0)
			biHarPb = new BiHarmonicMixedFormFalk_HHO<Dim>(mesh, testCase, hho, args.Solver.BiHarReconstructBoundary, saveMatrixBlocks);
		else if (args.Problem.Scheme.compare("g") == 0)
			biHarPb = new BiHarmonicMixedFormGlowinski_HHO<Dim>(mesh, testCase, hho, saveMatrixBlocks);
		else
			Utils::FatalError("Unknown scheme '" + args.Problem.Scheme + "'. Check -sch parameter. Possible values are 'f' and 'g'.");

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

			int blockSizeForBlockSolver = args.Solver.BlockSize != -1 ? args.Solver.BlockSize : faceBasis->Size();
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
				if (Utils::ProgramArgs.Solver.SolverCode.compare("uamg") == 0 || Utils::ProgramArgs.Solver.PreconditionerCode.compare("uamg") == 0)
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


			cout << "-------------------------------------" << endl;
			cout << "-     Solve biharmonic problem      -" << endl;
			cout << "-------------------------------------" << endl;


			Vector discreteExactSolution;
			double exactSolutionL2Norm = 0;
			if (testCase->ExactSolution)
			{
				discreteExactSolution = biHarPb->DiffPb().ReconstructSpace.Project(testCase->ExactSolution);
				exactSolutionL2Norm = biHarPb->DiffPb().ReconstructSpace.L2Norm(discreteExactSolution);
			}

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

			//biHarPb->DiffPb().ExportReconstructedVectorToGMSH(biHarPb->Solve2ndDiffProblem(lambda_f, false), out, "u_f", args.Actions.Export.VisuTolerance, args.Actions.Export.VisuMaxRefinements);
			//biHarPb->DiffPb().ExportReconstructedVectorToGMSH(biHarPb->Solve1stDiffProblemWithZeroSource(u_boundary_f), out, "lambda_f", args.Actions.Export.VisuTolerance, args.Actions.Export.VisuMaxRefinements);


			DenseMatrix A; // computed explicitly only if explicit solver or export requested
			if (args.Solver.BiHarmonicSolverCode.compare("lu") == 0  || 
				(args.Solver.BiHarmonicSolverCode.compare("cg") == 0 && (args.Solver.BiHarmonicPreconditionerCode.compare("j") == 0 || args.Solver.BiHarmonicPreconditionerCode.compare("bj") == 0)) ||
				args.Actions.Export.LinearSystem)
			{
				cout << "Computation of the matrix..." << endl;
				A = biHarPb->Matrix();
				//cout << "Matrix: " << endl << A << endl << endl;

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
			if (args.Solver.BiHarmonicSolverCode.compare("cg") == 0)
				biHarSolver = new BiHarmonicCG(biHarPb, args.Solver.Restart);
			else if (args.Solver.BiHarmonicSolverCode.compare("gd") == 0)
				biHarSolver = new BiHarmonicGradientDescent(biHarPb);
			else if (args.Solver.BiHarmonicSolverCode.compare("lu") == 0)
				biHarSolver = new EigenLU();
			else
				Utils::FatalError("Unknown biharmonic solver '" + args.Solver.BiHarmonicSolverCode + "'");

			cout << "Solver: " << *biHarSolver << endl;

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
					biHarIterSolver->OnNewSolution = [&biHarPb, &testCase, &theta_f, &u_boundary_f, &b, &discreteExactSolution , &exactSolutionL2Norm](IterationResult& result, const Vector& theta_0)
					{
						Vector reconstructedLap, reconstructedSolution;
						std::tie(reconstructedLap, reconstructedSolution) = biHarPb->ComputeSolution(theta_f + theta_0);
						result.L2Error = biHarPb->DiffPb().ReconstructSpace.L2Norm(reconstructedSolution - discreteExactSolution) / exactSolutionL2Norm;

						//Vector lambda2 = biHarPb->Solve1stDiffProblemWithFSource(theta_f + theta_0);
						//Vector u_boundary = biHarPb->Solve2ndDiffProblem(lambda2, true);
						/* // same thing
						Vector lambda0 = biHarPb->Solve1stDiffProblemWithZeroSource(theta);
						Vector u_boundary_0 = biHarPb->Solve2ndDiffProblem(lambda0, true);
						Vector u_boundary_b = u_boundary_f + u_boundary_0;*/
						//cout << "                                                                ||u_boundary|| = " << std::scientific << sqrt(biHarPb->L2InnerProdOnBoundary(u_boundary_b, u_boundary_b)) << endl;
						//result.BoundaryL2Norm = sqrt(biHarPb->L2InnerProdOnBoundary(u_boundary, u_boundary));

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

				if (args.Solver.BiHarmonicSolverCode.compare("cg") == 0)
				{
					BiHarmonicCG* cg = static_cast<BiHarmonicCG*>(biHarIterSolver);
					if (args.Problem.Scheme.compare("g") == 0)
					{
						BiHarmonicMixedFormGlowinski_HHO<Dim>* gloScheme = static_cast<BiHarmonicMixedFormGlowinski_HHO<Dim>*>(biHarPb);
						if (args.Solver.BiHarmonicPreconditionerCode.compare("j") == 0)
						{
							cout << "Preconditioner: Jacobi" << endl << endl;
							DenseBlockJacobiPreconditioner* p = new DenseBlockJacobiPreconditioner(1);
							p->Setup(A);
							cg->Precond = p;
						}
						else if (args.Solver.BiHarmonicPreconditionerCode.compare("bj") == 0)
						{
							cout << "Preconditioner: block Jacobi" << endl << endl;
							DenseBlockJacobiPreconditioner* p = new DenseBlockJacobiPreconditioner(hho->nFaceUnknowns);
							p->Setup(A);
							cg->Precond = p;
						}
						else if (args.Solver.BiHarmonicPreconditionerCode.compare("p") == 0)
						{
							cout << "Preconditioner: patch" << endl << endl;
							BiharPatchPreconditioner<Dim>* p = new BiharPatchPreconditioner<Dim>(*gloScheme, args.Solver.NeighbourhoodDepth);
							p->Setup();
							cg->Precond = p;
						}
						else if (args.Solver.BiHarmonicPreconditionerCode.compare("dp") == 0)
						{
							cout << "Preconditioner: diagonal patch" << endl << endl;
							BiharPatchPreconditioner<Dim>* p = new BiharPatchPreconditioner<Dim>(*gloScheme, args.Solver.NeighbourhoodDepth, true);
							p->Setup();
							cg->Precond = p;
						}
						else if (args.Solver.BiHarmonicPreconditionerCode.compare("no") == 0)
						{
							cout << "Preconditioner: none" << endl << endl;
							IdentityPreconditioner* p = new IdentityPreconditioner();
							cg->Precond = p;
						}
						else
							Utils::FatalError("Unknown preconditioner. Check -bihar-prec " + args.Solver.BiHarmonicPreconditionerCode + ".");
					}
					else if (args.Problem.Scheme.compare("f") == 0)
					{
						cout << "Preconditioner: none" << endl << endl;
						IdentityPreconditioner* p = new IdentityPreconditioner();
						cg->Precond = p;
					}
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
				biHarPb->DiffPb().ExportReconstructedVectorToGMSH(solPb2, out, "solPb2");


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

			theta_f += theta_0;

			Vector reconstructedLap, reconstructedSolution;
			std::tie(reconstructedLap, reconstructedSolution) = biHarPb->ComputeSolution(theta_f);

			//-----------------------------//
			//       Solution export       //
			//-----------------------------//

			if (args.Actions.Export.SolutionVectors)
				out.ExportVector(reconstructedSolution, "solutionHigherOrder");

			if (args.Actions.Export.MeshToMatlab)
			{
				mesh->ExportToMatlab(args.OutputDirectory);
				mesh->ExportToMatlab2(args.OutputDirectory + "/mesh.m");
			}

			if (args.Actions.Export.SolutionToGMSH && args.Discretization.Mesher.compare("gmsh") == 0)
				biHarPb->DiffPb().ExportReconstructedVectorToGMSH(reconstructedSolution, out, "solution", args.Actions.Export.VisuTolerance, args.Actions.Export.VisuMaxRefinements);

			//----------------------//
			//       L2 error       //
			//----------------------//

			cout.precision(2);

			if (testCase->ExactSolution)
			{
				Vector error = discreteExactSolution - reconstructedSolution;

				if (args.Actions.Export.ErrorToGMSH && args.Discretization.Mesher.compare("gmsh") == 0)
					biHarPb->DiffPb().ExportReconstructedVectorToGMSH(error, out, "error", args.Actions.Export.VisuTolerance, args.Actions.Export.VisuMaxRefinements);

				if (args.Actions.Export.AbsErrorToGMSH && args.Discretization.Mesher.compare("gmsh") == 0)
					biHarPb->DiffPb().ExportReconstructedVectorToGMSH(error, out, "abs_error", args.Actions.Export.VisuTolerance, args.Actions.Export.VisuMaxRefinements, true);

				double solutionError = biHarPb->DiffPb().ReconstructSpace.L2Norm(error) / exactSolutionL2Norm;
				cout << endl << "L2 Error (solution) = " << std::scientific << solutionError << endl;
			}

			if (testCase->MinusLaplacianOfSolution)
			{
				Vector discreteExact = biHarPb->DiffPb().ReconstructSpace.Project(testCase->MinusLaplacianOfSolution);
				double error = biHarPb->DiffPb().ReconstructSpace.RelativeL2Error(reconstructedLap, discreteExact);
				cout << "L2 Error (laplacian) = " << std::scientific << error << endl;
			}

			if (testCase->MinusLaplacianOfSolution_Dirichlet)
			{
				Vector discreteExact = biHarPb->DiffPb().BoundarySpace.Project(testCase->MinusLaplacianOfSolution_Dirichlet);
				double error = biHarPb->DiffPb().BoundarySpace.RelativeL2Error(theta_f, discreteExact);
				cout << "L2 Error (theta) = " << std::scientific << error << endl;
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