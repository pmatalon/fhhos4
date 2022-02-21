#pragma once
#include <iomanip>
#include "../ProgramArguments.h"
#include "../HHO/BiHarmonicMixedFormFalk_HHO.h"
#include "../HHO/BiHarmonicMixedFormGlowinski_HHO.h"
#include "../TestCases/BiHarmonic/BiHarTestCaseFactory.h"
#include "../Mesher/MeshFactory.h"
#include "../Solver/SolverFactory.h"
#include "../Solver/Krylov/BiHarmonicCG.h"
#include "../Solver/BiHarmonic/BiHarmonicGradientDescent.h"
#include "../Utils/ExportModule.h"

// Bi-harmonic equation in mixed form with mixed (homogeneous) Dirichlet-Neumann BC

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

		cout << "Mesh storage > " << Utils::MemoryString(mesh->MemoryUsage()) << endl;

		// Export source
		if (args.Actions.Export.SourceToGMSH && args.Discretization.Mesher.compare("gmsh") == 0)
			dynamic_cast<GMSHMesh<Dim>*>(mesh)->ExportToGMSH(testCase->SourceFunction, args.OutputDirectory + "/source", "source");

		// Export exact solution
		if (args.Actions.Export.ExactSolutionToGMSH && testCase->ExactSolution && args.Discretization.Mesher.compare("gmsh") == 0)
			dynamic_cast<GMSHMesh<Dim>*>(mesh)->ExportToGMSH(testCase->ExactSolution, args.OutputDirectory + "/exsol", "exact solution");

		//----------------------//
		//       Assembly       //
		//----------------------//

		if (args.Problem.Scheme.compare("f") == 0)
		{
			auto FullNeumann = BoundaryConditions::HomogeneousNeumannEverywhere();
			mesh->SetBoundaryConditions(&FullNeumann);
		}
		else if (args.Problem.Scheme.compare("g") == 0)
		{
			auto Dirichlet = BoundaryConditions::HomogeneousDirichletEverywhere();
			mesh->SetBoundaryConditions(&Dirichlet);
		}
		else
			Utils::FatalError("Unknown scheme '" + args.Problem.Scheme + "'. Check -sch parameter. Possible values are 'f' and 'g'.");

		FunctionalBasis<Dim>* reconstructionBasis = new FunctionalBasis<Dim>(args.Discretization.ElemBasisCode, args.Discretization.PolyDegree, args.Discretization.UsePolynomialSpaceQ);
		FunctionalBasis<Dim>* cellBasis = new FunctionalBasis<Dim>(args.Discretization.ElemBasisCode, args.Discretization.PolyDegree - 1, args.Discretization.UsePolynomialSpaceQ);
		FunctionalBasis<Dim - 1>* faceBasis = new FunctionalBasis<Dim - 1>(args.Discretization.FaceBasisCode, args.Discretization.PolyDegree - 1, args.Discretization.UsePolynomialSpaceQ);

		HHOParameters<Dim>* hho = new HHOParameters<Dim>(mesh, args.Discretization.Stabilization, reconstructionBasis, cellBasis, faceBasis, args.Discretization.OrthogonalizeElemBasesCode, args.Discretization.OrthogonalizeFaceBasesCode);

		bool saveMatrixBlocks = args.Solver.SolverCode.compare("uamg") == 0 || args.Solver.SolverCode.compare("fcguamg") == 0;
		BiHarmonicMixedForm_HHO<Dim>* biHarPb;
		if (args.Problem.Scheme.compare("f") == 0)
			biHarPb = new BiHarmonicMixedFormFalk_HHO<Dim>(mesh, testCase, hho, args.Solver.BiHarReconstructBoundary, /*args.Actions.EnforceDirichletBC,*/ saveMatrixBlocks);
		else if (args.Problem.Scheme.compare("g") == 0)
			biHarPb = new BiHarmonicMixedFormGlowinski_HHO<Dim>(mesh, testCase, hho, args.Solver.BiHarReconstructBoundary, saveMatrixBlocks);
		else
			Utils::FatalError("Unknown scheme '" + args.Problem.Scheme + "'. Check -sch parameter. Possible values are 'f' and 'g'.");

		cout << endl;
		cout << "------------------------------------------------------" << endl;
		cout << "-          Assembly of diffusion problem             -" << endl;
		cout << "------------------------------------------------------" << endl;
		Timer assemblyTimer;
		assemblyTimer.Start();

		biHarPb->Setup();

		/*if (args.Problem.Scheme.compare("f") == 0 && args.Actions.EnforceDirichletBC)
		{
			Mesh<Dim>* meshLast = MeshFactory<Dim>::BuildMesh(args, testCase);
			if (args.Discretization.Mesher.compare("gmsh") == 0)
				GMSHMesh<Dim>::CloseGMSH();

			BiHarmonicMixedFormFalk_HHO<Dim>* falkScheme = static_cast<BiHarmonicMixedFormFalk_HHO<Dim>*>(biHarPb);
			falkScheme->SetupLastPb(meshLast);
		}*/

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


			/*if (args.Problem.Scheme.compare("f") == 0 && args.Actions.EnforceDirichletBC)
			{
				cout << "Enforce Dirichlet BC to the solution, so setup of last solver..." << endl << endl;

				BiHarmonicMixedFormFalk_HHO<Dim>* falkScheme = static_cast<BiHarmonicMixedFormFalk_HHO<Dim>*>(biHarPb);

				Solver* lastSolver = SolverFactory<Dim>::CreateSolver(args, falkScheme->LastPb(), blockSizeForBlockSolver, out);
				lastSolver->Setup(falkScheme->LastPb()->A);

				falkScheme->SetLastPbSolver(lastSolver);
			}*/


			cout << "-------------------------------------" << endl;
			cout << "-     Solve bi-harmonic problem     -" << endl;
			cout << "-------------------------------------" << endl;

			// Solve 1st problem with f as source
			Vector theta0 = biHarPb->FindCompatibleTheta();
			Vector delta = biHarPb->Solve1stDiffProblem(theta0);
			// Solve 2nd problem and extract the boundary (or normal derivative)
			Vector u_boundary0 = biHarPb->Solve2ndDiffProblem(delta, true);

			// A(theta) -> -boundary(u) is s.p.d.
			// Note that because of the minus sign, u_boundary0 = -A0(theta0)
			// We want to solve A(theta1) = u_boundary0, in order to have
			// thetaFinal = theta0+theta1 with
			//         A(thetaFinal) = A0(theta0) + A(theta1) = -u_boundary0 + u_boundary0 = 0.
			Vector b = u_boundary0;


			//-------------------------------------//
			//            Create solver            //
			//-------------------------------------//

			Solver* biHarSolver = nullptr;
			if (args.Solver.BiHarmonicSolverCode.compare("cg") == 0)
				biHarSolver = new BiHarmonicCG<Dim>(biHarPb);
			else if (args.Solver.BiHarmonicSolverCode.compare("gd") == 0)
				biHarSolver = new BiHarmonicGradientDescent<Dim>(biHarPb);
			else if (args.Solver.BiHarmonicSolverCode.compare("lu") == 0)
				biHarSolver = new EigenLU();
			else
				Utils::FatalError("Unknown bi-harmonic solver '" + args.Solver.BiHarmonicSolverCode + "'");

			cout << "Solver: " << *biHarSolver << endl << endl;

			DenseMatrix A; // computed explicitly only if direct solver

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
					biHarIterSolver->OnNewSolution = [&biHarPb, &testCase, &theta0, /*&u_f,*/ &b](IterationResult& result, const Vector& theta)
					{
						Vector reconstructedSolution = biHarPb->ComputeSolution(theta0 + theta);
						//Vector reconstructedSolution = biHarPb->ComputeSolution(theta0+theta) + u_f;
						/*Vector delta = biHarPb->Solve1stDiffProblemWithZeroSource(theta);
						Vector u_0 = biHarPb->Solve2ndDiffProblem(delta, false);
						Vector reconstructedSolution = u_f + u_0;*/

						result.L2Error = biHarPb->DiffPb().L2Error(testCase->ExactSolution, reconstructedSolution);

						Vector lambda = biHarPb->Solve1stDiffProblemWithZeroSource(theta);
						Vector Atheta = -biHarPb->Solve2ndDiffProblem(lambda, true);
						cout << "                                                                                " << std::setprecision(8) << 0.5 * theta.dot(Atheta) - b.dot(theta) << endl;
					};
				}
				// Export iteration results
				if (args.Actions.Export.IterationResiduals || args.Actions.Export.IterationL2Errors)
				{
					if (args.Actions.Export.IterationResiduals)
						out.CleanFile("iteration_residuals");
					if (args.Actions.Export.IterationL2Errors)
						out.CleanFile("iteration_l2errors");

					biHarIterSolver->OnIterationEnd = [&args, &out](const IterationResult& result)
					{
						if (args.Actions.Export.IterationResiduals)
							out.ExportNewVectorValue(result.NormalizedResidualNorm, "iteration_residuals");
						if (args.Actions.Export.IterationL2Errors && result.L2Error != -1)
							out.ExportNewVectorValue(result.L2Error, "iteration_l2errors");
					};
				}

				// Give the Laplacian matrix so that the Work Units are computed with respect to it
				biHarIterSolver->Setup(biHarPb->DiffPb().A);
			}
			else // direct solver
			{
				cout << "Computation of the matrix..." << endl;
				A = biHarPb->Matrix();
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

				out.ExportMatrix(A, "matrix");
			}
			//--------------------------------------------------------------------------//

			if (!biHarIterSolver)
			{
				Utils::Empty(A);
				cout << "Solve linear system..." << endl;
			}

			Vector theta = biHarSolver->Solve(b);

			delete biHarSolver;

			cout << "Compute solution..." << endl;
			Vector reconstructedSolution = biHarPb->ComputeSolution(theta0 + theta);

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