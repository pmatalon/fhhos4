#pragma once
#include <iomanip>
#include "../ProgramArguments.h"
#include "../Discretizations/HHO/BiHarmonicMixedFormFalk_HHO.h"
#include "../Discretizations/HHO/BiHarmonicMixedFormGlowinski_HHO.h"
#include "../TestCases/BiHarmonic/BiHarTestCaseFactory.h"
#include "../Mesher/MeshFactory.h"
#include "../FunctionalBasis/FunctionalBasisFactory.h"
#include "../Solver/SolverFactory.h"
#include "../Solver/BiHarmonic/BiHarmonicCG.h"
#include "../Solver/BiHarmonic/BiHarmonicFCG.h"
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

		cout << "Bihar. stabilization: " << args.Discretization.BiharStabilization << endl;


		mesh->RenumberFacesInOrder_IntNeumDirich();

		// Export mesh
		if (args.Actions.Export.MeshToMatlab)
		{
			mesh->ExportToMatlab(args.OutputDirectory);
			mesh->ExportToMatlab2(args.OutputDirectory + "/mesh.m");
		}

		int k = args.Discretization.PolyDegree - 1;
		int reconstructDegree = k + 1;
		int faceDegree = k;
		int cellDegree = k + args.Discretization.RelativeCellPolyDegree;

		FunctionalBasis<Dim>* reconstructionBasis = FunctionalBasisFactory<Dim>::Create(args.Discretization.ElemBasisCode, reconstructDegree, args.Discretization.UsePolynomialSpaceQ);
		FunctionalBasis<Dim>* cellBasis = FunctionalBasisFactory<Dim>::Create(args.Discretization.ElemBasisCode, cellDegree, args.Discretization.UsePolynomialSpaceQ);
		FunctionalBasis<Dim - 1>* faceBasis = FunctionalBasisFactory<Dim-1>::Create(args.Discretization.FaceBasisCode, faceDegree, args.Discretization.UsePolynomialSpaceQ);

		// double maxNBoundaryFacesAsDouble = (double)cellBasis->Size() / (double)faceBasis->Size();
		// int maxNBoundaryFaces = cellBasis->Size() / faceBasis->Size();
		// cout << "Number of boundary faces per element: " << endl;
		// cout << "\ttheoretical maximum = " << maxNBoundaryFaces << " (" << maxNBoundaryFacesAsDouble << ")" << endl;
		// int largestNBoundaryFaces;
		// Element<Dim>* largestBdryElem = mesh->ElementWithTheMostBoundaryFaces(largestNBoundaryFaces);
		// cout << "\tfound in the mesh   = " << largestNBoundaryFaces << endl;
		// if (cellBasis->Size() < largestNBoundaryFaces * faceBasis->Size())
		// {
		// 	Utils::Warning("The largest number of boundary faces in one element is greater than the maximum value authorized to ensure a well-posed problem.");
		// 	if (args.Discretization.RelativeCellPolyDegree == 0)
		// 		Utils::Warning("A solution might be to increase the polynomial degree in the cells to k+1 (argument -kc 1).");
		// }

		// for (Element<Dim>* e : mesh->BoundaryElements)
		// {
		// 	if (e->HasCoplanarBoundaryFaces())
		// 		Utils::Error("Coplanar boundary faces found in element " + to_string(e->Number) + ", which prevents the problem to be well-posed.");
		// }

		HHOParameters<Dim>* hho = new HHOParameters<Dim>(mesh, args.Discretization.Stabilization, reconstructionBasis, cellBasis, faceBasis, args.Discretization.OrthogonalizeElemBasesCode, args.Discretization.OrthogonalizeFaceBasesCode);

		bool saveMatrixBlocks = args.Solver.SolverCode.compare("uamg") == 0 || args.Solver.PreconditionerCode.compare("uamg") == 0;
		BiHarmonicMixedForm_HHO<Dim>* biHarPb;
		if (args.Problem.Scheme.compare("f") == 0)
			biHarPb = new BiHarmonicMixedFormFalk_HHO<Dim>(mesh, testCase, hho, args.Solver.BiHarReconstructBoundary, saveMatrixBlocks);
		else if (args.Problem.Scheme.compare("g") == 0)
			biHarPb = new BiHarmonicMixedFormGlowinski_HHO<Dim>(mesh, testCase, hho, args.Discretization.BiharStabilization, saveMatrixBlocks);
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

		Timer lapSolverSetupTimer;

		if (args.Actions.SolveLinearSystem)
		{
			cout << endl;
			cout << "----------------------------------------------" << endl;
			cout << "-           Setup Laplacian solver           -" << endl;
			cout << "----------------------------------------------" << endl;

			int blockSizeForBlockSolver = args.Solver.BlockSize != -1 ? args.Solver.BlockSize : faceBasis->Size();
			Solver* diffSolver = SolverFactory<Dim>::CreateSolver(args, &biHarPb->DiffPb(), blockSizeForBlockSolver, out);


			cout << "Solver: " << *diffSolver << endl;
			IterativeSolver* iterativeSolver = dynamic_cast<IterativeSolver*>(diffSolver);
			if (iterativeSolver)
			{
				lapSolverSetupTimer.Start();
				if (Utils::ProgramArgs.Solver.SolverCode.compare("uamg") == 0 || Utils::ProgramArgs.Solver.PreconditionerCode.compare("uamg") == 0)
					iterativeSolver->Setup(biHarPb->DiffPb().A, biHarPb->DiffPb().A_T_T, biHarPb->DiffPb().A_T_ndF, biHarPb->DiffPb().A_ndF_ndF);
				else
					diffSolver->Setup(biHarPb->DiffPb().A);
				lapSolverSetupTimer.Stop();
			}
			else
			{
				lapSolverSetupTimer.Start();
				diffSolver->Setup(biHarPb->DiffPb().A);
				lapSolverSetupTimer.Stop();
			}
			cout << "Laplacian solver setup time: CPU = " << lapSolverSetupTimer.CPU() << ", elapsed = " << lapSolverSetupTimer.Elapsed() << endl << endl;

			biHarPb->SetDiffSolver(diffSolver);


			cout << "-------------------------------------" << endl;
			cout << "-     Solve biharmonic problem      -" << endl;
			cout << "-------------------------------------" << endl;

			cout << endl;

			Vector discreteExactSolution;
			double exactSolutionL2Norm = 0;
			if ((args.Actions.ComputeErrors || args.Solver.ComputeIterL2Error || args.Actions.Export.IterationL2Errors) && testCase->ExactSolution)
			{
				cout << "Projection of the exact solution on the discrete space..." << endl;
				discreteExactSolution = biHarPb->DiffPb().ReconstructSpace.Project(testCase->ExactSolution);
				exactSolutionL2Norm = biHarPb->DiffPb().ReconstructSpace.L2Norm(discreteExactSolution);
			}

			cout << "Computation of the right-hand side of the system..." << endl;

			biHarPb->SetDiffSolverTolerance(args.Solver.Tolerance);

			// Solve 1st problem with f as source
			Vector theta_f = biHarPb->FindCompatibleTheta();
			Vector lambda_f = biHarPb->Solve1stDiffProblem(theta_f);
			// Solve 2nd problem and extract the boundary (or normal derivative)
			Vector u_boundary_f = biHarPb->Solve2ndDiffProblem(lambda_f, true);

			// A(theta) -> -boundary(u) is s.p.d.
			// Note that because of the minus sign, u_boundary0 = -Af(theta_f)
			// We want to solve A(theta0) = u_boundary_f, in order to have
			// theta = theta_f+theta0 with
			//         A(theta) = Af(theta_f) + A(theta0) = -u_boundary_f + u_boundary_f = 0.
			Vector g_N = biHarPb->DiffPb().BoundarySpace.InnerProdWithBasis(testCase->NeumannBC.NeumannFunction);
			Vector b = u_boundary_f - g_N;

			//biHarPb->DiffPb().ExportReconstructedVectorToGMSH(biHarPb->Solve2ndDiffProblem(lambda_f, false), out, "u_f", args.Actions.Export.VisuTolerance, args.Actions.Export.VisuMaxRefinements);
			//biHarPb->DiffPb().ExportReconstructedVectorToGMSH(biHarPb->Solve1stDiffProblemWithZeroSource(u_boundary_f), out, "lambda_f", args.Actions.Export.VisuTolerance, args.Actions.Export.VisuMaxRefinements);

			//---------------------------//
			//     Explicit matrix A     //
			//---------------------------//

			DenseMatrix A; // computed explicitly only if explicit solver or export requested
			if (args.Solver.BiHarmonicSolverCode.compare("lu") == 0  || args.Solver.BiHarmonicSolverCode.compare("ch") == 0 ||
				(args.Solver.BiHarmonicSolverCode.compare("cg") == 0 && (args.Solver.BiHarmonicPreconditionerCode.compare("j") == 0 || args.Solver.BiHarmonicPreconditionerCode.compare("bj") == 0)) ||
				args.Actions.Export.LinearSystem || 
				args.Actions.UnitTests)
			{
				cout << "Explicit computation of the matrix..." << endl;
				A = biHarPb->Matrix();
				if (Utils::ProgramArgs.Actions.PrintDebug)
				{
					if (Dim <= 2)
						cout << "Matrix: " << std::scientific << std::setprecision(1) << endl << A << endl << endl;
					else
						cout << "Matrix: " << std::scientific << std::setprecision(1) << endl << A.topLeftCorner(3 * hho->nFaceUnknowns, 3 * hho->nFaceUnknowns) << endl << endl;
				}

				if (args.Actions.Export.LinearSystem)
				{
					cout << "Export linear system..." << endl;
					out.ExportMatrix(A, "matrix");
					out.ExportVector(b, "b");
				}

				if (args.Actions.UnitTests)
				{
					auto n = A.rows();
					//A = (A + A.transpose()) / 2.0;
					
					cout << "symmetry = " << (A-A.transpose()).norm() << endl;

					//double det = A.determinant();

					cout << "Computing eigenvalues..." << endl;
					Eigen::EigenSolver<DenseMatrix> es(A);
					Eigen::VectorXcd eigenvalues = es.eigenvalues();
					Eigen::MatrixXcd eigenvectors = es.eigenvectors();
					//cout << "Eigenvalues:" << endl << eigenvalues << endl;
					//Eigen::VectorXcd kernelVectorC = eigenvectors.col(n - 1);
					//cout << "---------------------" << endl << "Last eigenvector" << endl << kernelVector << endl;
					//cout << "---------------------" << endl << "Preceding eigenvector" << endl << eigenvectors.col(n - 2) << endl;
					//cout << "---------------------" << endl << "Preceding eigenvector" << endl << eigenvectors.col(n - 3) << endl;

					BiHarmonicMixedFormGlowinski_HHO<Dim>* biHarPb2 = static_cast<BiHarmonicMixedFormGlowinski_HHO<Dim>*>(biHarPb);

					int nZeroEigenvalues = 0;
					for (int i = 0; i < n; i++)
					{
						if (eigenvalues(i).real() < 1e-10)
							nZeroEigenvalues++;
					}

					if (nZeroEigenvalues > 0)
						Utils::Error(to_string(nZeroEigenvalues) + " zero eigenvalues found (amongst " + to_string(n) + "). The problem is not well-posed!");

					// Theta for the problematic element
					/*DenseMatrix theta = biHarPb2->DiffPb().HHOElement(largestBdryElem)->SolveCellUnknownsMatrix();
					cout << "theta: " << endl << theta << endl;
					// Theta_bF
					int nFaceUnknowns = hho->FaceBasis->Size();
					DenseMatrix theta_bF(theta.rows(), largestNBoundaryFaces* nFaceUnknowns);
					int i = 0;
					vector<Face<Dim>*> boundaryFaces = largestBdryElem->BoundaryFaces();
					for (int i=0; i< boundaryFaces.size(); i++)
					{
						Face<Dim>* f = boundaryFaces[i];
						theta_bF.middleCols(i * nFaceUnknowns, nFaceUnknowns) = theta.middleCols(largestBdryElem->LocalNumberOf(f) * nFaceUnknowns, nFaceUnknowns);
					}
					cout << "theta_bF: " << endl << theta_bF << endl;
					cout << endl;

					Interface<Dim> interf(boundaryFaces);
					list<set<Face<Dim>*>> coplanarSubsets = interf.CoplanarSubsets(Utils::NumericalZero);
					for (set<Face<Dim>*> subset : coplanarSubsets)
					{
						Vector kernelTheta = Vector::Zero(largestBdryElem->Faces.size() * nFaceUnknowns);
						auto iter = subset.begin();
						Face<Dim>* f1 = *iter;
						iter++;
						Face<Dim>* f2 = *iter;
						for (int i=0; i< nFaceUnknowns; i++)
							kernelTheta[largestBdryElem->LocalNumberOf(f1) * nFaceUnknowns + i] = 1;
						for (int i = 0; i < nFaceUnknowns; i++)
							kernelTheta[largestBdryElem->LocalNumberOf(f2) * nFaceUnknowns + i] = -f1->Measure() / f2->Measure();
						cout << "kernelTheta: " << endl << kernelTheta.transpose() << endl;
						cout << "theta * kernelTheta: " << endl << (theta * kernelTheta).transpose() << endl;
					}
					//FaceCollapsingStatus status = static_cast<PolyhedralMesh<Dim>*>.TryCollapseCoplanarFaces(interf);*/

					biHarPb2->print = true;
					for (int i = n - 1; i >= 0; i--)
					{
						if (eigenvalues(i).real() > 1e-10)
							continue;

						Vector kernelVector = eigenvectors.col(i).real();

						// Find the problematic element
						Element<Dim>* problematicElem = nullptr;
						for (int k=0; k<kernelVector.size(); ++k)
						{
							if (kernelVector(k) > 1e-5)
							{
								int boundaryFaceNumber = k / hho->nFaceUnknowns;
								int faceNumber = hho->nInteriorFaces + boundaryFaceNumber;
								auto it = std::find_if(mesh->Faces.begin(), mesh->Faces.end(), [&](Face<Dim>* f) { return f->Number == faceNumber; });
								Face<Dim>* f = *it;
								problematicElem = f->Element1;
								cout << "faceNumber: " << faceNumber << endl;
								cout << "problematic element: " << problematicElem->Number << endl;
								break;
							}
						}



						Vector allFacesValues = biHarPb->DiffPb().SkeletonSpace.ZeroVector();
						allFacesValues.tail(kernelVector.rows()) = kernelVector;
						GMSHMesh<Dim>::ExportToGMSH_Faces(static_cast<PolyhedralMesh<Dim>*>(mesh), hho->FaceBasis, allFacesValues, out.GetFilePathPrefix(), "kernelFaces_" + to_string(i));


						
						
						
						//cout << "kernelVector: " << endl << kernelVector.transpose() << endl;
						
						//Vector lambdaKernel = biHarPb->Solve1stDiffProblem_Homogeneous(kernelVector);
						//cout << "lambdaKernel: " << endl << lambdaKernel.transpose() << endl;
						//cout << lambdaKernel.norm() << endl;

						//DiffPb().ExportReconstructedVectorToGMSH(lambda, out, "lambda");
						//Vector solKernel = biHarPb->Solve2ndDiffProblem(lambdaKernel, false);
						//biHarPb->DiffPb().ExportReconstructedVectorToGMSH(solKernel, out, "solKernel");


						//Vector zero = -biHarPb->Solve2ndDiffProblem_Homogeneous(lambdaKernel);
						//cout << "Supposed to be 0: " << zero.norm() << endl;
						//break;
					}
					biHarPb2->print = false;
				}
			}

			//-------------------------------------//
			//            Create solver            //
			//-------------------------------------//

			cout << endl;

			Solver* biHarSolver = nullptr;
			if (args.Solver.BiHarmonicSolverCode.compare("cg") == 0)
			{
				ToleranceStrategy stgy = ToleranceStrategy::Fixed;
				if (Utils::ProgramArgs.Actions.Option2 == 1)
					stgy = ToleranceStrategy::DynamicFixedStep;
				else if (Utils::ProgramArgs.Actions.Option2 == 2)
					stgy = ToleranceStrategy::DynamicVariableStep;
				biHarSolver = new BiHarmonicCG(biHarPb, stgy, args.Solver.Tolerance2, 1e-3, args.Solver.Restart);
			}
			else if (args.Solver.BiHarmonicSolverCode.compare("fcg") == 0)
			{
				ToleranceStrategy stgy = ToleranceStrategy::Fixed;
				if (Utils::ProgramArgs.Actions.Option2 == 1)
					stgy = ToleranceStrategy::DynamicFixedStep;
				else if (Utils::ProgramArgs.Actions.Option2 == 2)
					stgy = ToleranceStrategy::DynamicVariableStep;
				biHarSolver = new BiHarmonicFCG(biHarPb, stgy, 1e-3, 4, args.Solver.Restart);
			}
			else if (args.Solver.BiHarmonicSolverCode.compare("gd") == 0)
				biHarSolver = new BiHarmonicGradientDescent(biHarPb);
			else if (args.Solver.BiHarmonicSolverCode.compare("lu") == 0)
				biHarSolver = new EigenLU();
			else if (args.Solver.BiHarmonicSolverCode.compare("ch") == 0)
				biHarSolver = new EigenCholesky();
			else
				Utils::FatalError("Unknown biharmonic solver '" + args.Solver.BiHarmonicSolverCode + "'");

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

				// Preconditioner
				if (args.Solver.BiHarmonicSolverCode.compare("cg") == 0 || args.Solver.BiHarmonicSolverCode.compare("fcg") == 0)
				{
					BiHarmonicCG* cg = static_cast<BiHarmonicCG*>(biHarIterSolver);
					if (args.Problem.Scheme.compare("g") == 0)
					{
						BiHarmonicMixedFormGlowinski_HHO<Dim>* gloScheme = static_cast<BiHarmonicMixedFormGlowinski_HHO<Dim>*>(biHarPb);

						ProgramArguments precSolverArgs;
						precSolverArgs.Solver.SolverCode = args.Solver.BiHarmonicPrecSolverCode;
						precSolverArgs.Solver.Tolerance = args.Solver.Tolerance;
						if (args.Solver.BiHarmonicPrecSolverTol > 0)
							precSolverArgs.Solver.Tolerance = args.Solver.BiHarmonicPrecSolverTol;
						precSolverArgs.Solver.MaxIterations = args.Solver.BiHarmonicPrecSolverMaxIter;
						Solver* precSolver = SolverFactory<Dim>::CreateSolver(precSolverArgs, blockSizeForBlockSolver, out);

						if (args.Solver.BiHarmonicPreconditionerCode.compare("j") == 0)
							cg->Precond = new DenseBlockJacobiPreconditioner(1);
						else if (args.Solver.BiHarmonicPreconditionerCode.compare("bj") == 0)
							cg->Precond = new DenseBlockJacobiPreconditioner(hho->nFaceUnknowns);
						else if (args.Solver.BiHarmonicPreconditionerCode.compare("s") == 0)
							cg->Precond = new BiharPatchPreconditioner<Dim>(*gloScheme, BiharPatchPreconditioner<Dim>::Type::SingleFaceNeighbourhood, 0, args.Solver.NeighbourhoodDepth, precSolver, false);
						else if (args.Solver.BiHarmonicPreconditionerCode.compare("ds") == 0)
							cg->Precond = new BiharPatchPreconditioner<Dim>(*gloScheme, BiharPatchPreconditioner<Dim>::Type::SingleFaceNeighbourhood, 0, args.Solver.NeighbourhoodDepth, precSolver, true);
						else if (args.Solver.BiHarmonicPreconditionerCode.compare("p") == 0)
							cg->Precond = new BiharPatchPreconditioner<Dim>(*gloScheme, BiharPatchPreconditioner<Dim>::Type::FacePatchNeighbourhood, args.Solver.PatchSize, args.Solver.NeighbourhoodDepth, precSolver, false);
						else if (args.Solver.BiHarmonicPreconditionerCode.compare("dp") == 0)
							cg->Precond = new BiharPatchPreconditioner<Dim>(*gloScheme, BiharPatchPreconditioner<Dim>::Type::FacePatchNeighbourhood, args.Solver.PatchSize, args.Solver.NeighbourhoodDepth, precSolver, true);
						else if (args.Solver.BiHarmonicPreconditionerCode.compare("no") == 0)
							cg->Precond = new IdentityPreconditioner();
						else
							Utils::FatalError("Unknown preconditioner. Check -bihar-prec " + args.Solver.BiHarmonicPreconditionerCode + ".");
					}
					else if (args.Problem.Scheme.compare("f") == 0)
						cg->Precond = new IdentityPreconditioner();
				}
			}

			cout << "Solver: " << *biHarSolver << endl << endl;

			Timer setupTimer;
			Timer solvingTimer;
			Timer totalTimer;
			totalTimer.Start();

			//------------------------------------//
			//            Setup solver            //
			//------------------------------------//

			setupTimer.Start();

			if (biHarIterSolver)
			{
				if (args.Solver.BiHarmonicSolverCode.compare("cg") == 0)
				{
					BiHarmonicCG* cg = static_cast<BiHarmonicCG*>(biHarIterSolver);
					cout << "Setup preconditioner..." << endl;
					cg->Precond->Setup(A);
				}
				else if (args.Solver.BiHarmonicSolverCode.compare("fcg") == 0)
				{
					BiHarmonicFCG* fcg = static_cast<BiHarmonicFCG*>(biHarIterSolver);
					cout << "Setup preconditioner..." << endl;
					fcg->Precond->Setup(A);
				}

				cout << "Setup solver..." << endl;
				// Give the Laplacian matrix so that the Work Units are computed with respect to it
				biHarIterSolver->Setup(biHarPb->DiffPb().A);
			}
			else // direct solver
			{
				cout << "Factorization..." << endl;
				biHarSolver->Setup(A);
			}

			setupTimer.Stop();

			cout << endl;

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

			
			//--------------------------------------------------------------------------//

			if (!biHarIterSolver)
			{
				Utils::Empty(A);
				cout << "Solve linear system..." << endl;
			}

			solvingTimer.Start();
			Vector theta_0 = biHarSolver->Solve(b);
			solvingTimer.Stop();

			totalTimer.Stop();

			cout << "Laplacian solver:" << endl << endl;
			int sizeTime = 12;
			cout << "        |   CPU time   | Elapsed time " << endl;
			cout << "---------------------------------------" << endl;
			cout << "Setup   | " << setw(sizeTime) << lapSolverSetupTimer.CPU() << " | " << setw(sizeTime) << lapSolverSetupTimer.Elapsed() << endl;
			cout << "        | " << setw(sizeTime - 2) << lapSolverSetupTimer.CPU().InSeconds() << " s | " << setw(sizeTime - 2) << lapSolverSetupTimer.Elapsed().InSeconds() << " s " << endl;
			cout << "---------------------------------------" << endl << endl;

			cout << "Biharmonic solver:" << endl << endl;
			SolverFactory<Dim>::PrintStats(biHarSolver, setupTimer, solvingTimer, totalTimer);

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

			if (args.Actions.Export.SolutionToGMSH && args.Discretization.Mesher.compare("gmsh") == 0)
				biHarPb->DiffPb().ExportReconstructedVectorToGMSH(reconstructedSolution, out, "solution", args.Actions.Export.VisuTolerance, args.Actions.Export.VisuMaxRefinements);

			//----------------------//
			//       L2 error       //
			//----------------------//

			if (args.Actions.ComputeErrors)
			{
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

				if (testCase->ExactSolution)
					cout << "h = " << mesh->H() << endl;
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