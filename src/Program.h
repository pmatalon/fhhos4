#pragma once
#include <iomanip>
#include "ProgramArguments.h"
#include "DG/Diffusion_DG.h"
#include "HHO/Diffusion_HHO.h"
#include "Mesher/AllMeshes.h"
#include "TestCases/TestCaseFactory.h"
#include "Utils/Timer.h"
#include "Solver/AllSolvers.h"
using namespace std;

class Program
{
public:
	Program() {}
	virtual void Start(ProgramArguments& args) = 0;
	virtual ~Program() {}
};

template <int Dim>
class ProgramDim : public Program
{
public:
	ProgramDim() : Program() {}

	void Start(ProgramArguments& args)
	{
		Utils::ProgramArgs = args;

		Timer totalTimer;
		totalTimer.Start();

#ifdef SMALL_INDEX
		cout << "Index type: int" << endl;
#else
		cout << "Index type: size_t" << endl;
#endif
		cout << "Shared memory parallelism: " << (BaseParallelLoop::GetDefaultNThreads() == 1 ? "sequential execution" : to_string(BaseParallelLoop::GetDefaultNThreads()) + " threads") << endl;

		GaussLegendre::Init();

		//-------------------------------------------------------------------------------------------//
		//   Test case defining the source function, boundary conditions and diffusion coefficient   //
		//-------------------------------------------------------------------------------------------//

		TestCase<Dim>* testCase = TestCaseFactory<Dim>::Create(args.Problem);
		
		//----------//
		//   Mesh   //
		//----------//

		cout << "-----------------------------------------------------------" << endl;
		cout << "-                   Mesh construction                     -" << endl;
		cout << "-----------------------------------------------------------" << endl;

		Mesh<Dim>::SetDirectories();
		GMSHMesh<Dim>::GMSHLogEnabled = args.Actions.GMSHLogEnabled;
		GMSHMesh<Dim>::UseCache = args.Actions.UseCache;

		Mesh<Dim>* mesh = BuildMesh(args, testCase);
		if (args.Discretization.Mesher.compare("gmsh") == 0)
			GMSHMesh<Dim>::CloseGMSH();

		cout << "Mesh storage > " << Utils::MemoryString(mesh->MemoryUsage()) << endl;


		mesh->SetDiffusionField(&testCase->DiffField);
		mesh->SetBoundaryConditions(&testCase->BC);

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
				if (!Utils::IsRefinementStrategy(args.Solver.MG.CoarseningStgy))
				{
					/*Mesh<Dim>* m = mesh;
					while (true)
					{
						cout << "Coarsening..." << endl;
						m->CoarsenMesh(args.Solver.MG.CoarseningStgy, args.Solver.MG.CoarseningFactor);
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
					mesh->CoarsenMesh(args.Solver.MG.CoarseningStgy, args.Solver.MG.FaceCoarseningStgy, args.Solver.MG.CoarseningFactor);
					/*cout << "Export..." << endl;
					mesh->CoarseMesh->ExportFacesToMatlab(args.OutputDirectory + "/coarse1.dat");
					mesh->CoarseMesh->ExportElementCentersToMatlab(args.OutputDirectory + "/elem_coarse1.m");*/
					cout << "Sanity check..." << endl;
					mesh->SanityCheck();
					// 2nd coarsening
					cout << "Coarsening..." << endl;
					mesh->CoarseMesh->CoarsenMesh(args.Solver.MG.CoarseningStgy, args.Solver.MG.FaceCoarseningStgy, args.Solver.MG.CoarseningFactor);
					/*cout << "Export..." << endl;
					mesh->CoarseMesh->CoarseMesh->ExportFacesToMatlab(args.OutputDirectory + "/coarse2.dat");
					mesh->CoarseMesh->CoarseMesh->ExportElementCentersToMatlab(args.OutputDirectory + "/elem_coarse2.m");*/
					cout << "Sanity check..." << endl;
					mesh->CoarseMesh->SanityCheck();
					// 3rd coarsening
					cout << "Coarsening..." << endl;
					mesh->CoarseMesh->CoarseMesh->CoarsenMesh(args.Solver.MG.CoarseningStgy, args.Solver.MG.FaceCoarseningStgy, args.Solver.MG.CoarseningFactor);
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

		//----------------------//
		//       Assembly       //
		//----------------------//

		Problem<Dim>* problem = nullptr;
		
		if (args.Discretization.Method.compare("dg") == 0)
		{
			FunctionalBasis<Dim>* basis = new FunctionalBasis<Dim>(args.Discretization.BasisCode, args.Discretization.PolyDegree, args.Discretization.UsePolynomialSpaceQ);
			problem = new Diffusion_DG<Dim>(mesh, testCase, args.OutputDirectory, basis, args.Discretization.PenalizationCoefficient);
		}
		else if (args.Discretization.Method.compare("hho") == 0)
		{
			FunctionalBasis<Dim>* reconstructionBasis = new FunctionalBasis<Dim>(args.Discretization.BasisCode, args.Discretization.PolyDegree, args.Discretization.UsePolynomialSpaceQ);
			FunctionalBasis<Dim>* cellBasis = new FunctionalBasis<Dim>(args.Discretization.BasisCode, args.Discretization.PolyDegree - 1, args.Discretization.UsePolynomialSpaceQ);
			FunctionalBasis<Dim - 1>* faceBasis = new FunctionalBasis<Dim - 1>(args.Discretization.BasisCode, args.Discretization.PolyDegree - 1, args.Discretization.UsePolynomialSpaceQ);

			HHOParameters<Dim>* hho = new HHOParameters<Dim>(mesh, args.Discretization.Stabilization, reconstructionBasis, cellBasis, faceBasis);

			bool saveMatrixBlocks = args.Solver.SolverCode.compare("camg") == 0 || args.Solver.SolverCode.compare("fcgcamg") == 0;
			problem = new Diffusion_HHO<Dim>(mesh, testCase, hho, args.Discretization.StaticCondensation, saveMatrixBlocks, args.OutputDirectory);
		}
		else
			Utils::FatalError("Unknown discretization.");

		cout << endl;
		cout << "----------------------------------------------------------" << endl;
		cout << "-                       Assembly                         -" << endl;
		cout << "----------------------------------------------------------" << endl;
		Timer assemblyTimer;
		assemblyTimer.Start();

		problem->Assemble(args.Actions);
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

			int blockSizeForBlockSolver = -1;
			if (args.Solver.BlockSize != -1)
				blockSizeForBlockSolver = args.Solver.BlockSize;
			else
			{
				if (args.Discretization.Method.compare("dg") == 0)
				{
					Diffusion_DG<Dim>* dgPb = static_cast<Diffusion_DG<Dim>*>(problem);
					blockSizeForBlockSolver = dgPb->Basis->Size();
				}
				else if (args.Discretization.Method.compare("hho") == 0)
				{
					Diffusion_HHO<Dim>* hhoPb = static_cast<Diffusion_HHO<Dim>*>(problem);
					blockSizeForBlockSolver = hhoPb->HHO->FaceBasis->Size();
				}
			}

			Solver* solver = CreateSolver(args, problem, blockSizeForBlockSolver);
			problem->SystemSolution = Solve(solver, problem, args.Solver.InitialGuessCode);

			if (args.Actions.ExportErrorToGMSH)
			{
				IterativeSolver* iterativeSolver = dynamic_cast<IterativeSolver*>(solver);
				if (iterativeSolver)
					problem->ExportErrorToGMSH(iterativeSolver->ExactSolution - problem->SystemSolution);
			}

			delete solver;

			//-----------------------------//
			//       Solution export       //
			//-----------------------------//

			if (args.Discretization.Method.compare("dg") == 0)
			{
				if (args.Actions.ExportSolutionVectors)
					problem->ExportSolutionVector();
			}
			else if (args.Discretization.Method.compare("hho") == 0)
			{
				Diffusion_HHO<Dim>* hhoPb = static_cast<Diffusion_HHO<Dim>*>(problem);

				if (args.Actions.ExportSolutionVectors || args.Actions.ExportSolutionToGMSH || testCase->ExactSolution)
				{
					cout << "----------------------------------------------------------" << endl;
					cout << "-                     Post-processing                    -" << endl;
					cout << "----------------------------------------------------------" << endl;

					hhoPb->ReconstructHigherOrderApproximation();
					if (args.Actions.ExportSolutionVectors)
					{
						if (args.Discretization.StaticCondensation)
							hhoPb->ExtractTraceSystemSolution();
						hhoPb->ExtractHybridSolution();
						hhoPb->ExportSolutionVector();
					}
				}

				if (args.Actions.ExportMeshToMatlab)
				{
					mesh->ExportToMatlab(args.OutputDirectory);
					mesh->ExportToMatlab2(args.OutputDirectory + "/mesh.m");
				}
			}

			if (args.Actions.ExportSolutionToGMSH && args.Discretization.Mesher.compare("gmsh") == 0)
				problem->ExportSolutionToGMSH();

			//----------------------//
			//       L2 error       //
			//----------------------//

			if (testCase->ExactSolution)
			{
				double error = problem->L2Error(testCase->ExactSolution);
				cout << endl << "L2 Error = " << std::scientific << error << endl;

				if (args.Discretization.Method.compare("hho") == 0)
				{
					Diffusion_HHO<Dim>* hhoPb = static_cast<Diffusion_HHO<Dim>*>(problem);
					hhoPb->AssertSchemeConvergence(error);
				}
			}
		}

		//--------------------------//
		//       Deallocation       //
		//--------------------------//
		
		delete mesh;
		if (args.Discretization.Method.compare("dg") == 0)
		{
			Diffusion_DG<Dim>* dgPb = static_cast<Diffusion_DG<Dim>*>(problem);
			delete dgPb->Basis;
		}
		else if (args.Discretization.Method.compare("hho") == 0)
		{
			Diffusion_HHO<Dim>* hhoPb = static_cast<Diffusion_HHO<Dim>*>(problem);
			delete hhoPb->HHO->CellBasis;
			delete hhoPb->HHO->FaceBasis;
			delete hhoPb->HHO->ReconstructionBasis;
			delete hhoPb->HHO;
		}
		delete problem;
		delete testCase;
		GaussLegendre::Free();

		totalTimer.Stop();
		cout << endl << "Total time: CPU = " << totalTimer.CPU() << ", elapsed = " << totalTimer.Elapsed() << endl;
	}

private:
	Mesh<Dim>* BuildMesh(ProgramArguments& args, TestCase<Dim>* testCase) { return nullptr; }

	Solver* CreateSolver(const ProgramArguments& args, Problem<Dim>* problem, int blockSize)
	{
		Solver* solver = nullptr;
		string GaussSeidelRelaxation1 = "The relaxation parameter of Gauss-Seidel is 1. Delete -relax arg to remove this warning, or use -s sor.";

		if (args.Solver.SolverCode.compare("lu") == 0)
			solver = new EigenSparseLU();
		else if (args.Solver.SolverCode.compare("j") == 0)
			solver = new BlockJacobi(1, args.Solver.RelaxationParameter);
		else if (args.Solver.SolverCode.compare("gs") == 0)
		{
			if (args.Solver.RelaxationParameter != 1)
				Utils::Warning(GaussSeidelRelaxation1);
			solver = new GaussSeidel(Direction::Forward);
		}
		else if (args.Solver.SolverCode.compare("rgs") == 0)
		{
			if (args.Solver.RelaxationParameter != 1)
				Utils::Warning(GaussSeidelRelaxation1);
			solver = new GaussSeidel(Direction::Backward);
		}
		else if (args.Solver.SolverCode.compare("sgs") == 0)
		{
			if (args.Solver.RelaxationParameter != 1)
				Utils::Warning(GaussSeidelRelaxation1);
			solver = new GaussSeidel(Direction::Symmetric);
		}
		else if (args.Solver.SolverCode.compare("sor") == 0)
			solver = new BlockSOR(1, args.Solver.RelaxationParameter, Direction::Forward);
		else if (args.Solver.SolverCode.compare("rsor") == 0)
			solver = new BlockSOR(1, args.Solver.RelaxationParameter, Direction::Backward);
		else if (args.Solver.SolverCode.compare("ssor") == 0)
			solver = new BlockSOR(1, args.Solver.RelaxationParameter, Direction::Symmetric);
		else if (args.Solver.SolverCode.compare("bj") == 0)
			solver = new BlockJacobi(blockSize, args.Solver.RelaxationParameter);
		else if (args.Solver.SolverCode.compare("bj23") == 0)
			solver = new BlockJacobi(blockSize, 2.0 / 3.0);
		else if (args.Solver.SolverCode.compare("bgs") == 0)
		{
			if (args.Solver.RelaxationParameter != 1)
				Utils::Warning(GaussSeidelRelaxation1);
			solver = new BlockSOR(blockSize, 1, Direction::Forward);
		}
		else if (args.Solver.SolverCode.compare("rbgs") == 0)
		{
			if (args.Solver.RelaxationParameter != 1)
				Utils::Warning(GaussSeidelRelaxation1);
			solver = new BlockSOR(blockSize, 1, Direction::Backward);
		}
		else if (args.Solver.SolverCode.compare("sbgs") == 0)
		{
			if (args.Solver.RelaxationParameter != 1)
				Utils::Warning(GaussSeidelRelaxation1);
			solver = new BlockSOR(blockSize, 1, Direction::Symmetric);
		}
		else if (args.Solver.SolverCode.compare("bsor") == 0)
			solver = new BlockSOR(blockSize, args.Solver.RelaxationParameter, Direction::Forward);
		else if (args.Solver.SolverCode.compare("rbsor") == 0)
			solver = new BlockSOR(blockSize, args.Solver.RelaxationParameter, Direction::Backward);
		else if (args.Solver.SolverCode.compare("sbsor") == 0)
			solver = new BlockSOR(blockSize, args.Solver.RelaxationParameter, Direction::Symmetric);
		else if (args.Solver.SolverCode.compare("eigencg") == 0)
			solver = new EigenCG();
		else if (args.Solver.SolverCode.rfind("cg", 0) == 0) // if SolverCode starts with "cg"
		{
			ConjugateGradient* cg = new ConjugateGradient();
			string preconditionerCode = args.Solver.SolverCode.substr(2, args.Solver.SolverCode.length() - 2);
			if (!preconditionerCode.empty())
			{
				ProgramArguments precondArgs = args;
				precondArgs.Solver.SolverCode = preconditionerCode;
				Solver* precondSolver = CreateSolver(precondArgs, problem, blockSize);
				cg->Precond = Preconditioner(dynamic_cast<IterativeSolver*>(precondSolver));
			}
			solver = cg;
		}
		else if (args.Solver.SolverCode.rfind("fcg", 0) == 0) // if SolverCode starts with "fcg"
		{
			FlexibleConjugateGradient* fcg = new FlexibleConjugateGradient(1);
			string preconditionerCode = args.Solver.SolverCode.substr(3, args.Solver.SolverCode.length() - 3);
			if (!preconditionerCode.empty())
			{
				ProgramArguments precondArgs = args;
				precondArgs.Solver.SolverCode = preconditionerCode;
				Solver* precondSolver = CreateSolver(precondArgs, problem, blockSize);
				fcg->Precond = Preconditioner(dynamic_cast<IterativeSolver*>(precondSolver));
			}
			solver = fcg;
		}
		else if (args.Solver.SolverCode.compare("agmg") == 0)
			solver = new AGMG();
		else if (args.Solver.SolverCode.compare("mg") == 0)
		{
			if (args.Discretization.StaticCondensation)
			{
				Diffusion_HHO<Dim>* hhoProblem = dynamic_cast<Diffusion_HHO<Dim>*>(problem);

				FunctionalBasis<Dim>* cellInterpolationBasis;
				if (args.Solver.MG.CellInterpolationCode == 1)
					cellInterpolationBasis = hhoProblem->HHO->ReconstructionBasis;
				else if (args.Solver.MG.CellInterpolationCode == 2)
					cellInterpolationBasis = hhoProblem->HHO->CellBasis;
				else
					assert(false);

				MultigridForHHO<Dim>* mg = new MultigridForHHO<Dim>(hhoProblem, args.Solver.MG.GMGProlong, cellInterpolationBasis, args.Solver.MG.WeightCode, args.Solver.MG.Levels);
				SetMultigridParameters(mg, args, blockSize);
				solver = mg;
			}
			else
				Utils::FatalError("The Multigrid for HHO (-s mg) only applicable on HHO discretization with static condensation.");
		}
		else if (args.Solver.SolverCode.compare("camg") == 0)
		{
			Diffusion_HHO<Dim>* hhoProblem = dynamic_cast<Diffusion_HHO<Dim>*>(problem);
			CondensedAMG* mg = new CondensedAMG(hhoProblem->HHO->nCellUnknowns, hhoProblem->HHO->nFaceUnknowns, 0.25, args.Solver.MG.CAMGFaceProlong, args.Solver.MG.CAMGProlong, args.Solver.MG.Levels);
			SetMultigridParameters(mg, args, blockSize);
			solver = mg;
		}
		else if (args.Solver.SolverCode.compare("aggregamg") == 0)
		{
			AggregAMG* mg = new AggregAMG(blockSize, 0.25, args.Solver.MG.Levels);
			SetMultigridParameters(mg, args, blockSize);
			mg->UseGalerkinOperator = 1;
			solver = mg;
		}
		else if (args.Solver.SolverCode.compare("hoaggregamg") == 0)
		{
			HighOrderAggregAMG* mg = new HighOrderAggregAMG(blockSize, 0.25, 2.0 / 3.0, 2);
			SetMultigridParameters(mg, args, blockSize);
			mg->UseGalerkinOperator = 1;
			solver = mg;
		}
		else
			Utils::FatalError("Unknown solver or not applicable.");

		IterativeSolver* iterativeSolver = dynamic_cast<IterativeSolver*>(solver);
		if (iterativeSolver)
		{
			iterativeSolver->Tolerance = args.Solver.Tolerance;
			iterativeSolver->MaxIterations = args.Solver.MaxIterations;
			iterativeSolver->PrintIterationResults = args.Solver.PrintIterationResults;
		}

		return solver;
	}

	void SetMultigridParameters(Multigrid* mg, const ProgramArguments& args, int blockSize)
	{
		mg->MatrixMaxSizeForCoarsestLevel = args.Solver.MG.MatrixMaxSizeForCoarsestLevel;
		mg->Cycle = args.Solver.MG.CycleLetter;
		mg->WLoops = args.Solver.MG.WLoops;
		mg->UseGalerkinOperator = args.Solver.MG.UseGalerkinOperator;
		mg->PreSmootherCode = args.Solver.MG.PreSmootherCode;
		mg->PostSmootherCode = args.Solver.MG.PostSmootherCode;
		mg->PreSmoothingIterations = args.Solver.MG.PreSmoothingIterations;
		mg->PostSmoothingIterations = args.Solver.MG.PostSmoothingIterations;
		mg->RelaxationParameter = args.Solver.RelaxationParameter;
		mg->BlockSizeForBlockSmoothers = blockSize;
		mg->CoarseLevelChangeSmoothingCoeff = args.Solver.MG.CoarseLevelChangeSmoothingCoeff;
		mg->CoarseLevelChangeSmoothingOperator = args.Solver.MG.CoarseLevelChangeSmoothingOperator;
		mg->CoarseningStgy = args.Solver.MG.CoarseningStgy;
		mg->FaceCoarseningStgy = args.Solver.MG.FaceCoarseningStgy;
		mg->CoarseningFactor = args.Solver.MG.CoarseningFactor;
		mg->ExportComponents = args.Actions.ExportMultigridComponents;

		// Coarse solver
		ProgramArguments argsCoarseSolver;
		argsCoarseSolver.Solver.SolverCode = args.Solver.MG.CoarseSolverCode;
		argsCoarseSolver.Solver.MaxIterations = 200;
		argsCoarseSolver.Solver.Tolerance = args.Solver.Tolerance;
		argsCoarseSolver.Solver.PrintIterationResults = false;
		if (Utils::EndsWith(args.Solver.MG.CoarseSolverCode, "aggregamg"))
		{
			argsCoarseSolver.Solver.MG.CoarseningStgy = CoarseningStrategy::AgglomerationCoarseningByFaceNeighbours;
			argsCoarseSolver.Solver.MG.CycleLetter = 'K';
		}
		mg->CoarseSolver = CreateSolver(argsCoarseSolver, nullptr, 1);
	}

	Vector Solve(Solver* solver, Problem<Dim>* problem, string initialGuessCode)
	{
		Timer setupTimer;
		Timer solvingTimer;
		Timer totalTimer;
		totalTimer.Start();

		cout << "Solver: " << *solver << endl << endl;

		Vector x;
		IterativeSolver* iterativeSolver = dynamic_cast<IterativeSolver*>(solver);
		if (iterativeSolver != nullptr)
		{
			iterativeSolver->ComputeExactSolution = Utils::ProgramArgs.Actions.ExportErrorToGMSH || problem->A.rows() <= 2000;

			setupTimer.Start();
			if (Utils::ProgramArgs.Solver.SolverCode.compare("camg") == 0 || Utils::ProgramArgs.Solver.SolverCode.compare("fcgcamg") == 0)
			{
				Diffusion_HHO<Dim>* hhoPb = static_cast<Diffusion_HHO<Dim>*>(problem);
				iterativeSolver->Setup(hhoPb->A, hhoPb->A_T_T, hhoPb->A_T_ndF, hhoPb->A_ndF_ndF);
			}
			else
				solver->Setup(problem->A);
			setupTimer.Stop();

			cout << "Solving..." << endl;
			solvingTimer.Start();
			x = iterativeSolver->Solve(problem->b, initialGuessCode);
			solvingTimer.Stop();
			cout << iterativeSolver->IterationCount << " iterations." << endl << endl;
		}
		else
		{
			setupTimer.Start();
			solver->Setup(problem->A);
			setupTimer.Stop();

			cout << "Solving..." << endl;
			solvingTimer.Start();
			x = solver->Solve(problem->b);
			solvingTimer.Stop();
			cout << endl;
		}
		totalTimer.Stop();

		int sizeTime = 12;
		int sizeWork = 8;
		int sizeMatVec = 8;

		auto oneFineMatVec = Cost::MatVec(problem->A);

		cout << "        |   CPU time   | Elapsed time ";
		if (iterativeSolver != nullptr)
			cout << "|  MFlops  |  MatVec  ";
		cout << endl;
		cout << "---------------------------------------";
		if (iterativeSolver != nullptr)
			cout << "---------------------";
		cout << endl;

		cout << "Setup   | " << setw(sizeTime) << setupTimer.CPU()                  << " | " << setw(sizeTime) << setupTimer.Elapsed();
		if (iterativeSolver != nullptr)
			cout << " | " << setw(sizeWork) << (iterativeSolver->SetupComputationalWork / (size_t)1e6) << " | " << setw(sizeMatVec) << (iterativeSolver->SetupComputationalWork / oneFineMatVec);
		cout << endl;
		cout << "        | " << setw(sizeTime-3) << setupTimer.CPU().InMilliseconds   << " ms | " << setw(sizeTime-3) << setupTimer.Elapsed().InMilliseconds << " ms ";
		if (iterativeSolver != nullptr)
			cout << "| " << setw(sizeWork) << " " << " | " << setw(sizeMatVec);
		cout << endl;
		cout << "---------------------------------------";
		if (iterativeSolver != nullptr)
			cout << "---------------------";
		cout << endl;

		cout << "Solving | " << setw(sizeTime) << solvingTimer.CPU()                  <<    " | " << setw(sizeTime) << solvingTimer.Elapsed();
		if (iterativeSolver != nullptr)
			cout << " | " << setw(sizeWork) << (iterativeSolver->SolvingComputationalWork / (size_t)1e6) << " | " << setw(sizeMatVec) << (iterativeSolver->SolvingComputationalWork / oneFineMatVec);
		cout << endl;
		cout << "        | " << setw(sizeTime-3) << solvingTimer.CPU().InMilliseconds << " ms | " << setw(sizeTime-3) << solvingTimer.Elapsed().InMilliseconds << " ms ";
		if (iterativeSolver != nullptr)
			cout << "| " << setw(sizeWork) << " " << " | " << setw(sizeMatVec);
		cout << endl;
		cout << "---------------------------------------";
		if (iterativeSolver != nullptr)
			cout << "---------------------";
		cout << endl;

		cout << "Total   | " << setw(sizeTime) << totalTimer.CPU()                  <<    " | " << setw(sizeTime) << totalTimer.Elapsed();
		if (iterativeSolver != nullptr)
			cout << " | " << setw(sizeWork) << ((iterativeSolver->SetupComputationalWork + iterativeSolver->SolvingComputationalWork)/(size_t)1e6) << " | " << setw(sizeMatVec) << ((iterativeSolver->SetupComputationalWork + iterativeSolver->SolvingComputationalWork) / oneFineMatVec);
		cout << endl;
		cout << "        | " << setw(sizeTime-3) << totalTimer.CPU().InMilliseconds   << " ms | " << setw(sizeTime-3) << totalTimer.Elapsed().InMilliseconds << " ms ";
		if (iterativeSolver != nullptr)
			cout << "| " << setw(sizeWork) << " " << " | " << setw(sizeMatVec);
		cout << endl;

		cout << endl;
		return x;
	}
};

template <>
Mesh<1>* ProgramDim<1>::BuildMesh(ProgramArguments& args, TestCase<1>* testCase)
{
	return new UniformMesh1D(args.Discretization.N);
}

template <>
Mesh<2>* ProgramDim<2>::BuildMesh(ProgramArguments& args, TestCase<2>* testCase)
{
	string geoCode = args.Problem.GeoCode;
	string mesher = args.Discretization.Mesher;
	BigNumber n = args.Discretization.N;
	BigNumber nx = args.Discretization.N;
	BigNumber ny = args.Discretization.Ny == -1 ? args.Discretization.N : args.Discretization.Ny;
	string meshCode = args.Discretization.MeshCode;
	double stretch = args.Discretization.Stretch;

	CoarseningStrategy refinementStgy = Utils::IsRefinementStrategy(args.Solver.MG.CoarseningStgy) ? args.Solver.MG.CoarseningStgy : CoarseningStrategy::GMSHSplittingRefinement;

	if (refinementStgy == CoarseningStrategy::BeyRefinement)
		Utils::FatalError("Bey's refinement method is only applicable to 3D tetrahedral meshes.");

	Mesh<2>* fineMesh = nullptr;
	//-------------------//
	//       Square      //
	//-------------------//
	if (geoCode.compare("square") == 0)
	{
		if (mesher.compare("inhouse") == 0)
		{
			if (Utils::IsRefinementStrategy(args.Solver.MG.CoarseningStgy))
				Utils::FatalError("Unmanaged refinement strategy.");

			if (meshCode.compare("cart") == 0)
				fineMesh = new Square_CartesianMesh(nx, ny);
			else if (meshCode.compare("cart-poly") == 0)
				fineMesh = new Square_CartesianPolygonalMesh(nx, ny);
			else if (meshCode.compare("stri") == 0)
				fineMesh = new Square_TriangularMesh(nx, ny);
			else if (meshCode.compare("quad") == 0)
				fineMesh = new Square_QuadrilateralMesh(nx, ny, stretch);
#ifdef CGAL_ENABLED
			else if (meshCode.compare("quad-poly") == 0)
				fineMesh = new Square_QuadrilateralAsPolygonalMesh(nx, ny, stretch);
#endif
			else if (meshCode.compare("tri") == 0)
				Utils::FatalError("The in-house mesher does not build unstructured meshes. Use '-mesher gmsh' or '-mesh stri' instead.");
			else
				Utils::FatalError("The requested mesh is not managed with this geometry.");
		}
#ifdef GMSH_ENABLED
		else if (mesher.compare("gmsh") == 0)
		{
			if (meshCode.compare("cart") == 0)
			{
				if (args.Solver.MG.CoarseningStgy == CoarseningStrategy::GMSHSplittingRefinement)
				{
					Mesh<2>* coarseMesh = new Square_GMSHCartesianMesh(n <= 16 ? 2 : 16);
					fineMesh = coarseMesh->RefineUntilNElements(nx*ny, refinementStgy);
				}
				else
					fineMesh = new Square_GMSHCartesianMesh(n);
			}
			else if (meshCode.compare("stri") == 0)
			{
				if (args.Solver.MG.CoarseningStgy == CoarseningStrategy::GMSHSplittingRefinement)
				{
					Mesh<2>* coarseMesh = new Square_GMSHTriangularMesh(n <= 16 ? 2 : 16);
					fineMesh = coarseMesh->RefineUntilNElements(2 * nx*ny, refinementStgy);
				}
				else
					fineMesh = new Square_GMSHTriangularMesh(n);
			}
			else if (meshCode.compare("tri") == 0)
			{
				if (args.Solver.MG.CoarseningStgy == CoarseningStrategy::GMSHSplittingRefinement)
				{
					Mesh<2>* coarseMesh = new Square_GMSHUnstructTriangularMesh(n <= 16 ? 2 : 16);
					fineMesh = coarseMesh->RefineUntilNElements(2 * nx*ny, refinementStgy);
				}
				else
					fineMesh = new Square_GMSHUnstructTriangularMesh(n);
			}
			else if (meshCode.compare("quad") == 0)
			{
				if (args.Solver.MG.CoarseningStgy == CoarseningStrategy::GMSHSplittingRefinement)
				{
					Mesh<2>* coarseMesh = new Square_GMSHQuadrilateralMesh(n <= 16 ? 2 : 16);
					fineMesh = coarseMesh->RefineUntilNElements(nx*ny, refinementStgy);
				}
				else
					fineMesh = new Square_GMSHQuadrilateralMesh(n);
			}
			else
				Utils::FatalError("The requested mesh is not managed with this geometry.");
		}
#endif // GMSH_ENABLED
		else
			Utils::FatalError("Unknown mesher. Check -mesher argument.");
	}
	//-------------------------------//
	//       Square 4 quadrants      //
	//-------------------------------//
	else if (geoCode.compare("square4quadrants") == 0)
	{
		bool with4quadrants = true;
		if (mesher.compare("inhouse") == 0)
		{
			if (Utils::IsRefinementStrategy(args.Solver.MG.CoarseningStgy))
				Utils::FatalError("Unmanaged refinement strategy.");

			if (meshCode.compare("cart") == 0)
				fineMesh = new Square_CartesianMesh(nx, ny, with4quadrants);
			else if (meshCode.compare("cart-poly") == 0)
				fineMesh = new Square_CartesianPolygonalMesh(nx, ny, with4quadrants);
			else if (meshCode.compare("stri") == 0)
				fineMesh = new Square_TriangularMesh(nx, ny, with4quadrants);
			else if (meshCode.compare("tri") == 0)
				Utils::FatalError("The in-house mesher does not build unstructured meshes. Use '-mesh stri' or '-mesher gmsh' instead.");
			else
				Utils::FatalError("The requested mesh is not managed with this geometry.");
		}
#ifdef GMSH_ENABLED
		else if (mesher.compare("gmsh") == 0)
		{
			if (meshCode.compare("cart") == 0)
			{
				if (args.Solver.MG.CoarseningStgy == CoarseningStrategy::GMSHSplittingRefinement)
				{
					Mesh<2>* coarseMesh = new Square4quadrants_GMSHCartesianMesh(n <= 16 ? 4 : 16);
					fineMesh = coarseMesh->RefineUntilNElements(nx*ny, refinementStgy);
				}
				else
					fineMesh = new Square4quadrants_GMSHCartesianMesh(n);
			}
			else if (meshCode.compare("tri") == 0)
			{
				if (args.Solver.MG.CoarseningStgy == CoarseningStrategy::GMSHSplittingRefinement)
				{
					Mesh<2>* coarseMesh = new Square4quadrants_GMSHUnstructTriangularMesh(n <= 16 ? 4 : 16);
					fineMesh = coarseMesh->RefineUntilNElements(2 * nx*ny, refinementStgy);
				}
				else
					fineMesh = new Square4quadrants_GMSHUnstructTriangularMesh(n);
			}
			else if (meshCode.compare("stri") == 0)
			{
				if (args.Solver.MG.CoarseningStgy == CoarseningStrategy::GMSHSplittingRefinement)
				{
					Mesh<2>* coarseMesh = new Square4quadrants_GMSHTriangularMesh(n <= 16 ? 4 : 16);
					fineMesh = coarseMesh->RefineUntilNElements(2 * nx*ny, refinementStgy);
				}
				else
					fineMesh = new Square4quadrants_GMSHTriangularMesh(n);
			}
			else if (meshCode.compare("quad") == 0)
			{
				if (args.Solver.MG.CoarseningStgy == CoarseningStrategy::GMSHSplittingRefinement)
				{
					Mesh<2>* coarseMesh = new Square4quadrants_GMSHQuadrilateralMesh(n <= 16 ? 4 : 16);
					fineMesh = coarseMesh->RefineUntilNElements(nx*ny, refinementStgy);
				}
				else
					fineMesh = new Square4quadrants_GMSHQuadrilateralMesh(n);
			}
			else
				Utils::FatalError("The requested mesh is not managed with this geometry.");
		}
#endif // GMSH_ENABLED
		else
			Utils::FatalError("Unknown mesher.");
	}
	//----------------------//
	//       GMSH file      //
	//----------------------//
	else
	{
		if (mesher.compare("inhouse") == 0)
			Utils::FatalError("The geometry is imported from a GMSH file, the mesher should be 'gmsh'. Use '-mesher gmsh' instead.");

#ifdef GMSH_ENABLED
		string filePath = geoCode;
		if (args.Solver.MG.CoarseningStgy == CoarseningStrategy::GMSHSplittingRefinement)
		{
			Mesh<2>* coarseMesh = new GMSHMesh<2>(testCase, filePath, args.Solver.MG.CoarseN);
			fineMesh = coarseMesh->RefineUntilNElements(2 * nx*ny, refinementStgy);
		}
		else
			fineMesh = new GMSHMesh<2>(testCase, filePath, n);
#endif // GMSH_ENABLED
	}
	
	if (!Utils::IsRefinementStrategy(args.Solver.MG.CoarseningStgy) && fineMesh->CoarseMesh)
		fineMesh->DeleteCoarseMeshes();

	return fineMesh;
}






template <>
Mesh<3>* ProgramDim<3>::BuildMesh(ProgramArguments& args, TestCase<3>* testCase)
{
	string geoCode = args.Problem.GeoCode;
	string mesher = args.Discretization.Mesher;
	BigNumber n = args.Discretization.N;
	BigNumber nx = args.Discretization.N;
	BigNumber ny = args.Discretization.Ny == -1 ? args.Discretization.N : args.Discretization.Ny;
	BigNumber nz = args.Discretization.Nz == -1 ? args.Discretization.N : args.Discretization.Nz;
	string meshCode = args.Discretization.MeshCode;
	CoarseningStrategy refinementStgy = args.Solver.MG.CoarseningStgy;

	Mesh<3>* fineMesh = nullptr;
	//------------------------//
	//          Cube          //
	//------------------------//
	if (geoCode.compare("cube") == 0)
	{
		if (mesher.compare("inhouse") == 0)
		{
			if (Utils::IsRefinementStrategy(args.Solver.MG.CoarseningStgy))
				Utils::FatalError("Unmanaged refinement strategy.");

			if (meshCode.compare("cart") == 0)
			{
				fineMesh = new Cube_CartesianMesh(nx, ny, nz);

				assert(fineMesh->Elements.size() == nx * ny*nz);
				if (nx == ny && ny == nz)
					assert(fineMesh->Faces.size() == 3 * n*n*(n + 1));
			}
			else if (meshCode.compare("stetra") == 0)
			{
				if (refinementStgy == CoarseningStrategy::StandardCoarsening || Utils::IsAlgebraic(args.Solver.SolverCode))
					fineMesh = new Cube_CartesianTetrahedralMesh(n);
				else
				{
					Mesh<3>* coarseMesh = new Cube_CartesianTetrahedralMesh(1);
					fineMesh = coarseMesh->RefineUntilNElements(6 * n*n*n, refinementStgy);
				}

				assert(fineMesh->Elements.size() == 6 * nx*ny*nz);
				if (nx == ny && ny == nz)
					assert(fineMesh->Faces.size() == 12 * n*n*n + 6 * n*n);
			}
			else if (meshCode.compare("tetra") == 0)
				Utils::FatalError("The in-house mesher does not build unstructured meshes. Use '-mesh stetra' or '-mesher gmsh'.");
			else
				Utils::FatalError("The requested mesh is not managed with this geometry.");
		}
#ifdef GMSH_ENABLED
		else if (mesher.compare("gmsh") == 0)
		{
			if (meshCode.compare("cart") == 0)
			{
				if (nx != ny || nx != nz)
					Utils::FatalError("Options -ny, -nz are not managed with this mesh");

				if (args.Solver.MG.CoarseningStgy == CoarseningStrategy::GMSHSplittingRefinement)
				{
					Mesh<3>* coarseMesh = new Cube_GMSHCartesianMesh(args.Solver.MG.CoarseN);
					fineMesh = coarseMesh->RefineUntilNElements(n*n*n, refinementStgy);

					assert(fineMesh->Elements.size() == n * n*n);
					assert(fineMesh->Faces.size() == 3 * n*n*(n + 1));
				}
				else
					fineMesh = new Cube_GMSHCartesianMesh(n);
			}
			else if (meshCode.compare("tetra") == 0)
			{
				if (nx != ny || nx != nz)
					Utils::FatalError("Options -ny, -nz are not managed with this mesh");

				if (Utils::IsRefinementStrategy(args.Solver.MG.CoarseningStgy))
				{
					Mesh<3>* coarseMesh = new Cube_GMSHTetrahedralMesh(args.Solver.MG.CoarseN);
					fineMesh = coarseMesh->RefineUntilNElements(6 * n*n*n, refinementStgy);

					assert(fineMesh->Elements.size() == 6 * n*n*n);
					assert(fineMesh->Faces.size() == 12 * n*n*n + 6 * n*n);
				}
				else
					fineMesh = new Cube_GMSHTetrahedralMesh(n);
			}
			else
				Utils::FatalError("The requested mesh is not managed with this geometry.");
		}
#endif // GMSH_ENABLED
	}
	//----------------------//
	//       GMSH file      //
	//----------------------//
	else
	{
		if (mesher.compare("inhouse") == 0)
			Utils::FatalError("The geometry is imported from a GMSH file, the mesher should be 'gmsh'. Use '-mesher gmsh' instead.");

#ifdef GMSH_ENABLED
		if (meshCode.compare("tetra") == 0)
		{
			string filePath = geoCode;
			if (nx != ny || nx != nz)
				Utils::FatalError("-ny, -ny not managed with this mesh");

			Mesh<3>* coarseMesh;
			if (Utils::IsRefinementStrategy(args.Solver.MG.CoarseningStgy))
			{
				if (refinementStgy == CoarseningStrategy::BeyRefinement)
					coarseMesh = new GMSHTetrahedralMesh(testCase, filePath, args.Solver.MG.CoarseN);
				else if (refinementStgy == CoarseningStrategy::GMSHSplittingRefinement)
					coarseMesh = new GMSHMesh<3>(testCase, filePath, args.Solver.MG.CoarseN);
				fineMesh = coarseMesh->RefineUntilNElements(6 * n*n*n, refinementStgy);
			}
			else
				fineMesh = new GMSHMesh<3>(testCase, filePath, n);
		}
		else
			Utils::FatalError("When the geometry is imported from a GMSH file, only unstructured tetrahedral meshing is allowed. Use '-mesh tetra' instead.");
#endif // GMSH_ENABLED
	}

	if (!Utils::IsRefinementStrategy(args.Solver.MG.CoarseningStgy) && fineMesh->CoarseMesh)
		fineMesh->DeleteCoarseMeshes();

	return fineMesh;
}
