#pragma once
#include "ProgramArguments.h"
#include "DG/Diffusion_DG.h"
#include "HHO/Diffusion_HHO.h"
#include "Mesh/AllMeshes.h"
#include "TestCases/TestCaseFactory.h"
#include "Utils/Action.h"
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
		Timer totalTimer;
		totalTimer.Start();

		GaussLegendre::Init();

		//----------//
		//   Mesh   //
		//----------//

		cout << "-----------------------------------------------------------" << endl;
		cout << "-                   Mesh construction                     -" << endl;
		cout << "-----------------------------------------------------------" << endl;
		Mesh<Dim>* mesh = BuildMesh(args);

		if ((args.Actions & Action::UnitTests) == Action::UnitTests)
		{
			// Unit tests
			Triangle::Test();
			Quadrilateral::Test();
			Polygon::Test();
			Tetrahedron::Test();
			TriangleIn3D::Test();

			mesh->SanityCheck();
			if (args.Discretization.N <= 2)
				cout << *mesh << endl << endl;

			/*if (args.Discretization.MeshCode.compare("gmsh-tri") == 0)
			{
				GMSHMesh<Dim>* gmshMesh = dynamic_cast<GMSHMesh<Dim>*>(mesh);
				gmshMesh->RenumberLikeMe();
			}*/
		}

		//-------------------------------------------------------------------------------------------//
		//   Test case defining the source function, boundary conditions and diffusion coefficient   //
		//-------------------------------------------------------------------------------------------//

		TestCase<Dim>* testCase = TestCaseFactory<Dim>::Create(args.Problem);

		mesh->SetDiffusionField(&testCase->DiffField);
		mesh->SetBoundaryConditions(&testCase->BC);

		//---------------------------//
		//   Coarsening unit tests   //
		//---------------------------//

		if ((args.Actions & Action::UnitTests) == Action::UnitTests && Dim == 2)
		{
			if (args.Solver.SolverCode.compare("mg") == 0 || args.Solver.SolverCode.compare("cgmg") == 0 || args.Solver.SolverCode.compare("fcgmg") == 0)
			{
				if (!Utils::IsRefinementStrategy(args.Solver.MG.CoarseningStgy))
				{
					mesh->ExportFacesToMatlab(args.OutputDirectory + "/fine.dat");
					mesh->ExportElementCentersToMatlab(args.OutputDirectory + "/elem_fine.m");

					// 1st coarsening
					mesh->CoarsenMesh(args.Solver.MG.CoarseningStgy);
					mesh->CoarseMesh->ExportFacesToMatlab(args.OutputDirectory + "/coarse1.dat");
					mesh->CoarseMesh->ExportElementCentersToMatlab(args.OutputDirectory + "/elem_coarse1.m");
					mesh->SanityCheck();
					// 2nd coarsening
					mesh->CoarseMesh->CoarsenMesh(args.Solver.MG.CoarseningStgy);
					mesh->CoarseMesh->CoarseMesh->ExportFacesToMatlab(args.OutputDirectory + "/coarse2.dat");
					mesh->CoarseMesh->CoarseMesh->ExportElementCentersToMatlab(args.OutputDirectory + "/elem_coarse2.m");
					mesh->SanityCheck();
					// 3rd coarsening
					mesh->CoarseMesh->CoarseMesh->CoarsenMesh(args.Solver.MG.CoarseningStgy);
					mesh->CoarseMesh->CoarseMesh->CoarseMesh->ExportFacesToMatlab(args.OutputDirectory + "/coarse3.dat");
					mesh->CoarseMesh->CoarseMesh->CoarseMesh->ExportElementCentersToMatlab(args.OutputDirectory + "/elem_coarse3.m");
					mesh->SanityCheck();
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

			problem = new Diffusion_HHO<Dim>(mesh, testCase, hho, args.Discretization.StaticCondensation, args.OutputDirectory);
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
			
		assemblyTimer.Stop();
		cout << endl << "Assembly time: CPU = " << assemblyTimer.CPU() << ", elapsed = " << assemblyTimer.Elapsed() << endl;

		//------------------------------------//
		//       Linear system solution       //
		//------------------------------------//

		if ((args.Actions & Action::SolveSystem) == Action::SolveSystem)
		{
			cout << endl;
			cout << "----------------------------------------------------------" << endl;
			cout << "-                 Linear system solution                 -" << endl;
			cout << "----------------------------------------------------------" << endl;

			int blockSizeForBlockSolver = 1;
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

			Solver* solver = CreateSolver(args, problem, blockSizeForBlockSolver);
			problem->SystemSolution = Solve(solver, problem->A, problem->b, args.Solver.InitialGuessCode);

			//-----------------------------//
			//       Solution export       //
			//-----------------------------//

			if (args.Discretization.Method.compare("dg") == 0)
			{
				if ((args.Actions & Action::ExtractSolution) == Action::ExtractSolution)
					problem->ExtractSolution();
			}
			else if (args.Discretization.Method.compare("hho") == 0)
			{
				Diffusion_HHO<Dim>* hhoPb = static_cast<Diffusion_HHO<Dim>*>(problem);
				if (args.Discretization.StaticCondensation && (args.Actions & Action::ExtractSolution) == Action::ExtractSolution)
					hhoPb->ExtractTraceSystemSolution();

				if (    (args.Actions & Action::ExtractSolution) == Action::ExtractSolution 
					|| ((args.Actions & Action::ComputeL2Error) == Action::ComputeL2Error && testCase->ExactSolution)
					||  (args.Actions & Action::ExportSolutionToGMSH) == Action::ExportSolutionToGMSH)
				{

					cout << "----------------------------------------------------------" << endl;
					cout << "-                     Post-processing                    -" << endl;
					cout << "----------------------------------------------------------" << endl;

					hhoPb->ReconstructHigherOrderApproximation();
					if ((args.Actions & Action::ExtractSolution) == Action::ExtractSolution)
					{
						hhoPb->ExtractHybridSolution();
						hhoPb->ExtractSolution();
					}

					if ((args.Actions & Action::ExportSolutionToGMSH) == Action::ExportSolutionToGMSH)
						hhoPb->ExportSolutionToGMSH();
				}

				if ((args.Actions & Action::ExportFaces) == Action::ExportFaces)
					mesh->ExportFacesToMatlab(args.OutputDirectory, true);
			}

			//----------------------//
			//       L2 error       //
			//----------------------//

			if ((args.Actions & Action::ComputeL2Error) == Action::ComputeL2Error && testCase->ExactSolution)
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
	Mesh<Dim>* BuildMesh(ProgramArguments& args) { return nullptr; }

	Solver* CreateSolver(const ProgramArguments& args, Problem<Dim>* problem, int blockSize)
	{
		Solver* solver = NULL;
		if (args.Solver.SolverCode.compare("mg") == 0 || args.Solver.SolverCode.compare("cgmg") == 0 || args.Solver.SolverCode.compare("fcgmg") == 0)
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

				MultigridForHHO<Dim>* mg = new MultigridForHHO<Dim>(hhoProblem, args.Solver.MG.ProlongationCode, cellInterpolationBasis, args.Solver.MG.WeightCode, args.Solver.MG.Levels);
				mg->MatrixMaxSizeForCoarsestLevel = args.Solver.MG.MatrixMaxSizeForCoarsestLevel;
				mg->Cycle = args.Solver.MG.CycleLetter;
				mg->WLoops = args.Solver.MG.WLoops;
				mg->UseGalerkinOperator = args.Solver.MG.UseGalerkinOperator;
				mg->PreSmootherCode = args.Solver.MG.PreSmootherCode;
				mg->PostSmootherCode = args.Solver.MG.PostSmootherCode;
				mg->PreSmoothingIterations = args.Solver.MG.PreSmoothingIterations;
				mg->PostSmoothingIterations = args.Solver.MG.PostSmoothingIterations;
				mg->RelaxationParameter = args.Solver.RelaxationParameter;
				mg->CoarseLevelChangeSmoothingCoeff = args.Solver.MG.CoarseLevelChangeSmoothingCoeff;
				mg->CoarseLevelChangeSmoothingOperator = args.Solver.MG.CoarseLevelChangeSmoothingOperator;
				mg->CoarseningStgy = args.Solver.MG.CoarseningStgy;
				mg->ExportMatrices = (args.Actions & Action::ExportMultigridMatrices) == Action::ExportMultigridMatrices;

				if (args.Solver.SolverCode.compare("mg") == 0)
					solver = mg;
				else if (args.Solver.SolverCode.compare("cgmg") == 0)
				{
					ConjugateGradient* cg = new ConjugateGradient();
					cg->Precond = Preconditioner(mg);
					solver = cg;
				}
				else if (args.Solver.SolverCode.compare("fcgmg") == 0)
				{
					FlexibleConjugateGradient* fcg = new FlexibleConjugateGradient(1);
					fcg->Precond = Preconditioner(mg);
					solver = fcg;
				}
			}
			else
				assert(false && "Multigrid only applicable on HHO discretization with static condensation.");
		}
		else if (args.Solver.SolverCode.compare("lu") == 0)
			solver = new EigenSparseLU();
		else if (args.Solver.SolverCode.compare("cg") == 0)
			solver = new ConjugateGradient();
		else if (args.Solver.SolverCode.compare("eigencg") == 0)
			solver = new EigenCG();
		else if (args.Solver.SolverCode.compare("fcg") == 0)
			solver = new FlexibleConjugateGradient();
		else if (args.Solver.SolverCode.compare("j") == 0)
			solver = new BlockJacobi(1, args.Solver.RelaxationParameter);
		else if (args.Solver.SolverCode.compare("sor") == 0 || args.Solver.SolverCode.compare("gs") == 0)
			solver = new BlockSOR(1, args.Solver.RelaxationParameter, Direction::Forward);
		else if (args.Solver.SolverCode.compare("rsor") == 0 || args.Solver.SolverCode.compare("rgs") == 0)
			solver = new BlockSOR(1, args.Solver.RelaxationParameter, Direction::Backward);
		else if (args.Solver.SolverCode.compare("bsor") == 0 || args.Solver.SolverCode.compare("bgs") == 0)
			solver = new BlockSOR(blockSize, args.Solver.RelaxationParameter, Direction::Forward);
		else if (args.Solver.SolverCode.compare("rbsor") == 0 || args.Solver.SolverCode.compare("rbgs") == 0)
			solver = new BlockSOR(blockSize, args.Solver.RelaxationParameter, Direction::Backward);
		else if (args.Solver.SolverCode.compare("bj") == 0)
			solver = new BlockJacobi(blockSize, args.Solver.RelaxationParameter);
		else if (args.Solver.SolverCode.compare("bj23") == 0)
			solver = new BlockJacobi(blockSize, 2.0/3.0);
#ifdef AGMG_ENABLED
		else if (args.Solver.solverCode.compare("agmg") == 0)
			solver = new AGMG(tolerance);
#endif // AGMG_ENABLED
		else
			assert(false && "Unknown solver or not applicable!");

		IterativeSolver* iterativeSolver = dynamic_cast<IterativeSolver*>(solver);
		if (iterativeSolver != nullptr)
		{
			iterativeSolver->Tolerance = args.Solver.Tolerance;
			iterativeSolver->MaxIterations = args.Solver.MaxIterations;
		}

		return solver;
	}

	Vector Solve(Solver* solver, const SparseMatrix& A, const Vector& b, string initialGuessCode)
	{
		Timer solverTimer;
		solverTimer.Start();

		cout << "Solver: " << *solver << endl << endl;

		Vector x;
		IterativeSolver* iterativeSolver = dynamic_cast<IterativeSolver*>(solver);
		if (iterativeSolver != nullptr)
		{
			iterativeSolver->ComputeExactSolution = A.rows() <= 2000;
			solver->Setup(A);
			cout << "Solving..." << endl;
			x = iterativeSolver->Solve(b, initialGuessCode);
			cout << iterativeSolver->IterationCount << " iterations." << endl;
		}
		else
		{
			solver->Setup(A);
			cout << "Solving..." << endl;
			x = solver->Solve(b);
		}

		delete solver;

		solverTimer.Stop();
		cout << "Solving time: CPU = " << solverTimer.CPU() << ", elapsed = " << solverTimer.Elapsed() << endl << endl;
		return x;
	}
};

template <>
Mesh<1>* ProgramDim<1>::BuildMesh(ProgramArguments& args)
{
	return new UniformMesh1D(args.Discretization.N);
}

template <>
Mesh<2>* ProgramDim<2>::BuildMesh(ProgramArguments& args)
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
			else if (meshCode.compare("quad-poly") == 0)
				fineMesh = new Square_QuadrilateralAsPolygonalMesh(nx, ny, stretch);
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
			Mesh<2>* coarseMesh = new GMSHMesh<2>(filePath, 2);
			fineMesh = coarseMesh->RefineUntilNElements(2 * nx*ny, refinementStgy);
		}
		else
			fineMesh = new GMSHMesh<2>(filePath, n);
#endif // GMSH_ENABLED
	}
	
	if (!Utils::IsRefinementStrategy(args.Solver.MG.CoarseningStgy) && fineMesh->CoarseMesh)
		fineMesh->DeleteCoarseMeshes();

	return fineMesh;
}






template <>
Mesh<3>* ProgramDim<3>::BuildMesh(ProgramArguments& args)
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
				if (refinementStgy == CoarseningStrategy::StandardCoarsening)
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
					Mesh<3>* coarseMesh = new Cube_GMSHCartesianMesh(2);
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
					Mesh<3>* coarseMesh = new Cube_GMSHTetrahedralMesh(2);
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
			if (refinementStgy == CoarseningStrategy::BeyRefinement)
				coarseMesh = new GMSHTetrahedralMesh(filePath, 2);
			else
				coarseMesh = new GMSHMesh<3>(filePath);
			fineMesh = coarseMesh->RefineUntilNElements(6 * n*n*n, refinementStgy);
		}
		else
			Utils::FatalError("When the geometry is imported from a GMSH file, only unstructured tetrahedral meshing is allowed. Use '-mesh tetra' instead.");
#endif // GMSH_ENABLED
	}

	if (!Utils::IsRefinementStrategy(args.Solver.MG.CoarseningStgy) && fineMesh->CoarseMesh)
		fineMesh->DeleteCoarseMeshes();

	return fineMesh;
}
