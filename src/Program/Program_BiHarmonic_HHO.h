#pragma once
#include <iomanip>
#include "../ProgramArguments.h"
#include "../HHO/BiHarmonicMixedForm_HHO.h"
#include "../TestCases/BiHarmonic/BiHarTestCaseFactory.h"
#include "../Mesher/MeshFactory.h"
#include "../Solver/SolverFactory.h"
#include "../Solver/Krylov/BiHarmonicCG.h"
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
		BiHarmonicMixedForm_HHO<Dim>* biHarPb = new BiHarmonicMixedForm_HHO<Dim>(mesh, testCase, hho, saveMatrixBlocks);

		cout << endl;
		cout << "------------------------------------------------------" << endl;
		cout << "-          Assembly of diffusion problem             -" << endl;
		cout << "------------------------------------------------------" << endl;
		Timer assemblyTimer;
		assemblyTimer.Start();

		biHarPb->AssembleDiffPb();
		cout << "System storage: " << Utils::MemoryString(Utils::MemoryUsage(biHarPb->DiffPb().A) + Utils::MemoryUsage(biHarPb->DiffPb().b)) << endl;

		assemblyTimer.Stop();
		cout << endl << "Assembly time: CPU = " << assemblyTimer.CPU() << ", elapsed = " << assemblyTimer.Elapsed() << endl;

		//------------------------------------//
		//       Linear system solution       //
		//------------------------------------//

		if (args.Actions.SolveLinearSystem)
		{
			cout << endl;
			cout << "----------------------------------------------------------" << endl;
			cout << "-                   Solve linear system                  -" << endl;
			cout << "----------------------------------------------------------" << endl;

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

			cout << "Solve with Conjugate Gradient" << endl;

			BiHarmonicCG<Dim> cg(testCase, biHarPb->DiffPb(), diffSolver);
			cg.Tolerance = args.Solver.Tolerance;
			cg.MaxIterations = args.Solver.MaxIterations;
			Vector theta = cg.Solve();

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