#pragma once
#include <cstdio>
#include <iostream>
#include <functional>
#include "DG/Poisson_DG.h"
#include "HHO/Poisson_HHO.h"
#include "FunctionalBasis/FunctionalBasis.h"
#include "Mesh/1D/CartesianGrid1D.h"
#include "Mesh/2D/CartesianGrid2D.h"
#include "Mesh/3D/CartesianGrid3D.h"
#include "Mesh/2D/CartesianPolygonalMesh2D.h"
#include "Utils/Action.h"
#include "Utils/DiffusionPartition.h"
#include "Solver/MultigridForHHO.h"
#include "Solver/BlockJacobi.h"
#include "Solver/EigenCG.h"
#include "Solver/AGMG.h"
using namespace std;

class Program
{
public:
	Program() {}
	virtual void Start(string solution, double kappa1, double kappa2, BigNumber n, string discretization, string basisCode, int polyDegree, bool fullTensorization, 
		int penalizationCoefficient, bool staticCondensation, Action action, 
		int nMultigridLevels, int wLoops, bool useGalerkinOperator, string preSmootherCode, string postSmootherCode, int nPreSmoothingIterations, int nPostSmoothingIterations, 
		string outputDirectory, string solverCode, double solverTolerance) = 0;
};

template <int Dim>
class ProgramDim : public Program
{
public:
	ProgramDim() : Program() {}

	void Start(string solution, double kappa1, double kappa2, BigNumber n, string discretization, string basisCode, int polyDegree, bool fullTensorization, 
		int penalizationCoefficient, bool staticCondensation, Action action, 
		int nMultigridLevels, int wLoops, bool useGalerkinOperator, string preSmootherCode, string postSmootherCode, int nPreSmoothingIterations, int nPostSmoothingIterations, 
		string outputDirectory, string solverCode, double solverTolerance)
	{
		//----------//
		//   Mesh   //
		//----------//

		Mesh<Dim>* mesh = BuildMesh(n);

		//mesh->CoarsenMesh(CoarseningStrategy::AgglomerationAndKeepFineFaces);
		//mesh = mesh->CoarseMesh;
		//cout << *mesh << endl << endl;
		//cout << "Coarse mesh" << endl << *(mesh->CoarseMesh) << endl << endl;

		//---------------------------------------------//
		//   Analytical solution and source function   //
		//---------------------------------------------//

		function<double(DomPoint)> exactSolution = NULL;
		SourceFunction* sourceFunction;

		function<bool(DomPoint)> isInPart1 = [](DomPoint p) { return p.X < 0.5; };
		DiffusionPartition diffusionPartition(isInPart1, kappa1, kappa2);

		if (Dim == 1)
		{
			if (solution.compare("sine") == 0)
			{
				exactSolution = [](DomPoint p)
				{
					double x = p.X;
					return sin(4 * M_PI * x) / (16 * pow(M_PI, 2));
				};
				sourceFunction = new SourceFunction1D([&diffusionPartition](double x) { return diffusionPartition.Coefficient(x) * sin(4 * M_PI * x); });
			}
			else if (solution.compare("poly") == 0)
			{
				exactSolution = [](DomPoint p)
				{
					double x = p.X;
					return x * (1 - x);
				};
				sourceFunction = new SourceFunction1D([&diffusionPartition](double x) { return diffusionPartition.Coefficient(x) * 2; });
				//sourceFunction = [](double x) { return (-1)*(-6 * x*pow(x - 1, 3) - 3 * pow(x, 3) * (2 * x - 2) - 18 * pow(x, 2) * pow(x - 1, 2)); };
			}
			else if (solution.compare("hetero") == 0)
			{
				exactSolution = [&diffusionPartition](DomPoint p)
				{
					double x = p.X;
					double alpha = diffusionPartition.Kappa1;
					double a1 = -1 / (2 * alpha);
					double a2 = -0.5;
					double b1 = (1 + 3 * alpha) / (2 * alpha*(1 + alpha));
					double b2 = -(alpha + 3) / (2 * (1 + alpha));
					if (diffusionPartition.IsInPart1(p))
						return 4 * a1 *pow(x, 2) + 2 * b1 * x;
					else
						return 4 * a2 * pow(x - 1, 2) + 2 * b2 * (x - 1);
				};
				sourceFunction = new SourceFunction1D([&diffusionPartition](double x) { return 4; });
			}
		}
		else if (Dim == 2)
		{
			if (solution.compare("sine") == 0)
			{
				exactSolution = [](DomPoint p)
				{
					double x = p.X;
					double y = p.Y;
					return sin(4 * M_PI * x)*sin(4 * M_PI * y);
				};
				sourceFunction = new SourceFunction2D([&diffusionPartition](double x, double y) { return diffusionPartition.Coefficient(x) * 2 * pow(4 * M_PI, 2) * sin(4 * M_PI * x)*sin(4 * M_PI * y); });
			}
			else if (solution.compare("poly") == 0)
			{
				exactSolution = [](DomPoint p)
				{
					double x = p.X;
					double y = p.Y;
					return x * (1 - x) * y*(1 - y);
				};
				sourceFunction = new SourceFunction2D([&diffusionPartition](double x, double y) { return diffusionPartition.Coefficient(x) * 2 * (y*(1 - y) + x * (1 - x)); });
			}
		}
		else if (Dim == 3)
		{
			if (solution.compare("sine") == 0)
			{
				exactSolution = [](DomPoint p)
				{
					double x = p.X;
					double y = p.Y;
					double z = p.Z;
					return sin(4 * M_PI * x)*sin(4 * M_PI * y)*sin(4 * M_PI * z);
				};
				sourceFunction = new SourceFunction3D([&diffusionPartition](double x, double y, double z) {  return diffusionPartition.Coefficient(x) * 3 * pow(4 * M_PI, 2) * sin(4 * M_PI * x)*sin(4 * M_PI * y)*sin(4 * M_PI * z); });
			}
			else if (solution.compare("poly") == 0)
			{
				exactSolution = [](DomPoint p)
				{
					double x = p.X;
					double y = p.Y;
					double z = p.Z;
					return x * (1 - x)*y*(1 - y)*z*(1 - z);
				};
				sourceFunction = new SourceFunction3D([&diffusionPartition](double x, double y, double z) { return diffusionPartition.Coefficient(x) * 2 * ((y*(1 - y)*z*(1 - z) + x * (1 - x)*z*(1 - z) + x * (1 - x)*y*(1 - y))); });
			}
		}

		//--------------------------------//
		//   Discretization and solving   //
		//--------------------------------//

		if (discretization.compare("dg") == 0)
		{
			Poisson_DG<Dim>* problem = new Poisson_DG<Dim>(solution, sourceFunction, diffusionPartition, outputDirectory);
			FunctionalBasis<Dim>* basis = new FunctionalBasis<Dim>(basisCode, polyDegree, fullTensorization);

			cout << endl;
			cout << "----------------------- Assembly -------------------------" << endl;
			problem->Assemble(mesh, basis, penalizationCoefficient, action);

			if ((action & Action::SolveSystem) == Action::SolveSystem)
			{
				cout << endl;
				cout << "------------------- Linear system resolution ------------------" << endl;

				Solver* solver = CreateSolver(solverCode, problem, solverTolerance, staticCondensation, nMultigridLevels, wLoops, useGalerkinOperator, preSmootherCode, postSmootherCode, nPreSmoothingIterations, nPostSmoothingIterations, basis->Size());
				cout << "Solver: " << *solver << endl << endl;
				solver->Setup(problem->A);
				problem->Solution = solver->Solve(problem->b);
				delete solver;
				if ((action & Action::ExtractSolution) == Action::ExtractSolution)
					problem->ExtractSolution();
				double error = L2::Error<Dim>(mesh, basis, problem->Solution, exactSolution);
				cout << "L2 Error = " << error << endl;
			}

			delete problem;
			delete basis;
		}
		else if (discretization.compare("hho") == 0)
		{
			FunctionalBasis<Dim>* reconstructionBasis = new FunctionalBasis<Dim>(basisCode, polyDegree, fullTensorization);
			FunctionalBasis<Dim>* cellBasis = new FunctionalBasis<Dim>(basisCode, polyDegree - 1, fullTensorization);
			FunctionalBasis<Dim-1>* faceBasis = new FunctionalBasis<Dim-1>(basisCode, polyDegree - 1, fullTensorization);

			Poisson_HHO<Dim>* problem = new Poisson_HHO<Dim>(mesh, solution, sourceFunction, reconstructionBasis, cellBasis, faceBasis, staticCondensation, outputDirectory);

			cout << endl;
			cout << "----------------------- Assembly -------------------------" << endl;
			problem->Assemble(action);

			if ((action & Action::ExportFaces) == Action::ExportFaces)
				mesh->ExportFacesToMatlab(outputDirectory);

			if ((action & Action::SolveSystem) == Action::SolveSystem)
			{
				cout << endl;
				cout << "------------------- Linear system resolution ------------------" << endl;

				Solver* solver = CreateSolver(solverCode, problem, solverTolerance, staticCondensation, nMultigridLevels, wLoops, useGalerkinOperator, preSmootherCode, postSmootherCode, nPreSmoothingIterations, nPostSmoothingIterations, faceBasis->Size());
				cout << "Solver: " << *solver << endl << endl;
				solver->Setup(problem->A);
				problem->Solution = solver->Solve(problem->b);
				delete solver;

				if (staticCondensation && (action & Action::ExtractSolution) == Action::ExtractSolution)
					problem->ExtractTraceSystemSolution();

				problem->ReconstructHigherOrderApproximation();
				if ((action & Action::ExtractSolution) == Action::ExtractSolution)
					problem->ExtractSolution();
				double error = L2::Error<Dim>(mesh, reconstructionBasis, problem->ReconstructedSolution, exactSolution);
				cout << "L2 Error = " << error << endl;
			}

			delete problem;
			delete reconstructionBasis;
			delete cellBasis;
			delete faceBasis;
		}
		delete mesh;
		delete sourceFunction;
	}

private:
	Mesh<Dim>* BuildMesh(int n) { return nullptr;  }

	Solver* CreateSolver(string solverCode, Problem* problem, double tolerance, bool staticCondensation, 
		int nMultigridLevels, int wLoops, bool useGalerkinOperator, string preSmootherCode, string postSmootherCode, int nPreSmoothingIterations, int nPostSmoothingIterations, int blockSize)
	{
		Solver* solver = NULL;
		if (solverCode.compare("mg") == 0)
		{
			if (staticCondensation)
			{
				Poisson_HHO<Dim>* hhoProblem = dynamic_cast<Poisson_HHO<Dim>*>(problem);
				MultigridForHHO<Dim>* mg = new MultigridForHHO<Dim>(hhoProblem, nMultigridLevels);
				mg->WLoops = wLoops;
				mg->UseGalerkinOperator = useGalerkinOperator;
				mg->PreSmootherCode = preSmootherCode;
				mg->PostSmootherCode = postSmootherCode;
				mg->PreSmoothingIterations = nPreSmoothingIterations;
				mg->PostSmoothingIterations = nPostSmoothingIterations;
				mg->ComputeExactSolution = hhoProblem->_mesh->Elements.size() <= (Dim == 2 ? 32 * 32 : 8 * 8);
				solver = mg;
			}
			else
				assert(false && "Multigrid only applicable on HHO discretization with static condensation.");
		}
		else if (solverCode.compare("lu") == 0)
			solver = new EigenSparseLU();
		else if (solverCode.compare("cg") == 0)
			solver = new EigenCG(tolerance);
		else if (solverCode.compare("bgs") == 0)
			solver = new BlockGaussSeidel(blockSize);
		else if (solverCode.compare("bj") == 0)
			solver = new BlockJacobi(blockSize);
		else if (solverCode.compare("agmg") == 0)
			solver = new AGMG(tolerance);
		else
			assert(false && "Unknown solver or not applicable!");

		IterativeSolver* iterativeSolver = dynamic_cast<IterativeSolver*>(solver);
		if (iterativeSolver != nullptr)
			iterativeSolver->Tolerance = tolerance;

		return solver;
	}
};

template <>
Mesh<1>* ProgramDim<1>::BuildMesh(int n)
{
	return new CartesianGrid1D(n);
}

template <>
Mesh<2>* ProgramDim<2>::BuildMesh(int n)
{
	return new CartesianGrid2D(n, n);
	//mesh = new CartesianPolygonalMesh2D(n, n);
}

template <>
Mesh<3>* ProgramDim<3>::BuildMesh(int n)
{
	return new CartesianGrid3D(n, n, n);
}