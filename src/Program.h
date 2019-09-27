#pragma once
#include <cstdio>
#include <iostream>
#include <functional>
#include "DG/Poisson_DG.h"
#include "HHO/Poisson_HHO.h"
#include "Mesh/1D/CartesianGrid1D.h"
#include "Mesh/2D/CartesianGrid2D.h"
#include "Mesh/3D/CartesianGrid3D.h"
#include "Mesh/2D/CartesianPolygonalMesh2D.h"
#include "Utils/Action.h"
#include "Utils/Timer.h"
#include "Utils/DiffusionPartition.h"
#include "Solver/ConjugateGradient.h"
#include "Solver/MultigridForHHO.h"
#include "Solver/BlockJacobi.h"
#include "Solver/EigenCG.h"
#include "Solver/AGMG.h"
#include "Utils/L2.h"
using namespace std;

class Program
{
public:
	Program() {}
	virtual void Start(string rhsCode, double kappa1, double kappa2, double anisotropyRatio, string partition, BigNumber n, string discretization, string basisCode, int polyDegree, bool usePolynomialSpaceQ,
		int penalizationCoefficient, bool staticCondensation, Action action, 
		int nMultigridLevels, int matrixMaxSizeForCoarsestLevel, int wLoops, bool useGalerkinOperator, string preSmootherCode, string postSmootherCode, int nPreSmoothingIterations, int nPostSmoothingIterations,
		CoarseningStrategy coarseningStgy, string initialGuessCode, string outputDirectory, string solverCode, double solverTolerance) = 0;
};

template <int Dim>
class ProgramDim : public Program
{
public:
	ProgramDim() : Program() {}

	void Start(string rhsCode, double kappa1, double kappa2, double anisotropyRatio, string partition, BigNumber n, string discretization, string basisCode, int polyDegree, bool usePolynomialSpaceQ,
		int penalizationCoefficient, bool staticCondensation, Action action, 
		int nMultigridLevels, int matrixMaxSizeForCoarsestLevel, int wLoops, bool useGalerkinOperator, string preSmootherCode, string postSmootherCode, int nPreSmoothingIterations, int nPostSmoothingIterations,
		CoarseningStrategy coarseningStgy, string initialGuessCode, string outputDirectory, string solverCode, double solverTolerance)
	{
		Timer totalTimer;
		totalTimer.Start();

		//----------//
		//   Mesh   //
		//----------//

		Mesh<Dim>* mesh = BuildMesh(n);

		//mesh->CoarsenMesh(CoarseningStrategy::Agglomeration);
		//mesh = mesh->CoarseMesh;
		//cout << *mesh << endl << endl;
		//cout << "Coarse mesh" << endl << *(mesh->CoarseMesh) << endl << endl;

		//--------------------------------------------//
		//   Diffusion heterogeneity and anisotropy   //
		//--------------------------------------------//

		double rotationAngleInDegrees = 0;
		if (discretization.compare("dg") == 0)
		{
			rotationAngleInDegrees = 0;
			anisotropyRatio = 1;
		}
		double rotationAngleInRandians = rotationAngleInDegrees * M_PI / 180;
		DimVector<Dim> anisotropyCoefficients1 = kappa1 * DimVector<Dim>::Ones(Dim);
		DimVector<Dim> anisotropyCoefficients2 = kappa2 * DimVector<Dim>::Ones(Dim);
		if (Dim > 1)
		{
			anisotropyCoefficients1[0] = anisotropyRatio * anisotropyCoefficients1[1];
			anisotropyCoefficients2[0] = anisotropyRatio * anisotropyCoefficients2[1];
		}
		Tensor<Dim> diffTensor1(anisotropyCoefficients1, rotationAngleInRandians);
		Tensor<Dim> diffTensor2(anisotropyCoefficients2, rotationAngleInRandians);

		DiffusionPartition<Dim> diffusionPartition(partition, &diffTensor1, &diffTensor2);
		
		mesh->SetDiffusionCoefficient(&diffusionPartition);

		//---------------------------------------------//
		//   Analytical solution and source function   //
		//---------------------------------------------//

		function<double(DomPoint)> exactSolution = NULL;
		SourceFunction* sourceFunction;

		if (Dim == 1)
		{
			if (rhsCode.compare("sine") == 0)
			{
				if (diffusionPartition.IsHomogeneous)
				{
					exactSolution = [](DomPoint p)
					{
						double x = p.X;
						return sin(4 * M_PI * x) / (16 * pow(M_PI, 2));
					};
				}
				sourceFunction = new SourceFunction1D([](double x) { return sin(4 * M_PI * x); });
			}
			else if (rhsCode.compare("poly") == 0)
			{
				if (diffusionPartition.IsHomogeneous)
				{
					exactSolution = [](DomPoint p)
					{
						double x = p.X;
						return x * (1 - x);
					};
				}
				sourceFunction = new SourceFunction1D([](double x) { return 2; });
			}
			else if (rhsCode.compare("heterog") == 0)
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
			if (rhsCode.compare("sine") == 0)
			{
				if (diffusionPartition.IsHomogeneous && rotationAngleInDegrees == 0)
				{
					exactSolution = [anisotropyCoefficients1](DomPoint p)
					{
						double x = p.X;
						double y = p.Y;
						double a = anisotropyCoefficients1[0];
						double b = anisotropyCoefficients1[1];
						return 2 / (a + b) * sin(4 * M_PI * x)*sin(4 * M_PI * y);
					};
				}
				sourceFunction = new SourceFunction2D([](double x, double y) { return 2 * pow(4 * M_PI, 2) * sin(4 * M_PI * x)*sin(4 * M_PI * y); });
			}
			else if (rhsCode.compare("poly") == 0)
			{
				if (diffusionPartition.IsHomogeneous && diffusionPartition.IsIsotropic)
				{
					exactSolution = [](DomPoint p)
					{
						double x = p.X;
						double y = p.Y;
						return x * (1 - x) * y*(1 - y);
					};
				}
				sourceFunction = new SourceFunction2D([](double x, double y) { return 2 * (y*(1 - y) + x * (1 - x)); });
			}
		}
		else if (Dim == 3)
		{
			if (rhsCode.compare("sine") == 0)
			{
				if (diffusionPartition.IsHomogeneous && diffusionPartition.IsIsotropic)
				{
					exactSolution = [](DomPoint p)
					{
						double x = p.X;
						double y = p.Y;
						double z = p.Z;
						return sin(4 * M_PI * x)*sin(4 * M_PI * y)*sin(4 * M_PI * z);
					};
				}
				sourceFunction = new SourceFunction3D([](double x, double y, double z) {  return 3 * pow(4 * M_PI, 2) * sin(4 * M_PI * x)*sin(4 * M_PI * y)*sin(4 * M_PI * z); });
			}
			else if (rhsCode.compare("poly") == 0)
			{
				if (diffusionPartition.IsHomogeneous && diffusionPartition.IsIsotropic)
				{
					exactSolution = [](DomPoint p)
					{
						double x = p.X;
						double y = p.Y;
						double z = p.Z;
						return x * (1 - x)*y*(1 - y)*z*(1 - z);
					};
				}
				sourceFunction = new SourceFunction3D([](double x, double y, double z) { return 2 * ((y*(1 - y)*z*(1 - z) + x * (1 - x)*z*(1 - z) + x * (1 - x)*y*(1 - y))); });
			}
		}

		//--------------------------------//
		//   Discretization and solving   //
		//--------------------------------//

		if (discretization.compare("dg") == 0)
		{
			Poisson_DG<Dim> problem(rhsCode, sourceFunction, &diffusionPartition, outputDirectory);
			FunctionalBasis<Dim> basis(basisCode, polyDegree, usePolynomialSpaceQ);

			cout << endl;
			cout << "----------------------- Assembly -------------------------" << endl;
			Timer assemblyTimer;
			assemblyTimer.Start();

			problem.Assemble(mesh, &basis, penalizationCoefficient, action);
			
			assemblyTimer.Stop();
			cout << endl << "Assembly time: CPU = " << assemblyTimer.CPU() << ", elapsed = " << assemblyTimer.Elapsed() << endl;

			if ((action & Action::SolveSystem) == Action::SolveSystem)
			{
				cout << endl;
				cout << "------------------- Linear system resolution ------------------" << endl;

				Solver* solver = CreateSolver(solverCode, &problem, solverTolerance, staticCondensation, nMultigridLevels, matrixMaxSizeForCoarsestLevel, wLoops, useGalerkinOperator, preSmootherCode, postSmootherCode, nPreSmoothingIterations, nPostSmoothingIterations, coarseningStgy, basis.Size());
				problem.Solution = Solve(solver, problem.A, problem.b, initialGuessCode);

				if ((action & Action::ExtractSolution) == Action::ExtractSolution)
					problem.ExtractSolution();
				if ((action & Action::ComputeL2Error) == Action::ComputeL2Error && exactSolution != NULL)
				{
					double error = L2::Error<Dim>(mesh, basis, problem.Solution, exactSolution);
					cout << "L2 Error = " << std::scientific << error << endl;
				}
			}
		}
		else if (discretization.compare("hho") == 0)
		{
			FunctionalBasis<Dim> reconstructionBasis(basisCode, polyDegree, usePolynomialSpaceQ);
			FunctionalBasis<Dim> cellBasis(basisCode, polyDegree - 1, usePolynomialSpaceQ);
			FunctionalBasis<Dim-1> faceBasis(basisCode, polyDegree - 1, usePolynomialSpaceQ);

			Poisson_HHO<Dim> problem(mesh, rhsCode, sourceFunction, &reconstructionBasis, &cellBasis, &faceBasis, staticCondensation, &diffusionPartition, outputDirectory);

			cout << endl;
			cout << "----------------------- Assembly -------------------------" << endl;
			Timer assemblyTimer;
			assemblyTimer.Start();

			problem.Assemble(action);

			assemblyTimer.Stop();
			cout << endl << "Assembly time: CPU = " << assemblyTimer.CPU() << ", elapsed = " << assemblyTimer.Elapsed() << endl;

			if ((action & Action::ExportFaces) == Action::ExportFaces)
				mesh->ExportFacesToMatlab(outputDirectory, true);

			if ((action & Action::SolveSystem) == Action::SolveSystem)
			{
				cout << endl;
				cout << "------------------- Linear system resolution ------------------" << endl;

				Solver* solver = CreateSolver(solverCode, &problem, solverTolerance, staticCondensation, nMultigridLevels, matrixMaxSizeForCoarsestLevel, wLoops, useGalerkinOperator, preSmootherCode, postSmootherCode, nPreSmoothingIterations, nPostSmoothingIterations, coarseningStgy, faceBasis.Size());
				problem.Solution = Solve(solver, problem.A, problem.b, initialGuessCode);

				if (staticCondensation && (action & Action::ExtractSolution) == Action::ExtractSolution)
					problem.ExtractTraceSystemSolution();

				if ((action & Action::ExtractSolution) == Action::ExtractSolution || ((action & Action::ComputeL2Error) == Action::ComputeL2Error && exactSolution != NULL))
				{
					problem.ReconstructHigherOrderApproximation();
					if ((action & Action::ExtractSolution) == Action::ExtractSolution)
						problem.ExtractSolution();
					if ((action & Action::ComputeL2Error) == Action::ComputeL2Error && exactSolution != NULL)
					{
						double error = L2::Error<Dim>(mesh, reconstructionBasis, problem.ReconstructedSolution, exactSolution);
						cout << "L2 Error = " << std::scientific << error << endl;
					}
				}
			}
		}
		delete mesh;
		delete sourceFunction;

		totalTimer.Stop();
		cout << endl << "Total time: CPU = " << totalTimer.CPU() << ", elapsed = " << totalTimer.Elapsed() << endl;
	}

private:
	Mesh<Dim>* BuildMesh(int n) { return nullptr; }

	Solver* CreateSolver(string solverCode, Problem* problem, double tolerance, bool staticCondensation, 
		int nMultigridLevels, int matrixMaxSizeForCoarsestLevel, int wLoops, bool useGalerkinOperator, string preSmootherCode, string postSmootherCode, int nPreSmoothingIterations, int nPostSmoothingIterations, CoarseningStrategy coarseningStgy, int blockSize)
	{
		Solver* solver = NULL;
		if (solverCode.compare("mg") == 0 || solverCode.compare("pcgmg") == 0 || solverCode.compare("mg2") == 0 || solverCode.compare("pcgmg2") == 0)
		{
			if (staticCondensation)
			{
				Poisson_HHO<Dim>* hhoProblem = dynamic_cast<Poisson_HHO<Dim>*>(problem);
				int algoNumber = solverCode.compare("mg") == 0 || solverCode.compare("pcgmg") == 0 ? 1 : 2;
				MultigridForHHO<Dim>* mg = new MultigridForHHO<Dim>(hhoProblem, algoNumber, nMultigridLevels);
				mg->MatrixMaxSizeForCoarsestLevel = matrixMaxSizeForCoarsestLevel;
				mg->WLoops = wLoops;
				mg->UseGalerkinOperator = useGalerkinOperator;
				mg->PreSmootherCode = preSmootherCode;
				mg->PostSmootherCode = postSmootherCode;
				mg->PreSmoothingIterations = nPreSmoothingIterations;
				mg->PostSmoothingIterations = nPostSmoothingIterations;
				mg->CoarseningStgy = coarseningStgy;

				if (solverCode.compare("mg") == 0 || solverCode.compare("mg2") == 0)
					solver = mg;
				else if (solverCode.compare("pcgmg") == 0 || solverCode.compare("pcgmg2") == 0)
				{
					ConjugateGradient* cg = new ConjugateGradient();
					cg->Precond = Preconditioner(mg);
					solver = cg;
				}
			}
			else
				assert(false && "Multigrid only applicable on HHO discretization with static condensation.");
		}
		else if (solverCode.compare("lu") == 0)
			solver = new EigenSparseLU();
		else if (solverCode.compare("cg") == 0)
			solver = new ConjugateGradient();
		else if (solverCode.compare("eigencg") == 0)
			solver = new EigenCG();
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

	Eigen::VectorXd Solve(Solver* solver, const SparseMatrix& A, const Eigen::VectorXd& b, string initialGuessCode)
	{
		Timer solverTimer;
		solverTimer.Start();

		cout << "Solver: " << *solver << endl << endl;

		Eigen::VectorXd x;
		IterativeSolver* iterativeSolver = dynamic_cast<IterativeSolver*>(solver);
		if (iterativeSolver != nullptr)
		{
			iterativeSolver->ComputeExactSolution = A.rows() <= 2000;
			solver->Setup(A);
			x = iterativeSolver->Solve(b, initialGuessCode);
			cout << iterativeSolver->IterationCount << " iterations." << endl;
		}
		else
		{
			solver->Setup(A);
			x = solver->Solve(b);
		}

		delete solver;

		solverTimer.Stop();
		cout << "Solving time: CPU = " << solverTimer.CPU() << ", elapsed = " << solverTimer.Elapsed() << endl << endl;
		return x;
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