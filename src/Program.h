#pragma once
#include "DG/Poisson_DG.h"
#include "HHO/Poisson_HHO.h"
#include "Mesh/1D/CartesianGrid1D.h"
#include "Mesh/2D/CartesianGrid2D.h"
#include "Mesh/3D/CartesianGrid3D.h"
#include "Mesh/2D/CartesianPolygonalMesh2D.h"
#include "Mesh/2D/TriangularMesh.h"
#include "Mesh/2D/QuadrilateralMesh.h"
#include "Utils/Action.h"
#include "Utils/Timer.h"
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
	virtual void Start(string rhsCode, double kappa1, double kappa2, double anisotropyRatio, string partition, 
		BigNumber n, string discretization, string meshCode, string stabilization, string basisCode, int polyDegree, bool usePolynomialSpaceQ,
		int penalizationCoefficient, bool staticCondensation, Action action, 
		int nMultigridLevels, int matrixMaxSizeForCoarsestLevel, int wLoops, bool useGalerkinOperator, string preSmootherCode, string postSmootherCode, int nPreSmoothingIterations, int nPostSmoothingIterations,
		CoarseningStrategy coarseningStgy, string initialGuessCode, string outputDirectory, string solverCode, double solverTolerance) = 0;
};

template <int Dim>
class ProgramDim : public Program
{
public:
	ProgramDim() : Program() {}

	void Start(string rhsCode, double kappa1, double kappa2, double anisotropyRatio, string partition, 
		BigNumber n, string discretization, string meshCode, string stabilization, string basisCode, int polyDegree, bool usePolynomialSpaceQ,
		int penalizationCoefficient, bool staticCondensation, Action action, 
		int nMultigridLevels, int matrixMaxSizeForCoarsestLevel, int wLoops, bool useGalerkinOperator, string preSmootherCode, string postSmootherCode, int nPreSmoothingIterations, int nPostSmoothingIterations,
		CoarseningStrategy coarseningStgy, string initialGuessCode, string outputDirectory, string solverCode, double solverTolerance)
	{
		Timer totalTimer;
		totalTimer.Start();

		GaussLegendre::Init();

		//----------//
		//   Mesh   //
		//----------//

		Mesh<Dim>* mesh = BuildMesh(n, meshCode);

		if (n <= 4)
			mesh->SanityCheck();
		if (n <= 2)
			cout << *mesh << endl << endl;

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

		DomFunction exactSolution = NULL;
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
			else if (rhsCode.compare("one") == 0)
			{
				if (diffusionPartition.IsHomogeneous && diffusionPartition.IsIsotropic)
				{
					exactSolution = [](DomPoint p) { return 1; };
				}
				sourceFunction = new SourceFunction2D([](double x, double y) { return 0; });
			}
			else if (rhsCode.compare("x") == 0)
			{
				if (diffusionPartition.IsHomogeneous && diffusionPartition.IsIsotropic)
				{
					exactSolution = [](DomPoint p) { return p.X; };
				}
				sourceFunction = new SourceFunction2D([](double x, double y) { return 0; });
			}
			else if (rhsCode.compare("kellogg") == 0)
			{
				if (diffusionPartition.IsIsotropic)
				{
					exactSolution = [](DomPoint p)
					{
						// Conversion from [0,1]x[0,1] to [-1,1]x[-1,1]
						double x = 2 * p.X - 1;
						double y = 2 * p.Y - 1;

						// Conversion to polar coordinates (r, t)
						double r = sqrt(x*x + y * y);
						double t = 0;
						if (x > 0 && y >= 0)
							t = atan(y / x);
						else if (x > 0 && y < 0)
							t = atan(y / x) + 2 * M_PI;
						else if (x < 0)
							t = atan(y / x) + M_PI;
						else if (x == 0 && y > 0)
							t = M_PI / 2;
						else if (x == 0 && y < 0)
							t = 3 * M_PI / 2;
						else
							assert(false);

						// Function in polar coordinates
						double eps = 0.1;
						double nu = M_PI / 4;
						double ksi = -14.9225565104455152;

						if (t >= 0 && t <= M_PI / 2)
							return pow(r, eps) * cos((M_PI / 2 - ksi)*eps) * cos((t - M_PI / 2 + nu)*eps);
						if (t >= M_PI / 2 && t <= M_PI)
							return pow(r, eps) * cos(nu*eps) * cos((t - M_PI + ksi)*eps);
						if (t >= M_PI && t <= 3 * M_PI / 2)
							return pow(r, eps) * cos(ksi*eps) * cos((t - M_PI - nu)*eps);
						if (t >= 3 * M_PI / 2 && t <= 2 * M_PI)
							return pow(r, eps) * cos((M_PI / 2 - nu)*eps) * cos((t - 3 * M_PI / 2 - ksi)*eps);

						assert(false);
					};
				}
				sourceFunction = new SourceFunction2D([](double x, double y) { return 0; });
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

		//-------------------------//
		//   Boundary conditions   //
		//-------------------------//

		function<BoundaryConditionType(DomPoint)> getBoundaryConditionType = [](DomPoint p)
		{
			return BoundaryConditionType::Dirichlet;
		};
		DomFunction dirichletBC = [exactSolution](DomPoint p)
		{
			if (exactSolution != nullptr)
				return exactSolution(p);
			return 0.0; 
		};
		DomFunction neumannBC = [](DomPoint p) { return 0; };

		BoundaryConditions bc(getBoundaryConditionType, dirichletBC, neumannBC);

		mesh->SetBoundaryConditions(&bc);

		//mesh->CoarsenMesh(CoarseningStrategy::Standard);
		//cout << *mesh << endl << endl;
		//cout << "Coarse mesh" << endl << *(mesh->CoarseMesh) << endl << endl;

		//--------------------------------//
		//   Discretization and solving   //
		//--------------------------------//
		
		if (discretization.compare("dg") == 0)
		{
			FunctionalBasis<Dim> basis(basisCode, polyDegree, usePolynomialSpaceQ);
			Poisson_DG<Dim> problem(mesh, rhsCode, sourceFunction, &diffusionPartition, outputDirectory, &basis, penalizationCoefficient);

			cout << endl;
			cout << "----------------------- Assembly -------------------------" << endl;
			Timer assemblyTimer;
			assemblyTimer.Start();

			problem.Assemble(action);
			
			assemblyTimer.Stop();
			cout << endl << "Assembly time: CPU = " << assemblyTimer.CPU() << ", elapsed = " << assemblyTimer.Elapsed() << endl;

			if ((action & Action::SolveSystem) == Action::SolveSystem)
			{
				cout << endl;
				cout << "------------------- Linear system resolution ------------------" << endl;

				Solver* solver = CreateSolver(solverCode, &problem, solverTolerance, staticCondensation, nMultigridLevels, matrixMaxSizeForCoarsestLevel, wLoops, useGalerkinOperator, preSmootherCode, postSmootherCode, nPreSmoothingIterations, nPostSmoothingIterations, coarseningStgy, basis.Size());
				problem.SystemSolution = Solve(solver, problem.A, problem.b, initialGuessCode);

				if ((action & Action::ExtractSolution) == Action::ExtractSolution)
					problem.ExtractSolution();
				if ((action & Action::ComputeL2Error) == Action::ComputeL2Error && exactSolution != NULL)
				{
					double error = L2::Error<Dim>(mesh, basis, problem.SystemSolution, exactSolution);
					cout << "L2 Error = " << std::scientific << error << endl;
				}
			}
		}
		else if (discretization.compare("hho") == 0)
		{
			FunctionalBasis<Dim> reconstructionBasis(basisCode, polyDegree, usePolynomialSpaceQ);
			FunctionalBasis<Dim> cellBasis(basisCode, polyDegree - 1, usePolynomialSpaceQ);
			FunctionalBasis<Dim-1> faceBasis(basisCode, polyDegree - 1, usePolynomialSpaceQ);

			HHOParameters<Dim> hho(mesh, stabilization, &reconstructionBasis, &cellBasis, &faceBasis);
			
			Poisson_HHO<Dim> problem(mesh, rhsCode, sourceFunction, &hho, staticCondensation, &diffusionPartition, &bc, outputDirectory);

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
				problem.SystemSolution = Solve(solver, problem.A, problem.b, initialGuessCode);

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
		GaussLegendre::Free();

		totalTimer.Stop();
		cout << endl << "Total time: CPU = " << totalTimer.CPU() << ", elapsed = " << totalTimer.Elapsed() << endl;
	}

private:
	Mesh<Dim>* BuildMesh(int n, string meshCode) { return nullptr; }

	Solver* CreateSolver(string solverCode, Problem<Dim>* problem, double tolerance, bool staticCondensation, 
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
Mesh<1>* ProgramDim<1>::BuildMesh(int n, string meshCode)
{
	return new CartesianGrid1D(n);
}

template <>
Mesh<2>* ProgramDim<2>::BuildMesh(int n, string meshCode)
{
	if (meshCode.compare("cart") == 0)
		return new CartesianGrid2D(n, n);
		//return new CartesianPolygonalMesh2D(n, n);
	else if(meshCode.compare("tri") == 0)
		return new TriangularMesh(n, n);
	else if (meshCode.compare("quad") == 0)
		return new QuadrilateralMesh(n, n);
	assert(false);
}

template <>
Mesh<3>* ProgramDim<3>::BuildMesh(int n, string meshCode)
{
	return new CartesianGrid3D(n, n, n);
}