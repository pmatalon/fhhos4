#pragma once
#include "ProgramArguments.h"
#include "DG/Diffusion_DG.h"
#include "HHO/Diffusion_HHO.h"
#include "Mesh/1D/UniformMesh1D.h"
#include "Mesh/2D/Square_CartesianMesh.h"
#include "Mesh/2D/Square_CartesianPolygonalMesh.h"
#include "Mesh/2D/Square_TriangularMesh.h"
#include "Mesh/2D/Square_QuadrilateralMesh.h"
#include "Mesh/2D/Square_QuadrilateralAsPolygonalMesh.h"
#include "Mesh/3D/Cube_CartesianMesh.h"
#include "Mesh/3D/Cube_CartesianTetrahedralMesh.h"
#ifdef GMSH_ENABLED
#include "Mesh/2D/Square_GMSHCartesianMesh.h"
#include "Mesh/2D/Square_GMSHTriangularMesh.h"
#include "Mesh/2D/Square4quadrants_GMSHUnstructTriangularMesh.h"
#include "Mesh/2D/Square_GMSHUnstructTriangularMesh.h"
#include "Mesh/2D/Square_GMSHQuadrilateralMesh.h"
#include "Mesh/3D/Cube_GMSHTetrahedralMesh.h"
#include "Mesh/3D/Cube_GMSHCartesianMesh.h"
#endif
#include "Utils/Action.h"
#include "Utils/Timer.h"
#include "Solver/ConjugateGradient.h"
#include "Solver/FlexibleConjugateGradient.h"
#include "Solver/MultigridForHHO.h"
#include "Solver/BlockJacobi.h"
#include "Solver/EigenCG.h"
#ifdef AGMG_ENABLED
#include "Solver/AGMG.h"
#endif
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

		cout << "----------------------- Mesh construction -------------------------" << endl;
		Mesh<Dim>* mesh = BuildMesh(args);

		if ((args.Actions & Action::UnitTests) == Action::UnitTests)
		{
			// Unit tests
			TriangleShape::Test();
			QuadrilateralShape::Test();
			PolygonalShape::Test();
			TetrahedronShape::Test();
			TriangleIn3DShape::Test();

			mesh->SanityCheck();
			if (args.Discretization.N <= 2)
				cout << *mesh << endl << endl;

			/*if (args.Discretization.MeshCode.compare("gmsh-tri") == 0)
			{
				GMSHMesh<Dim>* gmshMesh = dynamic_cast<GMSHMesh<Dim>*>(mesh);
				gmshMesh->RenumberLikeMe();
			}*/
		}

		//--------------------------------------------//
		//   Diffusion heterogeneity and anisotropy   //
		//--------------------------------------------//

		double rotationAngleInDegrees = 0;
		if (args.Discretization.Method.compare("dg") == 0)
		{
			rotationAngleInDegrees = 0;
			args.Problem.AnisotropyRatio = 1;
		}
		double rotationAngleInRandians = rotationAngleInDegrees * M_PI / 180;
		DimVector<Dim> anisotropyCoefficients1 = args.Problem.Kappa1 * DimVector<Dim>::Ones(Dim);
		DimVector<Dim> anisotropyCoefficients2 = args.Problem.Kappa2 * DimVector<Dim>::Ones(Dim);
		if (Dim > 1)
		{
			anisotropyCoefficients1[0] = args.Problem.AnisotropyRatio * anisotropyCoefficients1[1];
			anisotropyCoefficients2[0] = args.Problem.AnisotropyRatio * anisotropyCoefficients2[1];
		}
		Tensor<Dim> diffTensor1(anisotropyCoefficients1, rotationAngleInRandians);
		Tensor<Dim> diffTensor2(anisotropyCoefficients2, rotationAngleInRandians);

		DiffusionPartition<Dim> diffusionPartition(args.Problem.Partition, &diffTensor1, &diffTensor2);

		mesh->SetDiffusionCoefficient(&diffusionPartition);
		
		//---------------------------------------------//
		//   Analytical solution and source function   //
		//---------------------------------------------//

		DomFunction exactSolution = NULL;
		SourceFunction* sourceFunction;

		if (Dim == 1)
		{
			if (args.Problem.RHSCode.compare("sine") == 0)
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
			else if (args.Problem.RHSCode.compare("poly") == 0)
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
			else if (args.Problem.RHSCode.compare("heterog") == 0)
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
			else
				Utils::FatalError("'" + args.Problem.RHSCode + "' is unknown or not implemented in 1D. Check -rhs argument.");
		}
		else if (Dim == 2)
		{
			if (args.Problem.RHSCode.compare("sine") == 0)
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
			else if (args.Problem.RHSCode.compare("poly") == 0)
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
			else if (args.Problem.RHSCode.compare("exp") == 0)
			{
				if (diffusionPartition.IsHomogeneous && diffusionPartition.IsIsotropic)
				{
					exactSolution = [](DomPoint p)
					{
						double x = p.X;
						double y = p.Y;
						return exp(x*y*y);
					};
				}
				sourceFunction = new SourceFunction2D([](double x, double y) { return (-pow(y,4) - 2*x*(1+2*x*y*y))*exp(x*y*y); });
			}
			else if (args.Problem.RHSCode.compare("one") == 0)
			{
				if (diffusionPartition.IsHomogeneous && diffusionPartition.IsIsotropic)
				{
					exactSolution = [](DomPoint p) { return 1; };
				}
				sourceFunction = new SourceFunction2D([](double x, double y) { return 0; });
			}
			else if (args.Problem.RHSCode.compare("zero") == 0)
			{
				if (diffusionPartition.IsHomogeneous && diffusionPartition.IsIsotropic)
				{
					exactSolution = [](DomPoint p) { return 0; };
				}
				sourceFunction = new SourceFunction2D([](double x, double y) { return 0; });
			}
			else if (args.Problem.RHSCode.compare("x") == 0)
			{
				if (diffusionPartition.IsHomogeneous && diffusionPartition.IsIsotropic)
				{
					exactSolution = [](DomPoint p) { return p.X; };
				}
				sourceFunction = new SourceFunction2D([](double x, double y) { return 0; });
			}
			else if (args.Problem.RHSCode.compare("kellogg") == 0)
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
			else
				Utils::FatalError("'" + args.Problem.RHSCode + "' is unknown or not implemented in 2D. Check -rhs argument.");
		}
		else if (Dim == 3)
		{
			if (args.Problem.RHSCode.compare("sine") == 0)
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
			else if (args.Problem.RHSCode.compare("poly") == 0)
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
			else if (args.Problem.RHSCode.compare("exp") == 0)
			{
				if (diffusionPartition.IsHomogeneous && diffusionPartition.IsIsotropic)
				{
					exactSolution = [](DomPoint p)
					{
						double x = p.X;
						double y = p.Y;
						double z = p.Z;
						return exp(x*y*y*z*z*z);
					};
				}
				sourceFunction = new SourceFunction3D([](double x, double y, double z) { return -(pow(y, 4)*pow(z, 6) + 2 * x*pow(z, 3) + 4 * x*x*y*y*pow(z, 6) + 6 * x*y*y*z + 9 * x*x*pow(y, 4)*pow(z, 4))*exp(x*y*y*z*z*z); });
			}
			else
				Utils::FatalError("'" + args.Problem.RHSCode + "' is unknown or not implemented in 3D. Check -rhs argument.");
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

		if ((args.Actions & Action::UnitTests) == Action::UnitTests)
		{
			if (args.Solver.SolverCode.compare("mg") == 0 || args.Solver.SolverCode.compare("cgmg") == 0 || args.Solver.SolverCode.compare("fcgmg") == 0)
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
		}

		//----------------------//
		//       Assembly       //
		//----------------------//

		Problem<Dim>* problem = nullptr;
		
		if (args.Discretization.Method.compare("dg") == 0)
		{
			FunctionalBasis<Dim>* basis = new FunctionalBasis<Dim>(args.Discretization.BasisCode, args.Discretization.PolyDegree, args.Discretization.UsePolynomialSpaceQ);
			problem = new Diffusion_DG<Dim>(mesh, args.Problem.RHSCode, sourceFunction, &diffusionPartition, args.OutputDirectory, basis, args.Discretization.PenalizationCoefficient);
		}
		else if (args.Discretization.Method.compare("hho") == 0)
		{
			FunctionalBasis<Dim>* reconstructionBasis = new FunctionalBasis<Dim>(args.Discretization.BasisCode, args.Discretization.PolyDegree, args.Discretization.UsePolynomialSpaceQ);
			FunctionalBasis<Dim>* cellBasis = new FunctionalBasis<Dim>(args.Discretization.BasisCode, args.Discretization.PolyDegree - 1, args.Discretization.UsePolynomialSpaceQ);
			FunctionalBasis<Dim - 1>* faceBasis = new FunctionalBasis<Dim - 1>(args.Discretization.BasisCode, args.Discretization.PolyDegree - 1, args.Discretization.UsePolynomialSpaceQ);

			HHOParameters<Dim>* hho = new HHOParameters<Dim>(mesh, args.Discretization.Stabilization, reconstructionBasis, cellBasis, faceBasis);

			problem = new Diffusion_HHO<Dim>(mesh, args.Problem.RHSCode, sourceFunction, hho, args.Discretization.StaticCondensation, &diffusionPartition, &bc, args.OutputDirectory);
		}
		else
			Utils::FatalError("Unknown discretization.");

		cout << endl;
		cout << "----------------------- Assembly -------------------------" << endl;
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
			cout << "------------------- Linear system solution ------------------" << endl;

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

				if ((args.Actions & Action::ExtractSolution) == Action::ExtractSolution || ((args.Actions & Action::ComputeL2Error) == Action::ComputeL2Error && exactSolution != NULL))
				{
					hhoPb->ReconstructHigherOrderApproximation();
					if ((args.Actions & Action::ExtractSolution) == Action::ExtractSolution)
					{
						hhoPb->ExtractHybridSolution();
						hhoPb->ExtractSolution();
					}
				}

				if ((args.Actions & Action::ExportFaces) == Action::ExportFaces)
					mesh->ExportFacesToMatlab(args.OutputDirectory, true);
			}

			//----------------------//
			//       L2 error       //
			//----------------------//

			if ((args.Actions & Action::ComputeL2Error) == Action::ComputeL2Error && exactSolution != NULL)
			{
				double error = problem->L2Error(exactSolution);
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
		delete sourceFunction;
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
Mesh<1>* ProgramDim<1>::BuildMesh(ProgramArguments& args)
{
	return new UniformMesh1D(args.Discretization.N);
}

template <>
Mesh<2>* ProgramDim<2>::BuildMesh(ProgramArguments& args)
{
	BigNumber nx = args.Discretization.N;
	BigNumber ny = args.Discretization.Ny == -1 ? args.Discretization.N : args.Discretization.Ny;
	string meshCode = args.Discretization.MeshCode;
	double stretch = args.Discretization.Stretch;

	CoarseningStrategy refinementStgy = args.Solver.MG.CoarseningStgy;
	if (args.Solver.MG.CoarseningStgy != CoarseningStrategy::BeyRefinement && args.Solver.MG.CoarseningStgy != CoarseningStrategy::SplittingRefinement)
		refinementStgy = CoarseningStrategy::SplittingRefinement;

	Mesh<2>* fineMesh = nullptr;
	if (meshCode.compare("cart") == 0)
		fineMesh = new Square_CartesianMesh(nx, ny);
	else if (meshCode.compare("cart-poly") == 0)
		fineMesh = new Square_CartesianPolygonalMesh(nx, ny);
	else if (meshCode.compare("tri") == 0)
		fineMesh = new Square_TriangularMesh(nx, ny);
	else if (meshCode.compare("quad") == 0)
		fineMesh = new Square_QuadrilateralMesh(nx, ny, stretch);
	else if (meshCode.compare("quad-poly") == 0)
		fineMesh = new Square_QuadrilateralAsPolygonalMesh(nx, ny, stretch);
#ifdef GMSH_ENABLED
	else if (meshCode.compare("gmsh-cart") == 0)
	{
		Mesh<2>* coarseMesh = new Square_GMSHCartesianMesh();
		fineMesh = coarseMesh->RefineUntilNElements(nx*ny, refinementStgy);
	}
	else if (meshCode.compare("gmsh-tri") == 0)
	{
		Mesh<2>* coarseMesh = new Square_GMSHTriangularMesh();
		fineMesh = coarseMesh->RefineUntilNElements(2*nx*ny, refinementStgy);
	}
	else if (meshCode.compare("4quad-gmsh-uns-tri") == 0)
	{
		Mesh<2>* coarseMesh = new Square4quadrants_GMSHUnstructTriangularMesh();
		fineMesh = coarseMesh->RefineUntilNElements(2*nx*ny, refinementStgy);
	}
	else if (meshCode.compare("gmsh-uns-tri") == 0)
	{
		if (args.Solver.MG.CoarseningStgy == CoarseningStrategy::SplittingRefinement)
		{
			Mesh<2>* coarseMesh = new Square_GMSHUnstructTriangularMesh(8);
			fineMesh = coarseMesh->RefineUntilNElements(2 * nx*ny, refinementStgy);
		}
		else
			fineMesh = new Square_GMSHUnstructTriangularMesh(args.Discretization.N);
	}
	else if (meshCode.compare("gmsh-quad") == 0)
	{
		Mesh<2>* coarseMesh = new Square_GMSHQuadrilateralMesh();
		fineMesh = coarseMesh->RefineUntilNElements(nx*ny, refinementStgy);
	}
	else if (meshCode.compare("gmsh") == 0)
	{
		Mesh<2>* coarseMesh = new GMSHMesh<2>(args.Discretization.MeshFilePath);
		fineMesh = coarseMesh->RefineUntilNElements(2*nx*ny, refinementStgy);
	}
#endif // GMSH_ENABLED
	else
		assert(false);
	
	if (args.Solver.MG.CoarseningStgy != CoarseningStrategy::BeyRefinement && args.Solver.MG.CoarseningStgy != CoarseningStrategy::SplittingRefinement && fineMesh->CoarseMesh)
		fineMesh->DeleteCoarseMeshes();

	return fineMesh;
}

template <>
Mesh<3>* ProgramDim<3>::BuildMesh(ProgramArguments& args)
{
	BigNumber n = args.Discretization.N;
	BigNumber nx = args.Discretization.N;
	BigNumber ny = args.Discretization.Ny == -1 ? args.Discretization.N : args.Discretization.Ny;
	BigNumber nz = args.Discretization.Nz == -1 ? args.Discretization.N : args.Discretization.Nz;
	string meshCode = args.Discretization.MeshCode;
	CoarseningStrategy refinementStgy = args.Solver.MG.CoarseningStgy;

	if (meshCode.compare("cart") == 0)
	{
		Cube_CartesianMesh* fineMesh = new Cube_CartesianMesh(nx, ny, nz);

		assert(fineMesh->Elements.size() == nx*ny*nz);
		if (nx == ny && ny == nz)
			assert(fineMesh->Faces.size() == 3*n*n*(n+1));

		return fineMesh;
	}
	else if (meshCode.compare("tetra") == 0)
	{
		Mesh<3>* fineMesh;
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

		return fineMesh;
	}
#ifdef GMSH_ENABLED
	else if (meshCode.compare("gmsh-cart") == 0)
	{
		if (nx != ny || nx != nz)
			Utils::FatalError("-ny, -ny not managed with this mesh");

		Mesh<3>* coarseMesh = new Cube_GMSHCartesianMesh();
		Mesh<3>* fineMesh = coarseMesh->RefineUntilNElements(n*n*n, refinementStgy);

		assert(fineMesh->Elements.size() == n*n*n);
		assert(fineMesh->Faces.size() == 3*n*n*(n+1));

		return fineMesh;
	}
	else if (meshCode.compare("gmsh-tetra") == 0)
	{
		if (nx != ny || nx != nz)
			Utils::FatalError("-ny, -ny not managed with this mesh");

		Mesh<3>* coarseMesh = new Cube_GMSHTetrahedralMesh();
		Mesh<3>* fineMesh = coarseMesh->RefineUntilNElements(6*n*n*n, refinementStgy);

		assert(fineMesh->Elements.size() == 6*n*n*n);
		assert(fineMesh->Faces.size() == 12*n*n*n + 6*n*n);

		return fineMesh;
	}
	else if (meshCode.compare("gmsh") == 0)
	{
		if (nx != ny || nx != nz)
			Utils::FatalError("-ny, -ny not managed with this mesh");

		Mesh<3>* coarseMesh;
		if (refinementStgy == CoarseningStrategy::BeyRefinement)
			coarseMesh = new GMSHTetrahedralMesh(args.Discretization.MeshFilePath);
		else
			coarseMesh = new GMSHMesh<3>(args.Discretization.MeshFilePath);
		Mesh<3>* fineMesh = coarseMesh->RefineUntilNElements(6*n*n*n, refinementStgy);
		return fineMesh;
	}
#endif // GMSH_ENABLED
	assert(false);
}
