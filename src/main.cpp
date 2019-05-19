#include <cstdio>
#include <iostream>
#include <functional>
#include <getopt.h>
#include <Eigen/Core>
#include "DG/Poisson_DG.h"
#include "HHO/Poisson_HHO.h"
#include "FunctionalBasis/FunctionalBasis.h"
#include "Mesh/CartesianGrid1D.h"
#include "Mesh/CartesianGrid2D.h"
#include "Mesh/CartesianGrid3D.h"
#include "Mesh/CartesianShape.h"
#include "Utils/Action.h"
#include "Utils/DiffusionPartition.h"
#include "Solver/MultigridForHHO.h"
using namespace std;


void print_usage() {
	cout << "--------------------------------------------------------" << endl;
	cout << "Arguments:" << endl;
	cout << "-h:\t\t		help --> print usage" << endl;
	cout << "-d {1|2|3}:\t\t		space dimension (default: 1)" << endl;
	cout << "-k NUM:\t\t\t			diffusion coefficient k1 in the first part of the domain partition, while k2=1 in the second part (default: 1)" << endl;
	cout << "-s {sine|poly|hetero}:\t\t	analytical solution (default: sine)" << endl;
	cout << "\t\t\t\t\t					'sine'   = sine solution" << endl;
	cout << "\t\t\t\t\t					'poly'   = polynomial solution of global degree 2*d" << endl;
	cout << "\t\t\t\t\t					'hetero' = (1D only) heterogeneous diffusion-specific analytical solution" << endl;
	cout << "-n NUM:\t\t		number of subdivisions in each cartesian dimension (default: 5)" << endl;
	cout << "-t {dg|hho}:\t\t	discretization method (default: dg)" << endl;
	cout << "\t\t\t 'dg'  = Discontinuous Galerkin (Symmetric Interior Penalty)" << endl;
	cout << "\t\t\t 'hho' = Hybrid High Order" << endl;
	cout << "-b {monomials|legendre|bernstein}:	polynomial basis (default: monomials)" << endl;
	cout << "-p NUM:\t\t		max polynomial degree (default: 2)" << endl;
	cout << "-f:\t\t		full tensorization of the polynomials when d=2 or 3 (space Q) (default: false = space P)" << endl;
	cout << "-z NUM:\t\t\t	penalization coefficient (default: -1 = automatic)" << endl;
	cout << "-c:\t\t		static condensation (HHO only) (default: no static condensation)" << endl;
	cout << "-a {e|c|m|s}+\t\t	action (default: es): " << endl;
	cout <<	"\t\t\t 'e' = export system" << endl;
	cout << "\t\t\t 'c' = export all components of the matrix in separate files" << endl;
	cout << "\t\t\t 'm' = export mass matrix" << endl;
	cout << "\t\t\t 's' = solve system" << endl;
	cout << "\t\t\t 'v' = export solution vector (requires 's')" << endl;
	cout << "-o PATH:\t		output directory to export files (default: ./)" << endl;
	cout << "--------------------------------------------------------" << endl;
}

void argument_error(string msg)
{
	print_usage();
	cout << "Argument error: " << msg << endl;
	exit(EXIT_FAILURE);
}

int main(int argc, char* argv[])
{
	cout << "-------------------------- START ------------------------" << endl;
	cout << "Option -h for help." << endl;
	cout << "---------------------------------------------------------" << endl;
	Eigen::initParallel();

	int dimension = 1;
	string solution = "sine";
	double kappa1 = 1;
	double kappa2 = 1;
	BigNumber n = 5;
	string discretization = "dg";
	string basisCode = "monomials";
	int polyDegree = 2;
	bool fullTensorization = false;
	int penalizationCoefficient = -1;
	bool staticCondensation = false;
	string a = "es";
	string outputDirectory = ".";

	int option = 0;
	while ((option = getopt(argc, argv, "d:k:s:n:t:b:p:z:a:o:hfc")) != -1) 
	{
		switch (option) 
		{
			case 'h': print_usage(); exit(EXIT_SUCCESS);
				break;
			case 'd': dimension = atoi(optarg);
				if (dimension < 1 || dimension > 3)
					argument_error("dimension " + to_string(dimension) + "! Are you kidding?! Stop wasting my time.");
				break;
			case 's': solution = optarg;
				if (solution.compare("sine") != 0 && solution.compare("poly") != 0 && solution.compare("hetero") != 0)
					argument_error("unknown analytical solution '" + solution + "'. Check -s argument.");
				break;
			case 'k': kappa1 = atof(optarg);
				break;
			case 'n': n = stoul(optarg, nullptr, 0);
				break;
			case 't': discretization = optarg;
				if (discretization.compare("dg") != 0 && discretization.compare("hho") != 0)
					argument_error("unknown discretization '" + discretization + "'. Check -t argument.");
				break;
			case 'b': basisCode = optarg;
				if (basisCode.compare("monomials") != 0 && basisCode.compare("legendre") != 0 && basisCode.compare("bernstein") != 0 && basisCode.compare("hemker") != 0)
					argument_error("unknown polynomial basis '" + basisCode + "'. Check -b argument.");
				break;
			case 'p': polyDegree = atoi(optarg);
				break;
			case 'f': fullTensorization = true;
				break;
			case 'z': penalizationCoefficient = atoi(optarg);
				break;
			case 'c': staticCondensation = true;
				break;
			case 'a': a = optarg;
				break;
			case 'o': outputDirectory = optarg;
				break;
			default: print_usage();
				exit(EXIT_FAILURE);
		}
	}

	if (dimension != 1 && solution.compare("hetero") == 0)
		argument_error("-s hetero is only supported in 1D.");

	Action action = Action::LogAssembly;
	for (size_t i = 0; i < a.length(); i++)
	{
		if (a[i] == 'e')
			action |= Action::ExtractSystem;
		else if (a[i] == 'c')
			action |= Action::ExtractComponentMatrices;
		else if (a[i] == 'm')
			action |= Action::ExtractMassMatrix;
		else if (a[i] == 's')
			action |= Action::SolveSystem;
		else if (a[i] == 'v')
			action |= Action::ExtractSolution;
		else
			argument_error("unknown action '" + to_string(a[i]) + "'. Check -a argument.");
	}

	function<double(DomPoint)> exactSolution = NULL;
	SourceFunction* sourceFunction;

	function<bool(DomPoint)> isInPart1 = [](DomPoint p) { return p.X < 0.5; };
	DiffusionPartition diffusionPartition(isInPart1, kappa1, kappa2);

	//------------//
	//     1D     //
	//------------//

	if (dimension == 1)
	{
		Mesh<1>* mesh = new CartesianGrid1D(n);
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
				double b1 = (1 + 3 * alpha)/(2 * alpha*(1 + alpha));
				double b2 = -(alpha + 3) / (2 * (1 + alpha));
				if (diffusionPartition.IsInPart1(p))
					return 4 * a1 *pow(x, 2) + 2 * b1 * x;
				else
					return 4 * a2 * pow(x - 1, 2) + 2 * b2 * (x - 1);
			};
			sourceFunction = new SourceFunction1D([&diffusionPartition](double x) { return 4; });
		}

		if (discretization.compare("dg") == 0)
		{
			Poisson_DG<1>* problem = new Poisson_DG<1>(solution, sourceFunction, diffusionPartition, outputDirectory);
			FunctionalBasis<1>* basis = new FunctionalBasis<1>(basisCode, polyDegree);

			cout << "----------------------- Assembly -------------------------" << endl;
			problem->Assemble(mesh, basis, penalizationCoefficient, action);

			if ((action & Action::SolveSystem) == Action::SolveSystem)
			{
				cout << "------------------- Linear system resolution ------------------" << endl;
				problem->Solve();
				if ((action & Action::ExtractSolution) == Action::ExtractSolution)
					problem->ExtractSolution();
				double error = L2::Error<1>(mesh, basis, problem->Solution, exactSolution);
				cout << "L2 Error = " << error << endl;
			}

			delete problem;
			delete basis;
		}
		else
			argument_error("HHO in 1D not implemented.");

		delete mesh;
	}

	//------------//
	//     2D     //
	//------------//

	else if (dimension == 2)
	{
		Mesh<2>* mesh = new CartesianGrid2D(n);
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
				return x*(1 - x) * y*(1 - y);
			};
			sourceFunction = new SourceFunction2D([&diffusionPartition](double x, double y) { return diffusionPartition.Coefficient(x) * 2 * (y*(1 - y) + x*(1 - x)); });
		}

		if (discretization.compare("dg") == 0)
		{
			Poisson_DG<2>* problem = new Poisson_DG<2>(solution, sourceFunction, diffusionPartition, outputDirectory);
			FunctionalBasis<2>* basis = new FunctionalBasis<2>(basisCode, polyDegree, fullTensorization);

			cout << "----------------------- Assembly -------------------------" << endl;
			problem->Assemble(mesh, basis, penalizationCoefficient, action);

			if ((action & Action::SolveSystem) == Action::SolveSystem)
			{
				problem->Solve();
				if ((action & Action::ExtractSolution) == Action::ExtractSolution)
					problem->ExtractSolution();
				double error = L2::Error<2>(mesh, basis, problem->Solution, exactSolution);
				cout << "L2 Error = " << error << endl;
			}

			delete problem;
			delete basis;
		}
		else if (discretization.compare("hho") == 0)
		{
			FunctionalBasis<2>* reconstructionBasis = new FunctionalBasis<2>(basisCode, polyDegree, fullTensorization);
			FunctionalBasis<2>* cellBasis = new FunctionalBasis<2>(basisCode, polyDegree - 1, fullTensorization);
			FunctionalBasis<1>* faceBasis = new FunctionalBasis<1>(basisCode, polyDegree - 1, fullTensorization);

			Poisson_HHO<2>* problem = new Poisson_HHO<2>(mesh, solution, sourceFunction, reconstructionBasis, cellBasis, faceBasis, staticCondensation, outputDirectory);

			cout << "----------------------- Assembly -------------------------" << endl;
			problem->Assemble(action);

			/*for (auto f : mesh->Faces)
			{
				IntervalFace* face = dynamic_cast<IntervalFace*>(f);
				cout << *face << endl;
			}*/

			if ((action & Action::SolveSystem) == Action::SolveSystem)
			{
				cout << "------------------- Linear system resolution ------------------" << endl;
				if (staticCondensation && mesh->N > 2)
				{
					MultigridForHHO<2> solver(problem);
					cout << "Solver: " << solver << endl;
					//BlockGaussSeidel solver(problem->HHO.nLocalFaceUnknowns);
					solver.Setup(problem->A);
					solver.Tolerance = 1e-5;
					cout << endl;
					problem->Solution = solver.Solve(problem->b);

				}
				else
					problem->Solve();
				
				if ((action & Action::ExtractSolution) == Action::ExtractSolution)
					problem->ExtractTraceSystemSolution();

				problem->ReconstructSolution();
				if ((action & Action::ExtractSolution) == Action::ExtractSolution)
					problem->ExtractSolution();
				double error = L2::Error<2>(mesh, reconstructionBasis, problem->ReconstructedSolution, exactSolution);
				cout << "L2 Error = " << error << endl;
			}

			delete problem;
			delete reconstructionBasis;
			delete cellBasis;
			delete faceBasis;
		}
		delete mesh;
	}

	//------------//
	//     3D     //
	//------------//

	else if (dimension == 3)
	{
		Mesh<3>* mesh = new CartesianGrid3D(n);
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

		if (discretization.compare("dg") == 0)
		{
			cout << "----------------------- Assembly -------------------------" << endl;
			Poisson_DG<3>* problem = new Poisson_DG<3>(solution, sourceFunction, diffusionPartition, outputDirectory);
			FunctionalBasis<3>* basis = new FunctionalBasis<3>(basisCode, polyDegree, fullTensorization);

			problem->Assemble(mesh, basis, penalizationCoefficient, action);

			if ((action & Action::SolveSystem) == Action::SolveSystem)
			{
				cout << "------------------- Linear system resolution ------------------" << endl;
				problem->Solve();
				if ((action & Action::ExtractSolution) == Action::ExtractSolution)
					problem->ExtractSolution();
				double error = L2::Error<3>(mesh, basis, problem->Solution, exactSolution);
				cout << "L2 Error = " << error << endl;
			}

			delete problem;
			delete basis;
		}
		else if (discretization.compare("hho") == 0)
		{
			FunctionalBasis<3>* reconstructionBasis = new FunctionalBasis<3>(basisCode, polyDegree, fullTensorization);
			FunctionalBasis<3>* cellBasis = new FunctionalBasis<3>(basisCode, polyDegree - 1, fullTensorization);
			FunctionalBasis<2>* faceBasis = new FunctionalBasis<2>(basisCode, polyDegree - 1, fullTensorization);

			Poisson_HHO<3>* problem = new Poisson_HHO<3>(mesh, solution, sourceFunction, reconstructionBasis, cellBasis, faceBasis, staticCondensation, outputDirectory);

			cout << "----------------------- Assembly -------------------------" << endl;
			problem->Assemble(action);

			if ((action & Action::SolveSystem) == Action::SolveSystem)
			{
				cout << "------------------- Linear system resolution ------------------" << endl;
				problem->Solve();
				problem->ReconstructSolution();
				if ((action & Action::ExtractSolution) == Action::ExtractSolution)
					problem->ExtractSolution();
				double error = L2::Error<3>(mesh, reconstructionBasis, problem->ReconstructedSolution, exactSolution);
				cout << "L2 Error = " << error << endl;
			}

			delete problem;
			delete reconstructionBasis;
			delete cellBasis;
			delete faceBasis;
		}
		delete mesh;
	}

	delete sourceFunction;

	cout << "----------------- SUCCESSFUL TERMINATION ----------------" << endl;
    return EXIT_SUCCESS;
}
