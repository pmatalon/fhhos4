#include <cstdio>
#include <iostream>
#include <getopt.h>
#include <regex>
#include <Eigen/Core>
#include "Program.h"
using namespace std;


void print_usage() {
	cout << "-------------------------------------------------------------------------------------" << endl;
	cout << "Arguments:" << endl;
	cout << "-h                   : help --> print usage" << endl;
	cout << endl;
	cout << "------------------------- Problem definition -------------------------" << endl;
	cout << endl;
	cout << "-d {1|2|3}           : space dimension (default: 1)" << endl;
	cout << "-kappa NUM           : diffusion coefficient kappa1 in the first part of the domain partition," << endl;
	cout << "                       while kappa2=1 in the second part (default: 1 = homogeneous diffusion)" << endl;
	cout << "                       (DG only)" << endl;
	cout << "-s CODE              : analytical solution (default: 'sine')" << endl;
	cout << "                                 'sine'   = sine solution" << endl;
	cout << "                                 'poly'   = polynomial solution of global degree 2*d" << endl;
	cout << "                                 'hetero' = (DG 1D only) heterogeneous diffusion-specific analytical solution" << endl;
	cout << endl;
	cout << "------------------------- Discretization -------------------------" << endl;
	cout << endl;
	cout << "-n NUM               : number of subdivisions in each cartesian dimension (default: 5)" << endl;
	cout << "-t {dg|hho}          : discretization method (default: 'dg')" << endl;
	cout << "                                 'dg'  = Discontinuous Galerkin (Symmetric Interior Penalty)" << endl;
	cout << "                                 'hho' = Hybrid High Order" << endl;
	cout << "-b BASIS             : polynomial basis (default: monomials)" << endl;
	cout << "                                 'monomials'" << endl;
	cout << "                                 'legendre'" << endl;
	cout << "                                 'nlegendre' (normalized Legendre)" << endl;
	cout << "                                 'bernstein'" << endl;
	cout << "                                 'hemker'" << endl;
	cout << "-p NUM               : polynomial degree of approximation (default: 2)" << endl;
	cout << "-f                   : full tensorization of the polynomials when d=2 or 3 (i.e. space Q) (default: space P)" << endl;
	cout << "-z NUM               : penalization coefficient (default: -1 = automatic)" << endl;
	cout << "-c                   : static condensation (HHO only) (default: no static condensation)" << endl;
	cout << endl;
	cout << "------------------------- Linear solvers -------------------------" << endl;
	cout << endl;
	cout << "-v SOLVER            : linear solver (default: 'lu'): " << endl;
	cout << "                                 'lu'   = LU factorization (Eigen library)" << endl;
	cout << "                                 'cg'   = Conjugate gradient with Jacobi preconditioner (Eigen library)" << endl;
	cout << "                                 'bj'   = Block Jacobi: the block size is set to the number of DOFs per cell (DG) or face (HHO)" << endl;
	cout << "                                 'bgs'  = Block Gauss-Seidel: the block size is set to the number of DOFs per cell (DG) or face (HHO)" << endl;
	cout << "                                 'mg'   = Custom multigrid for HHO with static condensation" << endl;
	cout << "                                 'agmg' = Yvan Notay's AGMG solver" << endl;
	cout << "-cycle [V|W],NUM,NUM : multigrid cycle: V,1,0 or W,1,1 for example" << endl;
	cout << "-w NUM               : number of loops in the W-cycle" << endl;
	cout << "                                  1     - V-cycle (default)" << endl;
	cout << "                                 >1     - W-cycle" << endl;
	cout << "-l NUM               : number of multigrid levels (HHO with static cond. only):" << endl;
	cout << "                                  0     - (default) automatic coarsening until the matrix dimension reaches 100 or less" << endl;
	cout << "                                  other - fixed number of levels" << endl;
	cout << "-g {0|1}             : coarse grid operator for the multigrid" << endl;
	cout << "                                  0     - discretized operator" << endl;
	cout << "                                  1     - Galerkin operator (default)" << endl;
	cout << "-smoothers CODE,CODE : pre-smoother,post-smoother:" << endl;
	cout << "                                 'j'    = Jacobi" << endl;
	cout << "                                 'gs'   = Gauss-Seidel" << endl;
	cout << "                                 'rgs'  = Reverse Gauss-Seidel" << endl;
	cout << "                                 'bj'   = Block Jacobi: the block size is set to the number of DOFs per face" << endl;
	cout << "                                 'bgs'  = Block Gauss-Seidel: the block size is set to the number of DOFs per face" << endl;
	cout << "                                 'rbgs' = Reverse Block Gauss-Seidel" << endl;
	cout << endl;
	cout << "------------------------- Miscellaneous -------------------------" << endl;
	cout << endl;
	cout << "-threads NUM         : max number of threads used for parallelism" << endl;
	cout << "                                  0     - automatic (default)" << endl;
	cout << "                                  1     - sequential execution" << endl;
	cout << "                                  other - requested number of threads" << endl;
	cout << "-a {e|c|f|s|v}+      : action (default: 'es'): " << endl;
	cout << "                                 'e' = export system" << endl;
	cout << "                                 'c' = export all components of the matrix in separate files" << endl;
	cout << "                                 'f' = export faces for Matlab" << endl;
	cout << "                                 's' = solve system" << endl;
	cout << "                                 'v' = export solution vector (requires 's')" << endl;
	cout << "-o PATH              : output directory to export files (default: ./)" << endl;
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
	int nMultigridLevels = 0;
	int wLoops = 1;
	bool useGalerkinOperator = true;
	string preSmootherCode = "bgs";
	string postSmootherCode = "rbgs";
	int nPreSmoothingIterations = 1;
	int nPostSmoothingIterations = 1;
	string outputDirectory = ".";
	string solverCode = "lu";
	double solverTolerance = 1e-8;

	enum {
		OPT_Kappa = 1000,
		OPT_MGCycle,
		OPT_Smoothers,
		OPT_Threads
	};

	static struct option long_opts[] = {
		 { "help", no_argument, NULL, 'h' },
		 { "kappa", required_argument, NULL, OPT_Kappa },
		 { "cycle", required_argument, NULL, OPT_MGCycle },
		 { "smoothers", required_argument, NULL, OPT_Smoothers },
		 { "threads", required_argument, NULL, OPT_Threads },
		 { NULL, 0, NULL, 0 }
	};

	int long_index = 0;
	int option = 0;
	while ((option = getopt_long_only(argc, argv, "d:s:n:t:b:p:z:a:l:o:v:w:g:hfc", long_opts, &long_index)) != -1)
	{
		switch (option) 
		{
			case 'h': 
				print_usage(); 
				exit(EXIT_SUCCESS);
				break;
			case 'd': 
				dimension = atoi(optarg);
				if (dimension < 1 || dimension > 3)
					argument_error("dimension " + to_string(dimension) + "! Are you kidding?! Stop wasting my time.");
				break;
			case 's': 
				solution = optarg;
				if (solution.compare("sine") != 0 && solution.compare("poly") != 0 && solution.compare("hetero") != 0)
					argument_error("unknown analytical solution '" + solution + "'. Check -s argument.");
				break;
			case OPT_Kappa: 
				kappa1 = atof(optarg);
				break;
			case 'n': 
				n = stoul(optarg, nullptr, 0);
				break;
			case 't': 
				discretization = optarg;
				if (discretization.compare("dg") != 0 && discretization.compare("hho") != 0)
					argument_error("unknown discretization '" + discretization + "'. Check -t argument.");
				break;
			case 'b': 
				basisCode = optarg;
				if (basisCode.compare("monomials") != 0 && basisCode.compare("legendre") != 0 && basisCode.compare("nlegendre") != 0 && basisCode.compare("bernstein") != 0 && basisCode.compare("hemker") != 0)
					argument_error("unknown polynomial basis '" + basisCode + "'. Check -b argument.");
				break;
			case 'p': 
				polyDegree = atoi(optarg);
				break;
			case 'f': 
				fullTensorization = true;
				break;
			case 'z': 
				penalizationCoefficient = atoi(optarg);
				break;
			case 'c': 
				staticCondensation = true;
				break;
			case 'a': 
				a = optarg;
				break;
			case 'v': 
				solverCode = optarg;
				break;
			case OPT_MGCycle:
			{
				string s(optarg);
				regex pattern("^([vwVW]),([[:digit:]]+),([[:digit:]]+)$");
				smatch matches;

				if (std::regex_search(s, matches, pattern))
				{
					char cycleLetter = matches.str(1)[0];
					if (cycleLetter == 'v' || cycleLetter == 'V')
						wLoops = 1;
					else if (cycleLetter == 'w' || cycleLetter == 'W')
						wLoops = 2;
					nPreSmoothingIterations = stoi(matches.str(2));
					nPostSmoothingIterations = stoi(matches.str(3));
				}
				else
					argument_error("syntax error in the multigrid cycle. Check -cycle argument.");
				break;
			}
			case 'l': 
				nMultigridLevels = atoi(optarg);
				break;
			case 'w': 
				wLoops = atoi(optarg);
				if (wLoops < 1)
					argument_error("the number of loops in the W-cycle must be >= 1. Check -w argument.");
				break;
			case 'g': 
				useGalerkinOperator = atoi(optarg);
				break;
			case OPT_Smoothers:
			{
				string s(optarg);
				auto pos = s.find(",");
				preSmootherCode = s.substr(0, s.find(","));
				postSmootherCode = s.substr(pos + 1);
				break;
			}
			case OPT_Threads:
				BaseParallelLoop::SetDefaultNThreads(atoi(optarg));
				break;
			case 'o': 
				outputDirectory = optarg;
				break;
			default: print_usage();
				exit(EXIT_FAILURE);
		}
	}

	if (dimension != 1 && solution.compare("hetero") == 0)
		argument_error("-s hetero is only supported in 1D.");

	if (dimension == 1 && discretization.compare("hho") == 0 && polyDegree != 1)
		argument_error("HHO in 1D only exists for p = 1.");

	if (discretization.compare("hho") == 0 && polyDegree == 0)
		argument_error("HHO does not exist with p=0. Linear approximation at least (p >= 1).");

	if (solverCode.compare("mg") == 0 && discretization.compare("dg") == 0)
		argument_error("Multigrid only applicable on HHO discretization.");

	if (solverCode.compare("mg") == 0 && discretization.compare("hho") == 0 && !staticCondensation)
		argument_error("Multigrid only applicable if the static condensation is enabled.");

	Action action = Action::LogAssembly;
	for (size_t i = 0; i < a.length(); i++)
	{
		if (a[i] == 'e')
			action |= Action::ExtractSystem;
		else if (a[i] == 'c')
			action |= Action::ExtractComponentMatrices;
		else if (a[i] == 'f')
			action |= Action::ExportFaces;
		else if (a[i] == 's')
			action |= Action::SolveSystem;
		else if (a[i] == 'v')
			action |= Action::ExtractSolution;
		else
			argument_error("unknown action '" + to_string(a[i]) + "'. Check -a argument.");
	}

	Program* program = nullptr;
	if (dimension == 1)
		program = new ProgramDim<1>();
	else if (dimension == 2)
		program = new ProgramDim<2>();
	else if (dimension == 3)
		program = new ProgramDim<3>();

	program->Start(solution, kappa1, kappa2, n, discretization, basisCode, polyDegree, fullTensorization, penalizationCoefficient, staticCondensation, action, 
		nMultigridLevels, wLoops, useGalerkinOperator, preSmootherCode, postSmootherCode, nPreSmoothingIterations, nPostSmoothingIterations,
		outputDirectory, solverCode, solverTolerance);

	delete program;

	cout << "----------------- SUCCESSFUL TERMINATION ----------------" << endl;
    return EXIT_SUCCESS;
}
