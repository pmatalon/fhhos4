#include <getopt.h>
#include <regex>
#include "Program.h"
using namespace std;


void print_usage() {
	cout << "Arguments:" << endl;
	cout << endl;
	cout << "----------------------------------------------------------------------" << endl;
	cout << "                          Problem definition                          " << endl;
	cout << "----------------------------------------------------------------------" << endl;
	cout << endl;
	cout << "-d NUM" << endl;
	cout << "      Space dimension of the domain: 1, 2 or 3 (default: 2)." << endl;
	cout << endl;
	cout << "-heterog NUM" << endl;
	cout << "      Heterogeneity ratio. Constant diffusion coefficient in one part of the domain partition" << endl;
	cout << "      while equals to 1 in the second part." << endl;
	cout << "      Ex:       1             - homogeneous diffusion (default)" << endl;
	cout << "                0 < NUM < 1   - to be used with -rhs heterog" << endl;
	cout << "                1e4           - allows to set the order of magnitude of the ratio" << endl;
	cout << endl;
	cout << "-aniso NUM" << endl;
	cout << "      Anisotropy ratio valid across the whole domain." << endl;
	cout << "      Example: 1e4" << endl;
	cout << endl;
	cout << "-partition CODE" << endl;
	cout << "      Domain partition describing the heterogeneity pattern in case of heterogeneous diffusion." << endl;
	cout << "               halves     - the domain is split in two vertical halves" << endl;
	cout << "               chiasmus   - chiasmus shape (default)" << endl;
	cout << endl;
	cout << "-rhs CODE" << endl;
	cout << "      Right-hand side code determining the source function (default: sine)." << endl;
	cout << "      It also determines the analytical solution in the homogeneous isotropic case." << endl;
	cout << "               sine    - the source function and the analytical solution are a sine functions" << endl;
	cout << "               poly    - the source function is constant, the analytical solution is a polynomial of total degree 2*d" << endl;
	cout << "               one     - the source function is 0, the analytical solution is 1" << endl;
	cout << "               x       - the source function is 0, the analytical solution is x" << endl;
	cout << "               heterog - (1D only) heterogeneous diffusion-specific analytical solution" << endl;
	cout << "               kellogg - (2D only) heterogeneous diffusion-specific analytical solution (known benchmark)" << endl;
	cout << endl;
	cout << "----------------------------------------------------------------------" << endl;
	cout << "                             Discretization                           " << endl;
	cout << "----------------------------------------------------------------------" << endl;
	cout << endl;
	cout << "-n NUM" << endl;
	cout << "      Number of subdivisions in each cartesian dimension (default: 16)." << endl;
	cout << endl;
	cout << "-mesh CODE" << endl;
	cout << "      Type of mesh:" << endl;
	cout << "               cart   - Uniform Cartesian mesh (default)" << endl;
	cout << "               tri    - Trianglular mesh" << endl;
	cout << "               quad   - Quadrilateral mesh" << endl;
	cout << "-discr CODE" << endl;
	cout << "      Discretization method (default: hho)." << endl;
	cout << "               dg     - Discontinuous Galerkin (Symmetric Interior Penalty)" << endl;
	cout << "               hho    - Hybrid High Order" << endl;
	cout << endl;
	cout << "-stab CODE" << endl;
	cout << "      Stabilization term (only used in HHO)." << endl;
	cout << "               hho    - HHO stabilization term" << endl;
	cout << "               hdg    - HDG stabilization term" << endl;
	cout << endl;
	cout << "-b CODE" << endl;
	cout << "      Polynomial basis (default: legendre)." << endl;
	cout << "               monomials" << endl;
	cout << "               legendre" << endl;
	cout << "               nlegendre (normalized Legendre)" << endl;
	cout << "               bernstein" << endl;
	cout << "               hemker" << endl;
	cout << endl;
	cout << "-p NUM" << endl;
	cout << "      Polynomial degree of approximation (default: 1). In HHO, k = p-1." << endl;
	cout << endl;
	cout << "-poly-space CODE" << endl;
	cout << "      Polynomial space." << endl;
	cout << "               p - space P (default)" << endl;
	cout << "               q - space Q of tensor polynomials (when d=2 or 3)" << endl;
	cout << endl;
	cout << "-pen NUM" << endl;
	cout << "      Penalization coefficient in DG (default: -1 = automatic)." << endl;
	cout << endl;
	cout << "-no-static-cond" << endl;
	cout << "      Disables the static condensation in HHO." << endl;
	cout << endl;
	cout << "----------------------------------------------------------------------" << endl;
	cout << "                             Linear solver                            " << endl;
	cout << "----------------------------------------------------------------------" << endl;
	cout << endl;
	cout << "-s SOLVER" << endl;
	cout << "      Linear solver (default: lu)." << endl;
	cout << "              lu       - LU factorization (Eigen library)" << endl;
	cout << "              cg       - Conjugate Gradient, no preconditioner" << endl;
	cout << "              eigencg  - Conjugate Gradient (Eigen library) with diagonal preconditioner" << endl;
	cout << "              bj       - Block Jacobi: the block size is set to the number of DOFs per cell (DG) or face (HHO)" << endl;
	cout << "              bgs      - Block Gauss-Seidel: the block size is set to the number of DOFs per cell (DG) or face (HHO)" << endl;
	cout << "              mg       - Custom multigrid for HHO" << endl;
	cout << "              mg2      - Custom multigrid for HHO (identity for faces present on both grids)" << endl;
	cout << "              pcgmg    - Conjugate Gradient, preconditioned with the custom multigrid for HHO 'mg'" << endl;
	cout << "              pcgmg2   - Conjugate Gradient, preconditioned with the custom multigrid for HHO 'mg2'" << endl;
	cout << "              agmg     - Yvan Notay's AGMG solver" << endl;
	cout << endl;
	cout << "-cycle [V|W],NUM,NUM" << endl;
	cout << "      Multigrid cycle." << endl;
	cout << "      Examples: \"V,1,0\" or \"W,1,1\"." << endl;
	cout << endl;
	cout << "-w NUM" << endl;
	cout << "      Number of loops in the multigrid cycle." << endl;
	cout << "               1     - V-cycle (default)" << endl;
	cout << "              >1     - W-cycle" << endl;
	cout << endl;
	cout << "-l NUM" << endl;
	cout << "      Number of multigrid levels." << endl;
	cout << "              0     - (default) automatic coarsening until the matrix dimension reaches 100 or less" << endl;
	cout << "              other - fixed number of levels" << endl;
	cout << endl;
	cout << "-coarse-size NUM" << endl;
	cout << "      Matrix size limit below which the automatic coarsening stops and a direct solver is used (default: 1000)." << endl;
	cout << endl;
	cout << "-g {0|1}" << endl;
	cout << "      Coarse grid operator for the multigrid." << endl;
	cout << "              0     - discretized operator (default)" << endl;
	cout << "              1     - Galerkin operator" << endl;
	cout << endl;
	cout << "-smoothers CODE,CODE" << endl;
	cout << "      Pre-smoother,post-smoother: \"gs,rgs\" for example." << endl;
	cout << "              j    - Jacobi" << endl;
	cout << "              gs   - Gauss-Seidel" << endl;
	cout << "              rgs  - Reverse Gauss-Seidel" << endl;
	cout << "              bj   - Block Jacobi: the block size is set to the number of DOFs per face" << endl;
	cout << "              bgs  - Block Gauss-Seidel: the block size is set to the number of DOFs per face" << endl;
	cout << "              rbgs - Reverse Block Gauss-Seidel" << endl;
	cout << endl;
	cout << "-cs CODE" << endl;
	cout << "      Coarsening strategy of the multigrid. Default: standard coarsening." << endl;
	cout << "              s    - standard coarsening (merge colinear faces on the coarse mesh)" << endl;
	cout << "              a    - agglomeration coarsening (keep fine faces on the coarse mesh)" << endl;
	cout << endl;
	cout << "-initial-guess CODE" << endl;
	cout << "      Initial guess for the iterative solvers." << endl;
	cout << "              0     - zero vector (default)" << endl;
	cout << "              1     - all ones vector" << endl;
	cout << endl;
	cout << "----------------------------------------------------------------------" << endl;
	cout << "                             Miscellaneous                            " << endl;
	cout << "----------------------------------------------------------------------" << endl;
	cout << endl;
	cout << "-h, -help" << endl;
	cout << "      Prints usage." << endl;
	cout << endl;
	cout << "-threads NUM" << endl;
	cout << "      Max number of threads used for parallelism." << endl;
	cout << "              0     - automatic (default)" << endl;
	cout << "              1     - sequential execution" << endl;
	cout << "              other - requested number of threads" << endl;
	cout << endl;
	cout << "-a {CODE}+" << endl;
	cout << "      Action (default: sr)." << endl;
	cout << "              e   - export system" << endl;
	cout << "              c   - export all components of the matrix in separate files" << endl;
	cout << "              f   - export faces for Matlab" << endl;
	cout << "              s   - solve system" << endl;
	cout << "              v   - export solution vector (requires 's')" << endl;
	cout << "              r   - compute L2 error against the analytical solution if known" << endl;
	cout << endl;
	cout << "-o PATH" << endl;
	cout << "      Output directory to export files (default: ./)." << endl;
	cout << endl;
	cout << "--------------------------------------------------------" << endl;
}

void argument_error(string msg)
{
	cout << "Argument error: " << msg << endl;
	cout << "------------------------- FAILURE -------------------------" << endl;
	exit(EXIT_FAILURE);
}

int main(int argc, char* argv[])
{
	cout << "-------------------------- START --------------------------" << endl;
	cout << "Option -h for help." << endl;
	cout << "-----------------------------------------------------------" << endl;
	Eigen::initParallel();

	int dimension = 2;
	string rhsCode = "sine";
	double kappa1 = 1;
	double kappa2 = 1;
	double anisotropyRatio = 1;
	string partition = "chiasmus";
	BigNumber n = 16;
	string meshCode = "cart";
	string discretization = "hho";
	string stabilization = "hho";
	string basisCode = "legendre";
	int polyDegree = 1;
	bool usePolynomialSpaceQ = false;
	int penalizationCoefficient = -1;
	bool staticCondensation = true;
	string a = "sr";
	int nMultigridLevels = 0;
	int matrixMaxSizeForCoarsestLevel = 1000;
	int wLoops = 1;
	bool useGalerkinOperator = false;
	string preSmootherCode = "bgs";
	string postSmootherCode = "rbgs";
	int nPreSmoothingIterations = 1;
	int nPostSmoothingIterations = 1;
	CoarseningStrategy coarseningStgy = CoarseningStrategy::Standard;
	string outputDirectory = ".";
	string solverCode = "lu";
	double solverTolerance = 1e-8;
	string initialGuessCode = "0";

	enum {
		OPT_RightHandSide = 1000,
		OPT_Heterogeneity,
		OPT_Anisotropy,
		OPT_Partition,
		OPT_Discretization,
		OPT_Mesh,
		OPT_Stabilization,
		OPT_NoStaticCondensation,
		OPT_Penalization,
		OPT_PolySpace,
		OPT_MGCycle,
		OPT_CoarseMatrixSize,
		OPT_Smoothers,
		OPT_CoarseningStrategy,
		OPT_InitialGuess,
		OPT_Threads
	};

	static struct option long_opts[] = {
		 { "help", no_argument, NULL, 'h' },
		 { "rhs", required_argument, NULL, OPT_RightHandSide },
		 { "heterog", required_argument, NULL, OPT_Heterogeneity },
		 { "aniso", required_argument, NULL, OPT_Anisotropy },
		 { "partition", required_argument, NULL, OPT_Partition },
		 { "discr", required_argument, NULL, OPT_Discretization },
		 { "mesh", required_argument, NULL, OPT_Mesh },
		 { "stab", required_argument, NULL, OPT_Stabilization },
		 { "no-static-cond", required_argument, NULL, OPT_NoStaticCondensation },
		 { "pen", required_argument, NULL, OPT_Penalization },
		 { "poly-space", required_argument, NULL, OPT_PolySpace },
		 { "cycle", required_argument, NULL, OPT_MGCycle },
		 { "coarse-size", required_argument, NULL, OPT_CoarseMatrixSize },
		 { "smoothers", required_argument, NULL, OPT_Smoothers },
		 { "cs", required_argument, NULL, OPT_CoarseningStrategy },
		 { "initial-guess", required_argument, NULL, OPT_InitialGuess },
		 { "threads", required_argument, NULL, OPT_Threads },
		 { NULL, 0, NULL, 0 }
	};

	int long_index = 0;
	int option = 0;
	while ((option = getopt_long_only(argc, argv, "d:s:n:b:p:a:l:o:w:g:h", long_opts, &long_index)) != -1)
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
			case OPT_RightHandSide:
				rhsCode = optarg;
				if (   rhsCode.compare("sine")    != 0 
					&& rhsCode.compare("poly")    != 0 
					&& rhsCode.compare("one")     != 0
					&& rhsCode.compare("x")       != 0  
					&& rhsCode.compare("heterog") != 0 
					&& rhsCode.compare("kellogg") != 0)
					argument_error("unknown right-hand side code '" + rhsCode + "'. Check -rhs argument.");
				break;
			case OPT_Heterogeneity: 
				kappa1 = atof(optarg);
				break;
			case OPT_Anisotropy:
				anisotropyRatio = atof(optarg);
				break;
			case OPT_Partition:
				partition = optarg;
				if (partition.compare("halves") != 0 && partition.compare("chiasmus") != 0)
					argument_error("unknown partition '" + partition + "'. Check -partition argument.");
				break;
			case 'n': 
				n = stoul(optarg, nullptr, 0);
				break;
			case OPT_Discretization:
				discretization = optarg;
				if (discretization.compare("dg") != 0 && discretization.compare("hho") != 0)
					argument_error("unknown discretization '" + discretization + "'. Check -discr argument.");
				break;
			case OPT_Mesh:
				meshCode = optarg;
				if (   meshCode.compare("cart") != 0
					&& meshCode.compare("cart-poly") != 0
					&& meshCode.compare("tri") != 0
					&& meshCode.compare("gmsh-tri") != 0
					&& meshCode.compare("quad") != 0
					&& meshCode.compare("quad-poly") != 0)
					argument_error("unknown mesh code '" + meshCode + "'. Check -mesh argument.");
				break;
			case OPT_Stabilization:
				stabilization = optarg;
				if (stabilization.compare("hdg") != 0 && stabilization.compare("hho") != 0)
					argument_error("unknown stabilization code '" + stabilization + "'. Check -stab argument.");
				break;
			case 'b': 
				basisCode = optarg;
				if (basisCode.compare("monomials") != 0 && basisCode.compare("legendre") != 0 && basisCode.compare("nlegendre") != 0 && basisCode.compare("bernstein") != 0 && basisCode.compare("hemker") != 0)
					argument_error("unknown polynomial basis '" + basisCode + "'. Check -b argument.");
				break;
			case 'p': 
				polyDegree = atoi(optarg);
				break;
			case OPT_PolySpace:
			{
				string polySpace = optarg;
				if (polySpace.compare("p") != 0 && polySpace.compare("q") != 0)
					argument_error("unknown polynomial space '" + polySpace + "'. Check -poly-space argument.");
				if (polySpace.compare("q") == 0)
					usePolynomialSpaceQ = true;
				break;
			}
			case OPT_Penalization:
				penalizationCoefficient = atoi(optarg);
				break;
			case OPT_NoStaticCondensation:
				staticCondensation = false;
				break;
			case 'a': 
				a = optarg;
				break;
			case 's': 
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
			case OPT_CoarseMatrixSize:
				matrixMaxSizeForCoarsestLevel = atoi(optarg);
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
			case OPT_CoarseningStrategy:
			{
				string code = optarg;
				if (code.compare("s") != 0 && code.compare("a") != 0)
					argument_error("unknown coarsening strategy code '" + code + "'. Check -cs argument.");
				if (code.compare("a") == 0)
					coarseningStgy = CoarseningStrategy::Agglomeration;
				break;
			}
			case OPT_InitialGuess:
				initialGuessCode = optarg;
				if (initialGuessCode.compare("0") != 0 && initialGuessCode.compare("1") != 0)
					argument_error("unknown initial guess '" + initialGuessCode + "'. Check -initial-guess argument.");
			case OPT_Threads:
				BaseParallelLoop::SetDefaultNThreads(atoi(optarg));
				break;
			case 'o': 
				outputDirectory = optarg;
				break;
			default:
			{
				string character(1, option);
				argument_error("unknown option '" + character + "'.");
				exit(EXIT_FAILURE);
			}
		}
	}

	if (dimension != 1 && rhsCode.compare("heterog") == 0)
		argument_error("-rhs heterog is only supported in 1D.");

	if (rhsCode.compare("kellogg") == 0)
	{
		if (dimension != 2)
			argument_error("-rhs kellogg is only supported in 2D.");
		if (kappa1 != 1)
			cout << "Warning: -heterog argument is ignored due to -rhs kellogg" << endl;
		if (anisotropyRatio != 1)
			cout << "Warning: -haniso argument is ignored due to -rhs kellogg" << endl;
		if (partition.compare("chiasmus") != 0)
			cout << "Warning: -partition argument is ignored due to -rhs kellogg" << endl;

		partition = "chiasmus";
		anisotropyRatio = 1;
		kappa1 = 1;
		kappa2 = 161.4476387975881;
	}

	if (dimension > 1 && discretization.compare("dg") == 0 && polyDegree == 0)
		argument_error("In 2D/3D, DG is not a convergent scheme for p = 0.");

	if (dimension == 1 && discretization.compare("hho") == 0 && polyDegree != 1)
		argument_error("HHO in 1D only exists for p = 1.");

	if (discretization.compare("hho") == 0 && polyDegree == 0)
		argument_error("HHO does not exist with p = 0. Linear approximation at least (p >= 1).");

	if (meshCode.compare("tri") == 0 && dimension != 2)
		argument_error("The triangular mesh in only available in 2D.");

	if (meshCode.compare("quad") == 0 && dimension != 2)
		argument_error("The quadrilateral mesh in only available in 2D.");

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
		else if (a[i] == 'r')
			action |= Action::ComputeL2Error;
		else
		{
			string character(1, a[i]);
			argument_error("unknown action '" + character + "'. Check -a argument.");
		}
	}

	Program* program = nullptr;
	if (dimension == 1)
		program = new ProgramDim<1>();
	else if (dimension == 2)
		program = new ProgramDim<2>();
	else if (dimension == 3)
		program = new ProgramDim<3>();

	program->Start(rhsCode, kappa1, kappa2, anisotropyRatio, partition, 
		n, discretization, meshCode, stabilization, basisCode, polyDegree, usePolynomialSpaceQ, penalizationCoefficient, staticCondensation, action,
		nMultigridLevels, matrixMaxSizeForCoarsestLevel, wLoops, useGalerkinOperator, preSmootherCode, postSmootherCode, nPreSmoothingIterations, nPostSmoothingIterations,
		coarseningStgy, initialGuessCode, outputDirectory, solverCode, solverTolerance);

	delete program;

	cout << "----------------- SUCCESSFUL TERMINATION ----------------" << endl;
    return EXIT_SUCCESS;
}
