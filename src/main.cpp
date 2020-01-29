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
	cout << "-mesh CODE" << endl;
	cout << "      Type of mesh. (Default: cart)" << endl;
	cout << "               cart          - Unit square/cube discretized by an in-house uniform Cartesian mesh (default)" << endl;
	cout << "               tri           - (2D only) Unit square/cube discretized by an in-house uniform trianglular mesh" << endl;
	cout << "               quad          - (2D only) Unit square/cube discretized by an in-house uniform quadrilateral mesh" << endl;
	cout << "               gmsh-cart     - Unit square discretized by a uniform Cartesian mesh built by GMSH" << endl;
	cout << "               gmsh-tri      - (2D only) Unit square discretized by a uniform triangular mesh built by GMSH" << endl;
	cout << "               gmsh-uns-tri  - (2D only) Unit square discretized by an unstructured triangular mesh built by GMSH" << endl;
	cout << "               gmsh          - (2D only) .msh or .geo GMSH file given in the argument -file" << endl;
	cout << endl;
	cout << "-n NUM" << endl;
	cout << "      Number of subdivisions in each cartesian dimension of the unit square/cube (default: 16)." << endl;
	cout << endl;
	cout << "-file PATH" << endl;
	cout << "      GMSH file (.msh or .geo)." << endl;
	cout << endl;
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
	cout << "              pcgmg    - Conjugate Gradient, preconditioned with the custom multigrid for HHO 'mg'" << endl;
	cout << "              agmg     - Yvan Notay's AGMG solver" << endl;
	cout << endl;
	cout << "-initial-guess CODE" << endl;
	cout << "      Initial guess for the iterative solvers." << endl;
	cout << "              0     - zero vector (default)" << endl;
	cout << "              1     - all ones vector" << endl;
	cout << endl;
	cout << "-tol NUM" << endl;
	cout << "      Tolerance of the iterative solver (default: 1e-8)." << endl;
	cout << endl;
	cout << "-max-iter NUM" << endl;
	cout << "      Maximum number of iterations for the iterative solver (default: 200)." << endl;
	cout << endl;
	cout << "----------------------------------------------------------------------" << endl;
	cout << "                                 Multigrid                            " << endl;
	cout << "----------------------------------------------------------------------" << endl;
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
	cout << "      Coarsening strategy of the multigrid." << endl;
	cout << "              s    - standard coarsening (merge colinear faces on the coarse mesh)" << endl;
	cout << "              a    - agglomeration coarsening (keep fine faces on the coarse mesh)" << endl;
	cout << "              r    - fine meshes obtained by structured refinement of the coarse mesh" << endl;
	cout << endl;
	cout << "-prolong NUM" << endl;
	cout << "      How the prolongation operator is built." << endl;
	cout << "              1        - Interpolation from coarse faces to coarse cells" << endl;
	cout << "                         L2-projection on the fine faces" << endl;
	cout << "              2        - Interpolation from coarse faces to coarse cells" << endl;
	cout << "                         On faces present on both fine and coarse meshes, we keep the polynomials identical. Otherwise, L2-projection from the cell polynomials" << endl;
	cout << "              3        - Interpolation from coarse faces to coarse cells" << endl;
	cout << "                         Adjoint of the same interpolation on the fine mesh" << endl;
	cout << endl;
	cout << "-cell-interp NUM" << endl;
	cout << "      In the polongation, degree of the polynomial interpolated on the cells from the faces." << endl;
	cout << "              1        - degree k+1: recover cell unknowns by solving the local problem and apply the local reconstructor (default)" << endl;
	cout << "              2        - degree k  : recover cell unknowns by solving the local problem" << endl;
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
	cout << "              e   - Export system" << endl;
	cout << "              c   - export all Components of the matrix in separate files" << endl;
	cout << "              f   - export Faces for Matlab" << endl;
	cout << "              s   - Solve system" << endl;
	cout << "              m   - export Multigrid matrices" << endl;
	cout << "              v   - export solution Vector (requires 's')" << endl;
	cout << "              r   - compute L2 eRRor against the analytical solution (if known)" << endl;
	cout << "              t   - run unit Tests" << endl;
	cout << endl;
	cout << "-o PATH" << endl;
	cout << "      Output directory to export files (default: ./)." << endl;
	cout << endl;
	cout << "--------------------------------------------------------" << endl;
}

void argument_error(string msg)
{
	cout << Utils::BeginRed << "Argument error: " << msg << Utils::EndColor << endl;
	cout << "------------------------- FAILURE -------------------------" << endl;
	exit(EXIT_FAILURE);
}

int main(int argc, char* argv[])
{
	cout << "-------------------------- START --------------------------" << endl;
	cout << "Option -h for help." << endl;
	cout << "-----------------------------------------------------------" << endl;
	Eigen::initParallel();

	ProgramArguments args;

	enum {
		// Problem
		OPT_RightHandSide = 1000,
		OPT_Heterogeneity,
		OPT_Anisotropy,
		OPT_Partition,
		// Discretization
		OPT_Discretization,
		OPT_Mesh,
		OPT_MeshFilePath,
		OPT_Stabilization,
		OPT_NoStaticCondensation,
		OPT_Penalization,
		OPT_PolySpace,
		// Solver
		OPT_InitialGuess,
		OPT_Tolerance,
		OPT_MaxIterations,
		// Multigrid
		OPT_MGCycle,
		OPT_CellInterpCode,
		OPT_ProlongationCode,
		OPT_CoarseMatrixSize,
		OPT_Smoothers,
		OPT_CoarseningStrategy,
		// Misc
		OPT_Threads
	};

	static struct option long_opts[] = {
		 // Problem
		 { "rhs", required_argument, NULL, OPT_RightHandSide },
		 { "heterog", required_argument, NULL, OPT_Heterogeneity },
		 { "aniso", required_argument, NULL, OPT_Anisotropy },
		 { "partition", required_argument, NULL, OPT_Partition },
		 // Discretization
		 { "discr", required_argument, NULL, OPT_Discretization },
		 { "mesh", required_argument, NULL, OPT_Mesh },
		 { "file", required_argument, NULL, OPT_MeshFilePath },
		 { "stab", required_argument, NULL, OPT_Stabilization },
		 { "no-static-cond", required_argument, NULL, OPT_NoStaticCondensation },
		 { "pen", required_argument, NULL, OPT_Penalization },
		 { "poly-space", required_argument, NULL, OPT_PolySpace },
		 // Solver
		 { "initial-guess", required_argument, NULL, OPT_InitialGuess },
		 { "tol", required_argument, NULL, OPT_Tolerance },
		 { "max-iter", required_argument, NULL, OPT_MaxIterations },
		 // Multigrid
		 { "cycle", required_argument, NULL, OPT_MGCycle },
		 { "cell-interp", required_argument, NULL, OPT_CellInterpCode },
		 { "prolong", required_argument, NULL, OPT_ProlongationCode },
		 { "coarse-size", required_argument, NULL, OPT_CoarseMatrixSize },
		 { "smoothers", required_argument, NULL, OPT_Smoothers },
		 { "cs", required_argument, NULL, OPT_CoarseningStrategy },
		 // Misc
		 { "help", no_argument, NULL, 'h' },
		 { "threads", required_argument, NULL, OPT_Threads },
		 { NULL, 0, NULL, 0 }
	};

	int long_index = 0;
	int option = 0;
	while ((option = getopt_long_only(argc, argv, "d:s:n:b:p:a:l:o:w:g:h", long_opts, &long_index)) != -1)
	{
		switch (option) 
		{
			//-----------------//
			//     Problem     //
			//-----------------//

			case 'd': 
			{
				int dimension = atoi(optarg);
				if (dimension < 1 || dimension > 3)
					argument_error("dimension " + to_string(dimension) + "! Are you kidding?! Stop wasting my time.");
				args.Problem.Dimension = dimension;
				break;
			}
			case OPT_RightHandSide:
			{
				string rhsCode = optarg;
				if (   rhsCode.compare("sine") != 0
					&& rhsCode.compare("poly") != 0
					&& rhsCode.compare("one") != 0
					&& rhsCode.compare("x") != 0
					&& rhsCode.compare("heterog") != 0
					&& rhsCode.compare("kellogg") != 0)
					argument_error("unknown right-hand side code '" + rhsCode + "'. Check -rhs argument.");
				args.Problem.RHSCode = rhsCode;
				break;
			}
			case OPT_Heterogeneity: 
				args.Problem.Kappa1 = atof(optarg);
				break;
			case OPT_Anisotropy:
				args.Problem.AnisotropyRatio = atof(optarg);
				break;
			case OPT_Partition:
			{
				string partition = optarg;
				if (partition.compare("halves") != 0 && partition.compare("chiasmus") != 0)
					argument_error("unknown partition '" + partition + "'. Check -partition argument.");
				args.Problem.Partition = partition;
				break;
			}

			//--------------------//
			//   Discretization   //
			//--------------------//

			case 'n': 
				args.Discretization.N = stoul(optarg, nullptr, 0);
				break;
			case OPT_Discretization:
			{
				string discretization = optarg;
				if (discretization.compare("dg") != 0 && discretization.compare("hho") != 0)
					argument_error("unknown discretization '" + discretization + "'. Check -discr argument.");
				args.Discretization.Method = discretization;
				break;
			}
			case OPT_Mesh:
			{
				string meshCode = optarg;
				if (   meshCode.compare("cart") != 0
					&& meshCode.compare("cart-poly") != 0
					&& meshCode.compare("tri") != 0
					&& meshCode.compare("gmsh-cart") != 0
					&& meshCode.compare("gmsh-tri") != 0
					&& meshCode.compare("gmsh-quad") != 0
					&& meshCode.compare("gmsh-uns-tri") != 0
					&& meshCode.compare("gmsh-tetra") != 0
					&& meshCode.compare("gmsh") != 0
					&& meshCode.compare("quad") != 0
					&& meshCode.compare("quad-poly") != 0)
					argument_error("unknown mesh code '" + meshCode + "'. Check -mesh argument.");
				args.Discretization.MeshCode = meshCode;
				break;
			}
			case OPT_MeshFilePath:
				args.Discretization.MeshFilePath = optarg;
				break;
			case OPT_Stabilization:
			{
				string stabilization = optarg;
				if (stabilization.compare("hdg") != 0 && stabilization.compare("hho") != 0)
					argument_error("unknown stabilization code '" + stabilization + "'. Check -stab argument.");
				args.Discretization.Stabilization = stabilization;
				break;
			}
			case 'b': 
			{
				string basisCode = optarg;
				if (basisCode.compare("monomials") != 0 && basisCode.compare("legendre") != 0 && basisCode.compare("nlegendre") != 0 && basisCode.compare("bernstein") != 0 && basisCode.compare("hemker") != 0)
					argument_error("unknown polynomial basis '" + basisCode + "'. Check -b argument.");
				args.Discretization.BasisCode = basisCode;
				break;
			}
			case 'p': 
				args.Discretization.PolyDegree = atoi(optarg);
				break;
			case OPT_PolySpace:
			{
				string polySpace = optarg;
				if (polySpace.compare("p") != 0 && polySpace.compare("q") != 0)
					argument_error("unknown polynomial space '" + polySpace + "'. Check -poly-space argument.");
				if (polySpace.compare("q") == 0)
					args.Discretization.UsePolynomialSpaceQ = true;
				break;
			}
			case OPT_Penalization:
				args.Discretization.PenalizationCoefficient = atoi(optarg);
				break;
			case OPT_NoStaticCondensation:
				args.Discretization.StaticCondensation = false;
				break;

			//------------//
			//   Solver   //
			//------------//

			case 's': 
				args.Solver.SolverCode = optarg;
				break;
			case OPT_InitialGuess:
			{
				string initialGuessCode = optarg;
				if (initialGuessCode.compare("0") != 0 && initialGuessCode.compare("1") != 0)
					argument_error("unknown initial guess '" + initialGuessCode + "'. Check -initial-guess argument.");
				args.Solver.InitialGuessCode = initialGuessCode;
				break;
			}
			case OPT_Tolerance:
				args.Solver.Tolerance = atof(optarg);
				break;
			case OPT_MaxIterations:
				args.Solver.MaxIterations = atoi(optarg);
				if (args.Solver.MaxIterations < 1)
					argument_error("-max-iter argument must be > 0.");
				break;

			//---------------//
			//   Multigrid   //
			//---------------//

			case OPT_MGCycle:
			{
				string s(optarg);
				regex pattern("^([vwVW]),([[:digit:]]+),([[:digit:]]+)$");
				smatch matches;

				if (std::regex_search(s, matches, pattern))
				{
					char cycleLetter = matches.str(1)[0];
					if (cycleLetter == 'v' || cycleLetter == 'V')
						args.Solver.MG.WLoops = 1;
					else if (cycleLetter == 'w' || cycleLetter == 'W')
						args.Solver.MG.WLoops = 2;
					args.Solver.MG.PreSmoothingIterations = stoi(matches.str(2));
					args.Solver.MG.PostSmoothingIterations = stoi(matches.str(3));
				}
				else
					argument_error("syntax error in the multigrid cycle. Check -cycle argument.");
				break;
			}
			case OPT_CellInterpCode:
			{
				int cellInterp = atoi(optarg);
				if (cellInterp != 1 && cellInterp != 2)
					argument_error("check -cell-interp argument. Expecting 1 or 2.");
				args.Solver.MG.CellInterpolationCode = cellInterp;
				break;
			}
			case OPT_ProlongationCode:
			{
				int prolongationCode = atoi(optarg);
				if (prolongationCode != 1 && prolongationCode != 2 && prolongationCode != 3)
					argument_error("check -prolong argument. Expecting 1, 2 or 3.");
				args.Solver.MG.ProlongationCode = prolongationCode;
				break;
			}
			case 'l': 
				args.Solver.MG.Levels = atoi(optarg);
				break;
			case OPT_CoarseMatrixSize:
				args.Solver.MG.MatrixMaxSizeForCoarsestLevel = atoi(optarg);
				break;
			case 'w': 
				args.Solver.MG.WLoops = atoi(optarg);
				if (args.Solver.MG.WLoops < 1)
					argument_error("the number of loops in the W-cycle must be >= 1. Check -w argument.");
				break;
			case 'g': 
				args.Solver.MG.UseGalerkinOperator = atoi(optarg);
				break;
			case OPT_Smoothers:
			{
				string s(optarg);
				auto pos = s.find(",");
				args.Solver.MG.PreSmootherCode = s.substr(0, s.find(","));
				args.Solver.MG.PostSmootherCode = s.substr(pos + 1);
				break;
			}
			case OPT_CoarseningStrategy:
			{
				string coarseningStgyCode = optarg;
				if (coarseningStgyCode.compare("s") != 0 && coarseningStgyCode.compare("a") != 0 && coarseningStgyCode.compare("r") != 0)
					argument_error("unknown coarsening strategy code '" + coarseningStgyCode + "'. Check -cs argument.");
				args.Solver.MG.CoarseningStgyCode = coarseningStgyCode;
				break;
			}

			//----------------//
			//      Misc      //
			//----------------//

			case 'h':
				print_usage();
				exit(EXIT_SUCCESS);
				break;
			case 'a':
				args.ActionCodes = optarg;
				break;
			case OPT_Threads:
				BaseParallelLoop::SetDefaultNThreads(atoi(optarg));
				break;
			case 'o': 
				args.OutputDirectory = optarg;
				break;
			default:
			{
				string character(1, option);
				argument_error("unknown option '" + character + "'.");
			}
		}
	}

	if (args.Problem.Dimension != 1 && args.Problem.RHSCode.compare("heterog") == 0)
		argument_error("-rhs heterog is only supported in 1D.");

	if (args.Problem.RHSCode.compare("kellogg") == 0)
	{
		if (args.Problem.Dimension != 2)
			argument_error("-rhs kellogg is only supported in 2D.");
		if (args.Problem.Kappa1 != 1)
			cout << "Warning: -heterog argument is ignored due to -rhs kellogg" << endl;
		if (args.Problem.AnisotropyRatio != 1)
			cout << "Warning: -haniso argument is ignored due to -rhs kellogg" << endl;
		if (args.Problem.Partition.compare("chiasmus") != 0)
			cout << "Warning: -partition argument is ignored due to -rhs kellogg" << endl;

		args.Problem.Partition = "chiasmus";
		args.Problem.AnisotropyRatio = 1;
		args.Problem.Kappa1 = 1;
		args.Problem.Kappa2 = 161.4476387975881;
	}

	if (args.Problem.Dimension > 1 && args.Discretization.Method.compare("dg") == 0 && args.Discretization.PolyDegree == 0)
		argument_error("In 2D/3D, DG is not a convergent scheme for p = 0.");

	if (args.Problem.Dimension == 1 && args.Discretization.Method.compare("hho") == 0 && args.Discretization.PolyDegree != 1)
		argument_error("HHO in 1D only exists for p = 1.");

	if (args.Discretization.Method.compare("hho") == 0 && args.Discretization.PolyDegree == 0)
		argument_error("HHO does not exist with p = 0. Linear approximation at least (p >= 1).");

#ifndef GMSH_ENABLED
	if (args.Discretization.MeshCode.find("gmsh") != string::npos)
		argument_error("GMSH is disabled. Recompile with the cmake option -DENABLE_GMSH=ON to use GMSH meshes, or choose another argument for -mesh.");
#endif // GMSH_ENABLED

	if ((args.Discretization.MeshCode.compare("tri") == 0 || args.Discretization.MeshCode.compare("gmsh-tri") == 0 || args.Discretization.MeshCode.compare("gmsh-uns-tri") == 0) && args.Problem.Dimension != 2)
		argument_error("The triangular mesh in only available in 2D.");

	if ((args.Discretization.MeshCode.compare("quad") == 0 || args.Discretization.MeshCode.compare("gmsh-quad") == 0) && args.Problem.Dimension != 2)
		argument_error("The quadrilateral mesh in only available in 2D.");

	if (args.Discretization.MeshCode.compare("gmsh-tetra") == 0 && args.Problem.Dimension != 3)
		argument_error("The tetrahedral mesh in only available in 3D.");

	if (args.Discretization.MeshCode.compare("gmsh") == 0 && args.Discretization.MeshFilePath.compare("") == 0)
		argument_error("The GMSH file path is missing. Add the argument -file.");

#ifndef AGMG_ENABLED
	if (args.Solver.SolverCode.compare("agmg") == 0)
		argument_error("AGMG is disabled. Recompile with the cmake option -DENABLE_AGMG=ON, or choose another solver.");
#endif // AGMG_ENABLED

	if (args.Solver.SolverCode.compare("mg") == 0 && args.Discretization.Method.compare("dg") == 0)
		argument_error("Multigrid only applicable on HHO discretization.");

	if (args.Solver.SolverCode.compare("mg") == 0 && args.Discretization.Method.compare("hho") == 0 && !args.Discretization.StaticCondensation)
		argument_error("Multigrid only applicable if the static condensation is enabled.");


	args.Solver.MG.CoarseningStgy = CoarseningStrategy::Standard;
	if (args.Solver.MG.CoarseningStgyCode.compare("default") == 0)
	{
		if (args.Discretization.MeshCode.find("gmsh") == 0)
			args.Solver.MG.CoarseningStgy = CoarseningStrategy::StructuredRefinement;
		else
			args.Solver.MG.CoarseningStgy = CoarseningStrategy::Standard;
	}
	else if (args.Solver.MG.CoarseningStgyCode.compare("s") == 0)
		args.Solver.MG.CoarseningStgy = CoarseningStrategy::Standard;
	else if (args.Solver.MG.CoarseningStgyCode.compare("a") == 0)
		args.Solver.MG.CoarseningStgy = CoarseningStrategy::Agglomeration;
	else if (args.Solver.MG.CoarseningStgyCode.compare("r") == 0)
		args.Solver.MG.CoarseningStgy = CoarseningStrategy::StructuredRefinement;


	args.Actions = Action::LogAssembly;
	for (size_t i = 0; i < args.ActionCodes.length(); i++)
	{
		if (args.ActionCodes[i] == 'e')
			args.Actions |= Action::ExtractSystem;
		else if (args.ActionCodes[i] == 'c')
			args.Actions |= Action::ExtractComponentMatrices;
		else if (args.ActionCodes[i] == 'f')
			args.Actions |= Action::ExportFaces;
		else if (args.ActionCodes[i] == 's')
			args.Actions |= Action::SolveSystem;
		else if (args.ActionCodes[i] == 'm')
			args.Actions |= Action::ExportMultigridMatrices;
		else if (args.ActionCodes[i] == 'v')
			args.Actions |= Action::ExtractSolution;
		else if (args.ActionCodes[i] == 'r')
			args.Actions |= Action::ComputeL2Error;
		else if (args.ActionCodes[i] == 't')
			args.Actions |= Action::UnitTests;
		else
		{
			string character(1, args.ActionCodes[i]);
			argument_error("unknown action '" + character + "'. Check -a argument.");
		}
	}

	Program* program = nullptr;
	if (args.Problem.Dimension == 1)
		program = new ProgramDim<1>();
	else if (args.Problem.Dimension == 2)
		program = new ProgramDim<2>();
	else if (args.Problem.Dimension == 3)
		program = new ProgramDim<3>();

	program->Start(args);

	delete program;

	cout << "----------------- SUCCESSFUL TERMINATION ----------------" << endl;
    return EXIT_SUCCESS;
}
