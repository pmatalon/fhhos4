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
	cout << "-geo CODE" << endl;
	cout << "      Geometry: in-house test cases, or imported from GMSH." << endl;
	cout << "            segment           - Unit segment (1D)" << endl;
	cout << "            square            - Unit square" << endl;
	cout << "            square4quadrants  - Unit square divided into 4 quadrants" << endl;
	cout << "            cube              - Unit cube" << endl;
	cout << "            <file>            - GMSH .geo or .msh file" << endl;
	cout << "                                Use relative or absolute path, or simply the file name if the file is stored in the folder data/mesh/" << endl;
	cout << endl;
	cout << "-tc CODE" << endl;
	cout << "      Test case code. The test case defines the source function of the problem, the diffusion field, as well as available predefined boundary conditions." << endl;
	cout << "      In some, it also determines the analytical solution in the homogeneous isotropic case so that the L2 error can be computed." << endl;
	cout << "      The following test cases are predefined and available for the simple geometries listed in the -geo argument." << endl;
	cout << "               sine    - the source function and the analytical solution are a sine functions" << endl;
	cout << "               poly    - the source function is constant, the analytical solution is a polynomial of total degree 2*d" << endl;
	cout << "               zero    - the source function and the analytical solution are 0" << endl;
	cout << "               one     - the source function is 0, the analytical solution is 1" << endl;
	cout << "               x       - the source function is 0, the analytical solution is x" << endl;
	cout << "               heterog - (1D only) heterogeneous diffusion-specific analytical solution" << endl;
	cout << "               kellogg - (2D only) heterogeneous diffusion-specific analytical solution (known benchmark)" << endl;
	cout << "      When the geometry given is a GMSH file, a test case with same code as the file name must be defined." << endl;
	cout << endl;
	cout << "-bc CODE" << endl;
	cout << "      Boundary conditions, according to what the selected test case allows." << endl;
	cout << "            d          - Dirichlet (default)" << endl;
	cout << "            m          - Mixed Neumann-Dirichlet" << endl;
	cout << "            <other>    - Test case-specific boundary conditions" << endl;
	cout << endl;
	cout << "-heterog NUM" << endl;
	cout << "      Heterogeneity ratio. Has an effect only if multiple physical parts are defined in the geometry." << endl;
	cout << "      Ex:       1             - homogeneous diffusion (default)" << endl;
	cout << "                0 < NUM < 1   - to be used with -tc heterog" << endl;
	cout << "                1e4           - allows to set the order of magnitude of the ratio" << endl;
	cout << endl;
	cout << "-aniso NUM" << endl;
	cout << "      Anisotropy ratio (valid across the whole domain)." << endl;
	cout << "      Example: 1e4" << endl;
	cout << endl;
	cout << "-aniso-angle NUM" << endl;
	cout << "      Anisotropy angle expressed in degrees (valid across the whole domain)." << endl;
	cout << "      Example: 45" << endl;
	cout << endl;
	cout << "----------------------------------------------------------------------" << endl;
	cout << "                                  Mesh                                " << endl;
	cout << "----------------------------------------------------------------------" << endl;
	cout << endl;
	cout << "-mesh CODE" << endl;
	cout << "      Type of mesh (default: unstructured simplicial)." << endl;
	cout << "            cart      - Cartesian" << endl;
	cout << "            stri      - Structured triangular" << endl;
	cout << "            tri       - Unstructured triangular" << endl;
	cout << "            quad      - Unstructured quadrilateral (see also argument -stretch)" << endl;
	cout << "            stetra    - Structured tetrahedral (embedded in a Cartesian mesh)" << endl;
	cout << "            tetra     - Unstructured tetrahedral" << endl;
	cout << endl;
	cout << "-mesher CODE" << endl;
	cout << "      Mesher used: in-house or imported file from GMSH." << endl;
	cout << "            inhouse   - in-house (default for geometries 'segment', 'square', 'cube')" << endl;
	cout << "            gmsh      - GMSH (default for imported files)" << endl;
	cout << endl;
	cout << "-n NUM" << endl;
	cout << "      Number of subdivisions in each cartesian dimension of the unit square/cube (default: 16)." << endl;
	cout << "      If used in conjonction with -ny or -nz, then the value is used for -nx." << endl;
	cout << endl;
	cout << "-nx NUM" << endl;
	cout << "      Number of subdivisions in the x dimension of the unit square/cube (default: 16)." << endl;
	cout << endl;
	cout << "-ny NUM" << endl;
	cout << "      Number of subdivisions in the y dimension of the unit square/cube (default: same as -n or -nx). Used only for the meshes 'cart', 'tri', 'quad'." << endl;
	cout << endl;
	cout << "-nz NUM" << endl;
	cout << "      Number of subdivisions in the z dimension of the unit square/cube (default: same as -n or -nx). Used only for the meshes 'cart', 'tri', 'quad'." << endl;
	cout << endl;
	cout << "-stretch NUM" << endl;
	cout << "      0 <= NUM < 1   stretching factor of the quadrilateral element used in the 'quad' mesh of the 'inhouse' mesher." << endl;
	cout << "      0 yields square elements. Default is 0.5: the top right corner is moved to the right 0.5 times the element's width." << endl;
	cout << endl;
	cout << "----------------------------------------------------------------------" << endl;
	cout << "                             Discretization                           " << endl;
	cout << "----------------------------------------------------------------------" << endl;
	cout << endl;
	cout << "-discr CODE" << endl;
	cout << "      Discretization method (default: hho)." << endl;
	cout << "               dg     - Discontinuous Galerkin (Symmetric Interior Penalty)" << endl;
	cout << "               hho    - Hybrid High-Order" << endl;
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
	cout << "      Linear solver." << endl;
	cout << "              default   - 'lu' for small problems, 'eigencg' for bigger problems, 'mg' for HHO" << endl;
	cout << "              lu        - LU factorization (Eigen library)" << endl;
	cout << "              cg        - Conjugate Gradient, no preconditioner" << endl;
	cout << "              eigencg   - Conjugate Gradient (Eigen library) with diagonal preconditioner" << endl;
	cout << "              fcg       - Flexible Conjugate Gradient (truncation-restart: FCG(1)), no preconditioner" << endl;
	cout << "              [b]j      - [Block] Jacobi" << endl;
	cout << "              [r][b]gs  - [Reverse order] [Block] Gauss-Seidel" << endl;
	cout << "              [r][b]sor - [Reverse order] [Block] SOR (identical to gs, both can use argument -relax)" << endl;
	cout << "              mg        - Custom multigrid for HHO" << endl;
	cout << "              cgmg      - Conjugate Gradient, preconditioned with the custom multigrid for HHO 'mg'" << endl;
	cout << "              fcgmg     - Flexible Conjugate Gradient FCG(1), preconditioned with the custom multigrid for HHO 'mg' (meant to be used with K-cycle)" << endl;
	cout << "              agmg      - Yvan Notay's AGMG solver" << endl;
	cout << "      For the block solvers, the block size is set to the number of DOFs per cell (DG) or face (HHO)." << endl;
	cout << "      Jacobi, Gauss-Seidel and SOR can use argument -relax to change the relaxation parameter." << endl;
	cout << endl;
	cout << "-initial-guess CODE" << endl;
	cout << "      Initial guess for the iterative solvers." << endl;
	cout << "              0        - zero vector (default)" << endl;
	cout << "              1        - all ones vector" << endl;
	cout << "              rand     - random" << endl;
	cout << "              smooth   - initial guess which generates a smooth error w.r.t. the multigrid prolongation" << endl;
	cout << endl;
	cout << "-tol NUM" << endl;
	cout << "      Tolerance of the iterative solver (default: 1e-8)." << endl;
	cout << endl;
	cout << "-max-iter NUM" << endl;
	cout << "      Maximum number of iterations for the iterative solver (default: 200)." << endl;
	cout << endl;
	cout << "-relax NUM" << endl;
	cout << "      Relaxation parameter in (0,2) used for the Jacobi, Gauss-Seidel, SOR and their derived versions (default: 1)." << endl;
	cout << endl;
	cout << "----------------------------------------------------------------------" << endl;
	cout << "                                 Multigrid                            " << endl;
	cout << "----------------------------------------------------------------------" << endl;
	cout << endl;
	cout << "-cycle (V|W|K),NUM,NUM[,+-*NUM]" << endl;
	cout << "      Multigrid cycle." << endl;
	cout << "      Examples: \"V,1,0\", \"W,1,1\", \"K,1,1\"." << endl;
	cout << "      The third number defines the value of which the number of smoothing iterations is increased or decreased at each coarse level," << endl;
	cout << "      according to the arithmetic operation ('+', '-' or '*') that preceeds it." << endl;
	cout << "      Example: \"V,1,1,+1\" performs (1,1) smoothing steps at the finer level, then (2,2) on the next one, then (3,3), etc." << endl;
	cout << "      The K-cycle is usually meant to be used as a preconditioner for a Flexible Conjugate Gradient." << endl;
	cout << "      Example: \"-s fgcmg -cycle K,1,1\"." << endl;
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
	cout << "              0    - discretized operator (default)" << endl;
	cout << "              1    - Galerkin operator" << endl;
	cout << endl;
	cout << "-smoothers CODE,CODE" << endl;
	cout << "      Pre-smoother,post-smoother: \"gs,rgs\" for example." << endl;
	cout << "              j    - Jacobi" << endl;
	cout << "              gs   - Gauss-Seidel" << endl;
	cout << "              rgs  - Reverse Gauss-Seidel" << endl;
	cout << "              bj   - Block Jacobi: the block size is set to the number of DOFs per face" << endl;
	cout << "              bj23 - Block Jacobi with a damping factor of 2/3" << endl;
	cout << "              bgs  - Block Gauss-Seidel: the block size is set to the number of DOFs per face" << endl;
	cout << "              rbgs - Reverse Block Gauss-Seidel" << endl;
	cout << endl;
	cout << "-cs CODE" << endl;
	cout << "      Coarsening strategy of the multigrid." << endl;
	cout << "              s   - Standard coarsening (merge colinear faces on the coarse mesh)" << endl;
	cout << "              a   - Agglomeration coarsening (keep fine faces on the coarse mesh)" << endl;
	cout << "              l   - (Experimental) agglomeration coarsening by most collinear/coplanar faces (non-nested!)" << endl;
	cout << "              c   - (Experimental) agglomeration coarsening by closest center (non-nested!)" << endl;
	cout << "              g   - (Experimental) agglomeration coarsening by closest face (non-nested!)" << endl;
	cout << "              i   - (Experimental) agglomeration coarsening by largest interface (non-nested)" << endl;
	cout << "              p   - (Experimental) agglomeration coarsening by seed points (non-nested!)" << endl;
	cout << "              n   - (Experimental) agglomeration coarsening by face neighbours (non-nested!)" << endl;
	cout << "              v   - (Experimental) agglomeration coarsening by vertex neighbours (non-nested!)" << endl;
	cout << "              m   - Independant remeshing by GMSH with double the mesh size (non-nested!)." << endl;
	cout << "              f   - Face coarsening: the faces are coarsened and all kept on the coarse skeleton. Requires -g 1." << endl;
	cout << "              r   - Fine meshes obtained by structured refinement of the coarse mesh using GMSH's splitting method" << endl;
	cout << "              b   - Fine meshes obtained by structured refinement of the coarse mesh using the Bey method" << endl;
	cout << endl;
	cout << "-coarse-n NUM" << endl;
	cout << "      If a refinement strategy is used, sets the mesh size of the starting coarse mesh." << endl;
	cout << "      The chosen value will set the variable N defined at the beginning of the GMSH .geo file." << endl;
	cout << endl;
	cout << "-prolong NUM" << endl;
	cout << "      How the prolongation operator is built." << endl;
	cout << "              " << (unsigned)Prolongation::CellInterp_Trace << "  - ";
	cout <<                    "Step 1: Interpolation from coarse faces to coarse cells (refer to -cell-interp argument)" << endl;
	cout << "                   Step 2: Trace on the fine faces" << endl;
	cout << "              " << (unsigned)Prolongation::CellInterp_InjectAndTrace << "  - ";
	cout <<                    "Step 1: Interpolation from coarse faces to coarse cells (refer to -cell-interp argument)" << endl;
	cout << "                   Step 2: On faces present on both fine and coarse meshes, we keep the polynomials identical." << endl;
	cout << "                           On faces interior to coarse elements, trace of the cell polynomials." << endl;
	cout << "              " << (unsigned)Prolongation::CellInterp_Inject_Adjoint << "  - ";
	cout <<                    "Step 1: Interpolation from coarse faces to coarse cells (refer to -cell-interp argument)" << endl;
	cout << "                   Step 2: Canonical injection from coarse to fine cells" << endl;
	cout << "                   Step 3: Adjoint of the cell interpolation on the fine mesh" << endl;
	cout << "              " << (unsigned)Prolongation::Wildey << "  - ";
	cout <<                    "Algorithm from Wildey et al.: the coarse level is built by static condensation of the fine faces interior to coarse elements." << endl;
	cout << "                   The prolongation solves those condensed unknowns." << endl;
	cout << "                   To reproduce Wildey et al.'s algorithm, this option should be used with '-g 1 -cycle V,1,1,*2'." << endl;
	cout << "              " << (unsigned)Prolongation::FaceInject << "  - ";
	cout <<                    "Canonical injection from coarse faces to fine faces (implemented to be used with option -cs f)." << endl;
	cout << "              " << (unsigned)Prolongation::CellInterp_Inject_Trace << "  - ";
	cout <<                    "Same as 1, but another implementation:" << endl;
	cout << "                   Step 1: Interpolation from coarse faces to coarse cells (refer to -cell-interp argument)" << endl;
	cout << "                   Step 2: Canonical injection from coarse to fine cells" << endl;
	cout << "                   Step 3: Trace on the fine faces" << endl;
	cout << "              " << (unsigned)Prolongation::CellInterp_L2proj_Trace << "  - ";
	cout <<                    "Non-nested variant of " << (unsigned)Prolongation::CellInterp_Trace << " and " << (unsigned)Prolongation::CellInterp_Inject_Trace << ":" << endl;
	cout << "                   Step 1: Interpolation from coarse faces to coarse cells (refer to -cell-interp argument)" << endl;
	cout << "                   Step 2: L2-projection onto the fine cells" << endl;
	cout << "                   Step 3: Trace on the fine faces" << endl;
	cout << "              " << (unsigned)Prolongation::CellInterp_ApproxL2proj_Trace << "  - ";
	cout <<                    "Variant of " << (unsigned)Prolongation::CellInterp_L2proj_Trace << " where the L2-projection is not computed exactly but has the same approximation properties." << endl;
	cout << endl;
	cout << "-cell-interp NUM" << endl;
	cout << "      In the polongation, degree of the polynomial interpolated on the cells from the faces." << endl;
	cout << "              1   - degree k+1: recover cell unknowns by solving the local problem and apply the local reconstructor (default)" << endl;
	cout << "              2   - degree k  : recover cell unknowns by solving the local problem" << endl;
	cout << endl;
	cout << "-weight CODE" << endl;
	cout << "      In the polongation, weighting factor of the projection to the fine faces." << endl;
	cout << "              k   - proportional to the diffusion coefficient (default)" << endl;
	cout << "              a   - simple average (1/2)" << endl;
	cout << endl;
	cout << "----------------------------------------------------------------------" << endl;
	cout << "                             Miscellaneous                            " << endl;
	cout << "----------------------------------------------------------------------" << endl;
	cout << endl;
	cout << "-h, -help" << endl;
	cout << "      Prints usage." << endl;
	cout << endl;
	cout << "-threads NUM" << endl;
	cout << "      Max number of threads used for parallelism (default: 0)." << endl;
	cout << "              0     - automatic: usually number of cores x 2 (because of hyper-threading)" << endl;
	cout << "              1     - sequential execution" << endl;
	cout << "              other - requested number of threads" << endl;
	cout << endl;
	cout << "-export CODES" << endl;
	cout << "      Data export. The CODES must be comma-separated." << endl;
	cout << "              lsys    - linear system" << endl;
	cout << "              amat    - matrices decomposed by assembly terms in separate files (consistency, stabilization, mass, etc.)" << endl;
	cout << "              mesh    - mesh to be used in Matlab" << endl;
	cout << "              solvect - solution vector(s)" << endl;
	cout << "              solgmsh - solution files (.pos and .msh) to be used in GMSH for visualization" << endl;
	cout << "              mg      - Multigrid components (intergrid operator matrices, coarse meshes, etc.)" << endl;
	cout << endl;
	cout << "-o PATH" << endl;
	cout << "      Output directory to export files." << endl;
	cout << endl;
	cout << "-not-solve" << endl;
	cout << "      Do not solve the linear system." << endl;
	cout << endl;
	cout << "-no-cache" << endl;
	cout << "      Do not use the cached meshes." << endl;
	cout << endl;
	cout << "-gmsh-log" << endl;
	cout << "      Enable GMSH to log in the console." << endl;
	cout << endl;
	cout << "-ut" << endl;
	cout << "      Run unit tests." << endl;
	cout << endl;
	cout << "----------------------------------------------------------------------" << endl;
	cout << "                      Examples and typical use cases                  " << endl;
	cout << "----------------------------------------------------------------------" << endl;
	cout << endl;
	cout << "First example setting the geometry, the type of mesh and its granularity, and the degree of approximation." << endl;
	cout << "Homogeneous diffusion problem in the unit square, using a Cartesian mesh 16x16 and providing a quadratic approximation:" << endl;
	cout << "              -geo square -mesh cart -n 16 -p 2" << endl;
	cout << endl;
	cout << "Heterogeneous problem in the unit square. The domain partition describing the heterogeneity pattern is a chiasmus (4 quadrants)." << endl;
	cout << "The heterogeneity ratio between the two subdomains is set to 1e4." << endl;
	cout << "              -geo square4quadrants -heterog 1e4" << endl;
	cout << "Other use case: the heterogeneity ratio, as well as as non-homogeneous Dirichlet conditions, are set such that it corresponds to a Kellogg problem." << endl;
	cout << "              -geo square4quadrants -tc kellogg" << endl;
	cout << endl;
	cout << "Import a GMSH file describing the geometry (.geo) or the mesh (.msh)." << endl;
	cout << "You can use relative or absolute path, or simply the file name if the file is stored in /data/mesh/." << endl;
	cout << "              -geo file_in_data_mesh_folder.geo" << endl;
	cout << "              -geo /path/to/my/file.geo" << endl;
	cout << endl;
	cout << "Configure the parameter k the HHO discretization." << endl;
	cout << "The argument -p allows to set the desired polynomial degree of the approximation. In HHO, it corresponds to the degree after higher-order reconstruction, i.e. p=k+1." << endl;
	cout << "So setting k=2 is done by" << endl;
	cout << "              -p 3" << endl;
	cout << endl;
	cout << "Use DG-SIPG discretization instead of HHO:" << endl;
	cout << "              -discr dg" << endl;
	cout << endl;
	cout << "Choose the linear solver, for instance the Multigrid developed for HHO:" << endl;
	cout << "              -s mg" << endl;
	cout << endl;
	cout << "Set the Multigrid cycle to V(1,3), the pre- and post-smoothers to the Block-Jacobi iteration using 2/3 as relaxation parameter:" << endl;
	cout << "              -cycle V,1,3 -smoothers bj23,bj23" << endl;
	cout << endl;
	cout << "For unstructured meshes, use a coarsening strategy which also coarsens the faces (-cs n). This will build a hierarchy of non-nested meshes." << endl;
	cout << "To be compliant, you must set the prolongation operator to its non-nested version (which uses the L2-projection from the coarse cells onto the fine ones)." << endl;
	cout << "              -cs n -prolong 7" << endl;
	cout << endl;
	cout << "Use the multigrid as a preconditioner for a Flexible Conjugate Gradient, along with the K-cycle:" << endl;
	cout << "              -s fcgmg -cycle K,1,1" << endl;
	cout << endl;
	cout << "----------------------------------------------------------------------" << endl;
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
	Eigen::initParallel();
	CGALWrapper::Configure();

	bool defaultCycle = true;

	ProgramArguments args;
	args.OutputDirectory = FileSystem::RootPath() + "/out";
	FileSystem::CreateDirectoryIfNotExist(args.OutputDirectory);

	enum {
		// Problem
		OPT_Geometry = 1000,
		OPT_TestCase,
		OPT_RightHandSide, // deprecated
		OPT_BoundaryConditions,
		OPT_HeterogeneityRatio,
		OPT_AnisotropyRatio,
		OPT_AnisotropyAngle,
		// Mesh
		OPT_Mesh,
		OPT_Mesher,
		OPT_Nx,
		OPT_Ny,
		OPT_Nz,
		OPT_Stretch,
		// Discretization
		OPT_Discretization,
		OPT_Stabilization,
		OPT_NoStaticCondensation,
		OPT_Penalization,
		OPT_PolySpace,
		// Solver
		OPT_InitialGuess,
		OPT_Tolerance,
		OPT_MaxIterations,
		OPT_Relaxation,
		// Multigrid
		OPT_MGCycle,
		OPT_CellInterpCode,
		OPT_Weight,
		OPT_ProlongationCode,
		OPT_CoarseMatrixSize,
		OPT_Smoothers,
		OPT_CoarseningStrategy,
		OPT_CoarseN,
		// Misc
		OPT_Threads,
		OPT_Export,
		OPT_DoNotSolve,
		OPT_NoCache,
		OPT_UnitTests,
		OPT_GMSHLog
	};

	static struct option long_opts[] = {
		 // Problem
		 { "geo", required_argument, NULL, OPT_Geometry },
		 { "tc", required_argument, NULL, OPT_TestCase },
		 { "rhs", required_argument, NULL, OPT_RightHandSide }, // deprecated
		 { "bc", required_argument, NULL, OPT_BoundaryConditions },
		 { "heterog", required_argument, NULL, OPT_HeterogeneityRatio },
		 { "aniso", required_argument, NULL, OPT_AnisotropyRatio },
		 { "aniso-angle", required_argument, NULL, OPT_AnisotropyAngle },
		 // Mesh
		 { "mesh", required_argument, NULL, OPT_Mesh },
		 { "mesher", required_argument, NULL, OPT_Mesher },
		 { "nx", required_argument, NULL, OPT_Nx },
		 { "ny", required_argument, NULL, OPT_Ny },
		 { "stretch", required_argument, NULL, OPT_Stretch },
		 // Discretization
		 { "discr", required_argument, NULL, OPT_Discretization },
		 { "stab", required_argument, NULL, OPT_Stabilization },
		 { "no-static-cond", required_argument, NULL, OPT_NoStaticCondensation },
		 { "pen", required_argument, NULL, OPT_Penalization },
		 { "poly-space", required_argument, NULL, OPT_PolySpace },
		 // Solver
		 { "initial-guess", required_argument, NULL, OPT_InitialGuess },
		 { "tol", required_argument, NULL, OPT_Tolerance },
		 { "max-iter", required_argument, NULL, OPT_MaxIterations },
		 { "relax", required_argument, NULL, OPT_Relaxation },
		 // Multigrid
		 { "cycle", required_argument, NULL, OPT_MGCycle },
		 { "cell-interp", required_argument, NULL, OPT_CellInterpCode },
		 { "weight", required_argument, NULL, OPT_Weight },
		 { "prolong", required_argument, NULL, OPT_ProlongationCode },
		 { "coarse-size", required_argument, NULL, OPT_CoarseMatrixSize },
		 { "smoothers", required_argument, NULL, OPT_Smoothers },
		 { "cs", required_argument, NULL, OPT_CoarseningStrategy },
		 { "coarse-n", required_argument, NULL, OPT_CoarseN },
		 // Misc
		 { "help", no_argument, NULL, 'h' },
		 { "threads", required_argument, NULL, OPT_Threads },
		 { "export", required_argument, NULL, OPT_Export },
		 { "not-solve", no_argument, NULL, OPT_DoNotSolve },
		 { "no-cache", no_argument, NULL, OPT_NoCache },
		 { "ut", no_argument, NULL, OPT_UnitTests },
		 { "gmsh-log", no_argument, NULL, OPT_GMSHLog },
		 { NULL, 0, NULL, 0 }
	};

	//------------------------------------------//
	//            Get user arguments            //
	//------------------------------------------//

	int long_index = 0;
	int option = 0;
	while ((option = getopt_long_only(argc, argv, "s:n:b:p:l:o:w:g:h", long_opts, &long_index)) != -1)
	{
		switch (option) 
		{
			//-----------------//
			//     Problem     //
			//-----------------//

			case OPT_Geometry:
				args.Problem.GeoCode = optarg;
				if (!Utils::IsPredefinedGeometry(args.Problem.GeoCode))
					args.Discretization.Mesher = "gmsh";
				break;
			case OPT_TestCase:
			case OPT_RightHandSide:
				args.Problem.TestCaseCode = optarg;
				break;
			case OPT_BoundaryConditions:
				args.Problem.BCCode = optarg;
				break;
			case OPT_HeterogeneityRatio: 
				args.Problem.HeterogeneityRatio = atof(optarg);
				break;
			case OPT_AnisotropyRatio:
				args.Problem.AnisotropyRatio = atof(optarg);
				break;
			case OPT_AnisotropyAngle:
				args.Problem.AnisotropyAngle = atof(optarg) * M_PI / 180; // conversion degrees to radians
				break;

			//--------------------//
			//        Mesh        //
			//--------------------//

			case OPT_Mesh:
			{
				string meshCode = optarg;
				if (   meshCode.compare("cart") != 0
					&& meshCode.compare("cart-poly") != 0
					&& meshCode.compare("stri") != 0
					&& meshCode.compare("tri") != 0
					&& meshCode.compare("tetra") != 0
					&& meshCode.compare("stetra") != 0
					&& meshCode.compare("quad") != 0
					&& meshCode.compare("quad-poly") != 0)
					argument_error("unknown mesh code '" + meshCode + "'. Check -mesh argument.");
				args.Discretization.MeshCode = meshCode;
				break;
			}
			case OPT_Mesher:
				args.Discretization.Mesher = optarg;
				if (args.Discretization.Mesher.compare("inhouse") != 0 && args.Discretization.Mesher.compare("gmsh"))
					argument_error("unknown mesher '" + args.Discretization.Mesher + "'. Check -mesher argument.");
				break;
			case 'n':
			case OPT_Nx:
				args.Discretization.N = stoul(optarg, nullptr, 0);
				break;
			case OPT_Ny:
				args.Discretization.Ny = stoul(optarg, nullptr, 0);
				break;
			case OPT_Nz:
				args.Discretization.Nz = stoul(optarg, nullptr, 0);
				break;
			case OPT_Stretch:
				args.Discretization.Stretch = atof(optarg);
				break;

			//--------------------//
			//   Discretization   //
			//--------------------//

			case OPT_Discretization:
			{
				string discretization = optarg;
				if (discretization.compare("dg") != 0 && discretization.compare("hho") != 0)
					argument_error("unknown discretization '" + discretization + "'. Check -discr argument.");
				args.Discretization.Method = discretization;
				break;
			}
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
				if (   initialGuessCode.compare("0") != 0 
					&& initialGuessCode.compare("1") != 0 
					&& initialGuessCode.compare("rand") != 0
					&& initialGuessCode.compare("smooth") != 0)
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
			case OPT_Relaxation:
				args.Solver.RelaxationParameter = atof(optarg);
				break;

			//---------------//
			//   Multigrid   //
			//---------------//

			case OPT_MGCycle:
			{
				defaultCycle = false;
				string s(optarg);
				regex pattern("^([vwkVWK]),([[:digit:]]+),([[:digit:]]+)(,([*+-])([[:digit:]]+))?$");
				smatch matches;

				if (std::regex_search(s, matches, pattern))
				{
					char cycleLetter = matches.str(1)[0];
					if (cycleLetter == 'v' || cycleLetter == 'V')
					{
						args.Solver.MG.CycleLetter = 'V';
						args.Solver.MG.WLoops = 1;
					}
					else if (cycleLetter == 'w' || cycleLetter == 'W')
					{
						args.Solver.MG.CycleLetter = 'W';
						args.Solver.MG.WLoops = 2;
					}
					else if (cycleLetter == 'k' || cycleLetter == 'K')
					{
						args.Solver.MG.CycleLetter = 'K';
						args.Solver.MG.WLoops = 1;
					}
					args.Solver.MG.PreSmoothingIterations = stoi(matches.str(2));
					args.Solver.MG.PostSmoothingIterations = stoi(matches.str(3));
					if (!matches.str(4).empty())
					{
						args.Solver.MG.CoarseLevelChangeSmoothingOperator = matches.str(5)[0];
						args.Solver.MG.CoarseLevelChangeSmoothingCoeff = stoi(matches.str(6));
					}
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
			case OPT_Weight:
			{
				string weightCode = optarg;
				if (weightCode.compare("k") != 0 && weightCode.compare("a") != 0)
					argument_error("unknown weight code '" + weightCode + "'. Check -weight argument.");
				args.Solver.MG.WeightCode = weightCode;
				break;
			}
			case OPT_ProlongationCode:
			{
				int prolongationCode = atoi(optarg);
				if (prolongationCode < 1 || prolongationCode > 8)
					argument_error("unknown prolongation code. Check -prolong argument.");
				args.Solver.MG.ProlongationCode = static_cast<Prolongation>(prolongationCode);
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
				if (coarseningStgyCode.compare("s") == 0)
					args.Solver.MG.CoarseningStgy = CoarseningStrategy::StandardCoarsening;
				else if (coarseningStgyCode.compare("a") == 0)
					args.Solver.MG.CoarseningStgy = CoarseningStrategy::AgglomerationCoarsening;
				else if (coarseningStgyCode.compare("l") == 0)
					args.Solver.MG.CoarseningStgy = CoarseningStrategy::AgglomerationCoarseningByMostCoplanarFaces;
				else if (coarseningStgyCode.compare("c") == 0)
					args.Solver.MG.CoarseningStgy = CoarseningStrategy::AgglomerationCoarseningByClosestCenter;
				else if (coarseningStgyCode.compare("g") == 0)
					args.Solver.MG.CoarseningStgy = CoarseningStrategy::AgglomerationCoarseningByClosestFace;
				else if (coarseningStgyCode.compare("i") == 0)
					args.Solver.MG.CoarseningStgy = CoarseningStrategy::AgglomerationCoarseningByLargestInterface;
				else if (coarseningStgyCode.compare("p") == 0)
					args.Solver.MG.CoarseningStgy = CoarseningStrategy::AgglomerationCoarseningBySeedPoints;
				else if (coarseningStgyCode.compare("n") == 0)
					args.Solver.MG.CoarseningStgy = CoarseningStrategy::AgglomerationCoarseningByFaceNeighbours;
				else if (coarseningStgyCode.compare("v") == 0)
					args.Solver.MG.CoarseningStgy = CoarseningStrategy::AgglomerationCoarseningByVertexNeighbours;
				else if (coarseningStgyCode.compare("m") == 0)
					args.Solver.MG.CoarseningStgy = CoarseningStrategy::IndependentRemeshing;
				else if (coarseningStgyCode.compare("f") == 0)
					args.Solver.MG.CoarseningStgy = CoarseningStrategy::FaceCoarsening;
				else if (coarseningStgyCode.compare("r") == 0)
					args.Solver.MG.CoarseningStgy = CoarseningStrategy::GMSHSplittingRefinement;
				else if (coarseningStgyCode.compare("b") == 0)
					args.Solver.MG.CoarseningStgy = CoarseningStrategy::BeyRefinement;
				else
					argument_error("unknown coarsening strategy code '" + coarseningStgyCode + "'. Check -cs argument.");
				break;
			}
			case OPT_CoarseN:
				args.Solver.MG.CoarseN = stoul(optarg, nullptr, 0);
				break;

			//----------------//
			//      Misc      //
			//----------------//

			case 'h':
				print_usage();
				exit(EXIT_SUCCESS);
				break;
			case OPT_Threads:
				BaseParallelLoop::SetDefaultNThreads(atoi(optarg));
				break;
			case OPT_Export:
			{
				vector<string> exports = Utils::Explode(optarg, ',');
				for (string code : exports)
				{
					if (code.compare("lsys") == 0)
						args.Actions.ExportLinearSystem = true;
					else if (code.compare("amat") == 0)
						args.Actions.ExportAssemblyTermMatrices = true;
					else if (code.compare("mesh") == 0)
						args.Actions.ExportMeshToMatlab = true;
					else if (code.compare("solvect") == 0)
						args.Actions.ExportSolutionVectors = true;
					else if (code.compare("solgmsh") == 0)
						args.Actions.ExportSolutionToGMSH = true;
					else if (code.compare("mg") == 0)
						args.Actions.ExportMultigridComponents = true;
					else
						argument_error("unknown export option '" + code + "'.");
				}
				break;
			}
			case OPT_DoNotSolve:
				args.Actions.SolveLinearSystem = false;
				break;
			case OPT_NoCache:
				args.Actions.UseCache = false;
				break;
			case OPT_UnitTests:
				args.Actions.UnitTests = true;
				break;
			case OPT_GMSHLog:
				args.Actions.GMSHLogEnabled = true;
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

	//------------------------------------------//
	//             Problem dimension            //
	//------------------------------------------//

	if (args.Problem.Dimension == -1)
	{
		if (args.Problem.GeoCode.compare("segment") == 0)
			args.Problem.Dimension = 1;
		else if (args.Problem.GeoCode.compare("square") == 0 || args.Problem.GeoCode.compare("square4quadrants") == 0)
			args.Problem.Dimension = 2;
		else if (args.Problem.GeoCode.compare("cube") == 0)
			args.Problem.Dimension = 3;
#ifdef GMSH_ENABLED
		else
		{
			Mesh<2>::SetDirectories();
			args.Problem.Dimension = GMSHMesh<2>::GetDimension(args.Problem.GeoCode);
		}
#else
		else
			argument_error("Unknown geometry.");
#endif
	}

	GMSHMesh<3>::GMSHLogEnabled = true;

	// Test case
	if (args.Problem.TestCaseCode.compare("") == 0)
	{
		if (Utils::IsPredefinedGeometry(args.Problem.GeoCode))
			args.Problem.TestCaseCode = "default";
		else
			args.Problem.TestCaseCode = FileSystem::FileNameWithoutExtension(args.Problem.GeoCode);
	}

	//------------------------------------------//
	//               Heterogeneity              //
	//------------------------------------------//

	if (args.Problem.GeoCode.compare("square") == 0 && args.Problem.HeterogeneityRatio != 1)
		Utils::Warning("The geometry 'square' has only one physical part: -heterog argument is ignored. Use 'square4quadrants' instead to run an heterogeneous problem.");
	
	if (args.Problem.Dimension > 1 && args.Discretization.Method.compare("dg") == 0 && args.Discretization.PolyDegree == 0)
		argument_error("In 2D/3D, DG is not a convergent scheme for p = 0.");

	if (args.Discretization.Method.compare("dg") == 0 && args.Problem.BCCode.compare("d") != 0)
		argument_error("In DG, only Dirichlet conditions are implemented.");

	if (args.Discretization.Method.compare("dg") == 0 && args.Problem.AnisotropyRatio != 1)
		argument_error("In DG, anisotropy is not implemented.");

	if (args.Problem.Dimension == 1 && args.Discretization.Method.compare("hho") == 0 && args.Discretization.PolyDegree != 1)
		argument_error("HHO in 1D only exists for p = 1.");

	if (args.Discretization.Method.compare("hho") == 0 && args.Discretization.PolyDegree == 0)
		argument_error("HHO does not exist with p = 0. Linear approximation at least (p >= 1).");

	//------------------------------------------//
	//                   Mesh                   //
	//------------------------------------------//

#ifndef GMSH_ENABLED
	if (args.Discretization.Mesher.compare("gmsh") == 0)
		argument_error("GMSH is disabled. Recompile with the cmake option -DENABLE_GMSH=ON to use GMSH meshes, or choose another argument for -mesh.");
#endif // GMSH_ENABLED

	if (args.Discretization.MeshCode.compare("default") == 0)
	{
		if (args.Problem.Dimension == 1)
			args.Discretization.MeshCode = "cart";
		else if (args.Problem.Dimension == 2)
			args.Discretization.MeshCode = args.Discretization.Mesher.compare("inhouse") == 0 ? "stri" : "tri";
		else if (args.Problem.Dimension == 3)
			args.Discretization.MeshCode = args.Discretization.Mesher.compare("inhouse") == 0 ? "stetra" : "tetra";
	}

	if ((args.Discretization.MeshCode.compare("tri") == 0 || args.Discretization.MeshCode.compare("stri") == 0) && args.Problem.Dimension != 2)
		argument_error("Triangular mesh in only available in 2D.");

	if (args.Discretization.MeshCode.compare("quad") == 0 && args.Problem.Dimension != 2)
		argument_error("Quadrilateral mesh in only available in 2D.");

	if ((args.Discretization.MeshCode.compare("tetra") == 0 || args.Discretization.MeshCode.compare("stetra") == 0) && args.Problem.Dimension != 3)
		argument_error("Tetrahedral mesh in only available in 3D.");

	//------------------------------------------//
	//                  Solver                  //
	//------------------------------------------//

	if (args.Solver.SolverCode.compare("default") == 0)
	{
		if (args.Discretization.Method.compare("hho") == 0 && args.Discretization.StaticCondensation && args.Problem.Dimension > 1)
		{
			args.Solver.SolverCode = "mg";
			if ((args.Solver.MG.ProlongationCode == Prolongation::Wildey || args.Solver.MG.ProlongationCode == Prolongation::FaceInject) && !args.Solver.MG.UseGalerkinOperator)
			{
				Utils::Warning("The multigrid with prolongation code " + to_string((unsigned)args.Solver.MG.ProlongationCode) + " requires the Galerkin operator. Option -g 0 ignored.");
				args.Solver.MG.UseGalerkinOperator = true;
			}
		}
		else if ((args.Problem.Dimension == 2 && args.Discretization.N < 64) || (args.Problem.Dimension == 3 && args.Discretization.N < 16))
			args.Solver.SolverCode = "lu";
		else
			args.Solver.SolverCode = "eigencg";
	}

#ifndef AGMG_ENABLED
	if (args.Solver.SolverCode.compare("agmg") == 0)
		argument_error("AGMG is disabled. Recompile with the cmake option -DENABLE_AGMG=ON, or choose another solver.");
#endif // AGMG_ENABLED

	//------------------------------------------//
	//                Multigrid                 //
	//------------------------------------------//

	if (args.Solver.SolverCode.compare("mg") == 0 || args.Solver.SolverCode.compare("pcgmg") == 0)
	{
		if (args.Discretization.Method.compare("dg") == 0)
			argument_error("Multigrid only applicable on HHO discretization.");

		if (!args.Discretization.StaticCondensation)
			argument_error("Multigrid only applicable if the static condensation is enabled.");

		if (args.Solver.MG.ProlongationCode == Prolongation::Wildey && !args.Solver.MG.UseGalerkinOperator)
			argument_error("To use the prolongationCode " + to_string((unsigned)Prolongation::Wildey) + ", you must also use the Galerkin operator. To do so, add option -g 1.");

		if (args.Solver.MG.CoarseningStgy == CoarseningStrategy::FaceCoarsening && !args.Solver.MG.UseGalerkinOperator)
			argument_error("To use the face coarsening, you must also use the Galerkin operator. To do so, add option -g 1.");

		if (args.Solver.MG.CoarseningStgy == CoarseningStrategy::IndependentRemeshing && 
			args.Solver.MG.ProlongationCode != Prolongation::CellInterp_L2proj_Trace &&
			args.Solver.MG.ProlongationCode != Prolongation::CellInterp_ApproxL2proj_Trace &&
			args.Solver.MG.ProlongationCode != Prolongation::Default)
			argument_error("The coarsening by independent remeshing is only applicable with the non-nested versions of the multigrid (-prolong " + to_string((unsigned)Prolongation::CellInterp_L2proj_Trace) + " or " + to_string((unsigned)Prolongation::CellInterp_ApproxL2proj_Trace) + ").");

		if (args.Solver.MG.CoarseningStgy == CoarseningStrategy::None)
		{
			if (args.Problem.Dimension < 3)
			{
				if (args.Solver.MG.ProlongationCode == Prolongation::Default)
				{
					args.Solver.MG.CoarseningStgy = CoarseningStrategy::AgglomerationCoarseningByFaceNeighbours;
					args.Solver.MG.ProlongationCode = Prolongation::CellInterp_ApproxL2proj_Trace;
				}
				else if (Utils::RequiresNestedHierarchy(args.Solver.MG.ProlongationCode))
					args.Solver.MG.CoarseningStgy = CoarseningStrategy::GMSHSplittingRefinement;
				else
					args.Solver.MG.CoarseningStgy = CoarseningStrategy::AgglomerationCoarseningByFaceNeighbours;
			}
			else
			{
				if (args.Discretization.Mesher.compare("inhouse") == 0 && args.Discretization.MeshCode.compare("tetra") == 0)
					args.Solver.MG.CoarseningStgy = CoarseningStrategy::BeyRefinement;
				else
					args.Solver.MG.CoarseningStgy = CoarseningStrategy::IndependentRemeshing;
			}
		}

		if (args.Solver.MG.ProlongationCode == Prolongation::Default)
		{
			if (Utils::BuildsNestedMeshHierarchy(args.Solver.MG.CoarseningStgy))
				args.Solver.MG.ProlongationCode = Prolongation::CellInterp_Trace;
			else
				args.Solver.MG.ProlongationCode = Prolongation::CellInterp_ApproxL2proj_Trace;
		}

		if (defaultCycle)
		{
			args.Solver.MG.PreSmoothingIterations = 0;
			if (args.Problem.Dimension < 3)
				args.Solver.MG.PostSmoothingIterations = 3;
			else
				args.Solver.MG.PostSmoothingIterations = Utils::IsRefinementStrategy(args.Solver.MG.CoarseningStgy) ? 10 : 6;
		}
	}
	
	//------------------------------------------//
	//             Launch program               //
	//------------------------------------------//

	Program* program = nullptr;
	if (args.Problem.Dimension == 1)
		program = new ProgramDim<1>();
	else if (args.Problem.Dimension == 2)
		program = new ProgramDim<2>();
	else if (args.Problem.Dimension == 3)
		program = new ProgramDim<3>();

	try
	{
		program->Start(args);
	}
	catch (exception* e)
	{
		string error(e->what());
		Utils::FatalError("Unhandled exception: " + error);
	}

	delete program;

	cout << "----------------- SUCCESSFUL TERMINATION ----------------" << endl;
    return EXIT_SUCCESS;
}
