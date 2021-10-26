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
	cout << "      Test case code. By default, it is defined as the geometry." << endl;
	cout << "      The test case defines the source function of the problem, the diffusion field, as well as available predefined boundary conditions." << endl;
	cout << "      By default, the test case is defined as the geometry." << endl;
	cout << "      When the geometry given is a GMSH file, a test case with same code as the file name should be defined in the source code. Alternatively, you can use" << endl;
	cout << "               default - discontinuous source function, homogeneous Dirichlet conditions" << endl;
	cout << "      Additional specific test cases:" << endl;
	cout << "               heterog - ('segment' geometry only) heterogeneous diffusion-specific analytical solution" << endl;
	cout << "               kellogg - ('square4quadrants' geometry only) heterogeneous diffusion-specific analytical solution (known benchmark)" << endl;
	cout << endl;
	cout << "-source CODE" << endl;
	cout << "      Code allowing to change the default source function defined in the test case." << endl;
	cout << "      For some, it also determines the analytical solution in the homogeneous isotropic case so that the L2 error can be computed." << endl;
	cout << "      The following test cases are predefined and available for the simple geometries listed in the -geo argument." << endl;
	cout << "               sine    - the source function and the analytical solution are a sine functions" << endl;
	cout << "               poly    - the source function is constant, the analytical solution is a polynomial of total degree 2*d" << endl;
	cout << "               zero    - the source function and the analytical solution are 0" << endl;
	cout << "               one     - the source function is 0, the analytical solution is 1" << endl;
	cout << "               x       - the source function is 0, the analytical solution is x" << endl;
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
	cout << "-e-basis CODE" << endl;
	cout << "      Polynomial basis for the elements." << endl;
	cout << "               monomials" << endl;
	cout << "               legendre" << endl;
	cout << "               nlegendre (normalized Legendre)" << endl;
	cout << "               bernstein" << endl;
	cout << "               hemker" << endl;
	cout << endl;
	cout << "-f-basis CODE" << endl;
	cout << "      Polynomial basis for the faces. Same values as for -e-basis." << endl;
	cout << endl;
	cout << "-e-ogb NUM" << endl;
	cout << "      Orthogonalization of the local bases against each element. Default: 1." << endl;
	cout << "               0 - no orthogonalization" << endl;
	cout << "               1 - orthogonalization, without normalization" << endl;
	cout << "               2 - double orthogonalization, without normalization" << endl;
	cout << "               3 - orthonormalization" << endl;
	cout << "               4 - orthonormalization with double orthogonalization" << endl;
	cout << endl;
	cout << "-f-ogb NUM" << endl;
	cout << "      Orthogonalization of the local bases against each face. Same value as -e-ogb." << endl;
	cout << endl;
	cout << "-p NUM" << endl;
	cout << "      Polynomial degree of approximation (default: 1). In HHO, k = p-1." << endl;
	cout << endl;
	cout << "-k NUM" << endl;
	cout << "      Polynomial degree on the faces for the HHO scheme. Note that setting this parameter does the same as setting p (p = k+1)." << endl;
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
	cout << "      Linear solver for the solution of the system." << endl;
	cout << "      - Direct methods:" << endl;
	cout << "              lu           - LU factorization (Eigen library)" << endl;
	cout << "              ch           - Cholesky factorization (Eigen library)" << endl;
	cout << "      - Fixed-point iterative methods:" << endl;
	cout << "              [b]j         - [Block] Jacobi. Use argument -relax to change the relaxation parameter." << endl;
	cout << "              [b]j23       - [Block] Jacobi with 2/3 as relaxation parameter." << endl;
	cout << "              [r|s][b]gs   - [Reverse|Symmetric][Block] Gauss-Seidel" << endl;
	cout << "              [r|s][b]sor  - [Reverse|Symmetric][Block] SOR. Use argument -relax to change the relaxation parameter." << endl;
	cout << "        For the block solvers, the block size is set to the number of DOFs per cell (DG) or face (HHO)." << endl;
	cout << "      - Krylov methods:" << endl;
	cout << "              cg           - Conjugate Gradient" << endl;
	cout << "              eigencg      - Conjugate Gradient (Eigen library) with diagonal preconditioner" << endl;
	cout << "              fcg          - Flexible Conjugate Gradient (truncation-restart: FCG(1))" << endl;
	cout << "        'cg' and 'fcg' can be suffixed with another solver to use it as a preconditioner. Ex: fcgmg" << endl;
	cout << "      - Multigrid methods:" << endl;
	cout << "              agmg         - Yvan Notay's AGMG solver" << endl;
	cout << "              mg           - Custom multigrid for HHO" << endl;
	cout << "              p_mg         - p-Multigrid to be used on top of mg" << endl;
	cout << "              uamg         - Uncondensed AMG (for hybrid discretizations with static condensation)" << endl;
	cout << "              aggregamg    - In-house implementation of AGMG, without the outer Krylov iteration. Meant to be used with K-cycle and with FCG." << endl;
	cout << "              hoaggregamg  - Algebraic p-multigrid on top of aggregamg" << endl;
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
	cout << "      Relaxation parameter used for Jacobi, Gauss-Seidel, SOR and their derived versions (default: 1)." << endl;
	cout << endl;
	cout << "-block-size NUM" << endl;
	cout << "      Forces a block size for the block verions of Jacobi, Gauss-Seidel, SOR." << endl;
	cout << "      By default, the value is adapted to the space dimension and the polynomial order of the discretization." << endl;
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
	cout << "      Pre-smoother,post-smoother: \"bgs,rbgs\" for example." << endl;
	cout << "              j    - Jacobi" << endl;
	cout << "              gs   - Gauss-Seidel" << endl;
	cout << "              rgs  - Reverse Gauss-Seidel" << endl;
	cout << "              sgs  - Symmetric Gauss-Seidel" << endl;
	cout << "              bj   - Block Jacobi: the block size is set to the number of DOFs per face" << endl;
	cout << "              bj23 - Block Jacobi with a damping factor of 2/3" << endl;
	cout << "              bgs  - Block Gauss-Seidel: the block size is set to the number of DOFs per face" << endl;
	cout << "              rbgs - Reverse Block Gauss-Seidel" << endl;
	cout << "              sbgs - Symmetric Block Gauss-Seidel" << endl;
	cout << endl;
	cout << "-hp-cs CODE" << endl;
	cout << "      hp-coarsening strategy to build the coarse levels in the case of high-order." << endl;
	cout << "              h    - h only" << endl;
	cout << "              p    - p only" << endl;
	cout << "              p_h  - p, then h" << endl;
	cout << "              h_p  - h, then p" << endl;
	cout << "              hp_h - simultaneous hp, then h" << endl;
	cout << "              hp_p - simultaneous hp, then p" << endl;
	cout << "              p_hp - p, then simultaneous hp" << endl;
	cout << "              alt  - alternation of p and h, starting with p" << endl;
	cout << endl;
	cout << "-cs CODE" << endl;
	cout << "      Mesh coarsening strategy." << endl;
	cout << "      'mg':" << endl;
	cout << "              s    - Standard coarsening (only for structured meshes built by the in-house mesher)" << endl;
	cout << "              r    - Fine meshes obtained by structured refinement of the coarse mesh using GMSH's splitting method" << endl;
	cout << "              b    - Fine meshes obtained by Bey's tetrahedral refinements of coarse meshes" << endl;
	cout << "              m    - Independant remeshing by GMSH with double the mesh size (non-nested!)" << endl;
	cout << "              n    - Agglomeration coarsening (agglomerate all face neighbours)" << endl;
	cout << "              mn   - Multiple agglomeration coarsening (use -coarsening-factor)" << endl;
	cout << "              dpa  - Double pairwise aggregation" << endl;
	cout << "              mpa  - Multiple pairwise aggregation (use -coarsening-factor)" << endl;
	cout << "              vr   - (Experimental) Agglomeration coarsening by vertex removal" << endl;
	cout << "              mcf  - (Experimental) agglomeration coarsening by most collinear/coplanar faces" << endl;
	cout << "              cc   - (Experimental) agglomeration coarsening by closest center" << endl;
	cout << "              clf  - (Experimental) agglomeration coarsening by closest face" << endl;
	cout << "              li   - (Experimental) agglomeration coarsening by largest interface" << endl;
	cout << "              seed - (Experimental) agglomeration coarsening by seed points" << endl;
	cout << "              vn   - (Experimental) agglomeration coarsening by vertex neighbours" << endl;
	cout << "              f    - (Experimental) Face coarsening: the faces are coarsened and all kept on the coarse skeleton. Requires -g 1." << endl;
	cout << endl;
	cout << "-p-cs CODE" << endl;
	cout << "      p-coarsening strategy:" << endl;
	cout << "              -1    - decrement the degree by 1" << endl;
	cout << "              -2    - decrement the degree by 2 (default)" << endl;
	cout << "              /2    - divide the degree by 2" << endl;
	cout << "              =0    - go directly to the low order" << endl;
	cout << endl;
	cout << "-fcs CODE" << endl;
	cout << "      Face coarsening (default: c)." << endl;
	cout << "              c   - Collapse interfaces made of multiple faces" << endl;
	cout << "              n   - None (fine faces are kept as is on the coarse mesh)" << endl;
	cout << "              i   - Collapse interfaces made of multiple faces and try to aggregate interior faces to the boundary ones" << endl;
	//cout << "              z   - Aggregate all faces using A_F_F" << endl;
	cout << endl;
	cout << "-num-meshes NUM" << endl;
	cout << "      Fixes the number of meshes in the mesh hierarchy." << endl;
	cout << "      By default, the program tries to coarsen the mesh until the linear system reaches the size set by the argument -coarse-size." << endl;
	cout << endl;
	cout << "-coarse-solver CODE" << endl;
	cout << "      Any solver code (see -solver). Only purely algebraic solvers are allowed. Default: 'ch' (Cholesky factorization)." << endl;
	cout << endl;
	cout << "-bfc CODE" << endl;
	cout << "      Face collapsing method used at the domain boundaries or physical parts boundaries. Requires -cs n." << endl;
	cout << "              d   - Disabled" << endl;
	cout << "              c   - Collinear only" << endl;
	cout << "              p   - By pairs" << endl;
	cout << "              m   - Maximum" << endl;
	cout << endl;
	cout << "-rcm CODE" << endl;
	cout << "      Re-entrant corner management. Requires -cs n." << endl;
	cout << "              d   - Disabled" << endl;
	cout << "              f   - Agglomerate elements at re-entrant corners first" << endl;
	cout << endl;
	cout << "-coarse-n NUM" << endl;
	cout << "      If a refinement strategy is used, sets the mesh size of the coarse mesh." << endl;
	cout << "      The chosen value will set the variable N defined at the beginning of the GMSH .geo file." << endl;
	cout << endl;
	cout << "-coarsening-factor NUM" << endl;
	cout << "      Requested coarsening factor for the coarsening strategies compatible." << endl;
	cout << "              Independent remeshing         (-cs r  ): coarsening factor H/h, default 2" << endl;
	cout << "              Multiple pairwise aggregation (-cs mpa): coarsening factor #FineVariables/#CoarseVariables, default 3.5" << endl;
	cout << endl;
	cout << "-prolong NUM" << endl;
	cout << "      How the prolongation operator is built." << endl;
	cout << "      Values for the geometric multigrid for HHO:" << endl;
	cout << "              " << (unsigned)GMG_H_Prolongation::CellInterp_Trace << "  - ";
	cout <<                    "Step 1: Interpolation from coarse faces to coarse cells" << endl;
	cout << "                   Step 2: Trace on the fine faces" << endl;
	cout << "              " << (unsigned)GMG_H_Prolongation::CellInterp_InjectAndTrace << "  - ";
	cout <<                    "Step 1: Interpolation from coarse faces to coarse cells" << endl;
	cout << "                   Step 2: On faces present on both fine and coarse meshes, we keep the polynomials identical." << endl;
	cout << "                           On faces interior to coarse elements, trace of the cell polynomials." << endl;
	cout << "              " << (unsigned)GMG_H_Prolongation::CellInterp_Inject_Adjoint << "  - ";
	cout <<                    "Step 1: Interpolation from coarse faces to coarse cells" << endl;
	cout << "                   Step 2: Canonical injection from coarse to fine cells" << endl;
	cout << "                   Step 3: Adjoint of the cell interpolation on the fine mesh" << endl;
	cout << "              " << (unsigned)GMG_H_Prolongation::Wildey << "  - ";
	cout <<                    "Algorithm from Wildey et al.: the coarse level is built by static condensation of the fine faces interior to coarse elements." << endl;
	cout << "                   The prolongation solves those condensed unknowns." << endl;
	cout << "                   To reproduce Wildey et al.'s algorithm, this option should be used with '-g 1 -cycle V,1,1,*2'." << endl;
	cout << "              " << (unsigned)GMG_H_Prolongation::FaceInject << "  - ";
	cout <<                    "Canonical injection from coarse faces to fine faces (implemented to be used with option -cs f)." << endl;
	cout << "              " << (unsigned)GMG_H_Prolongation::CellInterp_Inject_Trace << "  - ";
	cout <<                    "Same as 1, but another implementation:" << endl;
	cout << "                   Step 1: Interpolation from coarse faces to coarse cells" << endl;
	cout << "                   Step 2: Canonical injection from coarse to fine cells" << endl;
	cout << "                   Step 3: Trace on the fine faces" << endl;
	cout << "              " << (unsigned)GMG_H_Prolongation::CellInterp_ExactL2proj_Trace << "  - ";
	cout <<                    "Non-nested variant of " << (unsigned)GMG_H_Prolongation::CellInterp_Trace << " and " << (unsigned)GMG_H_Prolongation::CellInterp_Inject_Trace << ":" << endl;
	cout << "                   Step 1: Interpolation from coarse faces to coarse cells" << endl;
	cout << "                   Step 2: Exact L2-projection onto the fine cells" << endl;
	cout << "                   Step 3: Trace on the fine faces" << endl;
	cout << "              " << (unsigned)GMG_H_Prolongation::CellInterp_ApproxL2proj_Trace << "  - ";
	cout <<                    "Variant of " << (unsigned)GMG_H_Prolongation::CellInterp_ExactL2proj_Trace << " where the L2-projection is not computed exactly but has the same approximation properties." << endl;
	cout << "              " << (unsigned)GMG_H_Prolongation::CellInterp_FinerApproxL2proj_Trace << "  - ";
	cout <<                    "Variant of " << (unsigned)GMG_H_Prolongation::CellInterp_ApproxL2proj_Trace << " where the L2-projection is better approximated (by subtriangulation of the elements)." << endl;
	cout << "      Values for UncondensedAMG:" << endl;
	cout << "              " << (unsigned)UAMGProlongation::ChainedCoarseningProlongations << "  - ";
	cout <<                    "chained coarsening prolongations (use -coarsening-prolong to select)" << endl;
	cout << "              " << (unsigned)UAMGProlongation::ReconstructionTrace << "  - ";
	cout <<                    "cell-reconstruction + injection + trace" << endl;
	cout << "              " << (unsigned)UAMGProlongation::FaceProlongation << "  - ";
	cout <<                    "face injection" << endl;
	cout << "              " << (unsigned)UAMGProlongation::ReconstructTraceOrInject << "  - ";
	cout <<                    "cell-reconstruction + trace for interior faces, injection for boundary faces" << endl;
	cout << "              " << (unsigned)UAMGProlongation::ReconstructSmoothedTraceOrInject << "  - ";
	cout <<                    "cell-reconstruction + trace + smoothing for interior faces, injection for boundary faces" << endl;
	cout << endl;
	cout << "-coarsening-prolong NUM" << endl;
	cout << "      Intermediate prolongation operator used to build the coarse levels in 'uamg'. Same possble values as -prolong." << endl;
	cout << endl;
	cout << "-face-prolong NUM" << endl;
	cout << "      Intermediate face prolongation operator used in 'uamg'." << endl;
	cout << "              " << (unsigned)UAMGFaceProlongation::BoundaryAggregatesInteriorAverage << "  - ";
	cout <<                     "1 non-zero for aggregated faces, average for interior non-aggregated faces" << endl;
	cout << "              " << (unsigned)UAMGFaceProlongation::BoundaryAggregatesInteriorZero << "  - ";
	cout <<                     "1 non-zero for boundary faces, nothing for interior faces" << endl;
	cout << "              " << (unsigned)UAMGFaceProlongation::FaceAggregates << "  - ";
	cout <<                     "1 non-zero per face (as in AGMG)" << endl;
	cout << endl;
	cout << "-p-prolong NUM" << endl;
	cout << "      How the p-prolongation operator for 'mg' is built." << endl;
	cout << "              " << (unsigned)GMG_P_Prolongation::Injection << "  - natural injection (default)" << endl;
	cout << "              " << (unsigned)GMG_P_Prolongation::H_Prolongation << "  - same as prolongation in h" << endl;
	cout << endl;
	cout << "-p-restrict NUM" << endl;
	cout << "      How the p-restriction operator for 'mg' is built." << endl;
	cout << "              " << (unsigned)GMG_P_Restriction::RemoveHigherOrders << "  - remove higher orders (default), which corresponds to the L2-proj. if the basis is hierarchical and orthonormalized" << endl;
	cout << "              " << (unsigned)GMG_P_Restriction::P_Transpose << "  - transpose of p-prolongation" << endl;
	cout << endl;
	cout << "-hp-config NUM" << endl;
	cout << "      Shortcut for the following set of arguments:" << endl;
	cout << "              1  -   -hp-cs h" << endl;
	cout << "              2  -   -hp-cs p_h  -p-cs -2 -p-prolong 1 -p-restrict 1" << endl;
	cout << "              3  -   -hp-cs p_h  -p-cs -2 -p-prolong 2 -p-restrict 2" << endl;
	cout << "              4  -   -hp-cs hp_h -p-cs -1" << endl;
	cout << "              5  -   -hp-cs alt  -p-cs -2" << endl;
	cout << endl;
	cout << "-subtri NUM" << endl;
	cout << "      If the approximated L2-projection is used in the multigrid, sets the number of subtriangulations of the fine elements." << endl;
	cout << endl;
	cout << "-disable-hor" << endl;
	cout << "      In the polongation of 'mg', disables the use of the higher-order reconstruction." << endl;
	cout << "      The local cell polynomials in the intermediary step are obtained by solving the local problems only." << endl;
	cout << endl;
	cout << "-disable-heterog-weight" << endl;
	cout << "      In the polongation of 'mg', disables the heterogeneous weighting." << endl;
	cout << "      Homogeneous weighting is used instead, i.e. non-weighted average." << endl;
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
	cout << "              lsys       - linear system" << endl;
	cout << "              amat       - matrices decomposed by assembly terms in separate files (consistency, stabilization, mass, etc.)" << endl;
	cout << "              mesh       - mesh to be used in Matlab" << endl;
	cout << "              solvect    - solution vector(s)" << endl;
	cout << "              solgmsh    - solution files (.pos and .msh) to be used in GMSH for visualization" << endl;
	cout << "              errgmsh    - error computed against the solution of an exact solver to be used in GMSH for visualization" << endl;
	cout << "              sourcegmsh - solution files (.pos and .msh) to be used in GMSH for visualization" << endl;
	cout << "              mg         - Multigrid components (intergrid operator matrices, coarse meshes, etc.)" << endl;
	cout << "              mgit       - Multigrid iteration vectors: for each iteration, right-hand side, solution before/after smoothing, etc." << endl;
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
#ifdef CGAL_ENABLED
	CGALWrapper::Configure();
#endif

	bool defaultCycle = true;
	bool defaultCoarseOperator = true;
	bool defaultCoarseSolver = true;

	ProgramArguments args;
	args.OutputDirectory = FileSystem::RootPath() + "/out";
	FileSystem::CreateDirectoryIfNotExist(args.OutputDirectory);

	enum {
		// Problem
		OPT_Geometry = 1000,
		OPT_TestCase,
		OPT_Source,
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
		OPT_HHO_K,
		OPT_Stabilization,
		OPT_ElemBasis,
		OPT_FaceBasis,
		OPT_OrthogonalizeElemBases,
		OPT_OrthogonalizeFaceBases,
		OPT_NoStaticCondensation,
		OPT_Penalization,
		OPT_PolySpace,
		// Solver
		OPT_InitialGuess,
		OPT_Tolerance,
		OPT_MaxIterations,
		OPT_Relaxation,
		OPT_BlockSize,
		// Multigrid
		OPT_MGCycle,
		OPT_DisableHigherOrderReconstruction,
		OPT_DisableHeterogeneousWeighting,
		OPT_MultigridProlongationCode,
		OPT_PProlongationCode,
		OPT_PRestrictionCode,
		OPT_HPConfig,
		OPT_CoarseningProlongationCode,
		OPT_FaceProlongationCode,
		OPT_CoarseMatrixSize,
		OPT_Smoothers,
		OPT_HP_CS,
		OPT_H_CS,
		OPT_P_CS,
		OPT_FaceCoarseningStrategy,
		OPT_NumberOfMeshes,
		OPT_CoarseSolver,
		OPT_BoundaryFaceCollapsing,
		OPT_ReEntrantCornerManagement,
		OPT_CoarseningFactor,
		OPT_CoarseN,
		OPT_ApproxL2ProjNSubtriangulations,
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
		 { "source", required_argument, NULL, OPT_Source },
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
		 { "k", required_argument, NULL, OPT_HHO_K },
		 { "stab", required_argument, NULL, OPT_Stabilization },
		 { "e-basis", required_argument, NULL, OPT_ElemBasis },
		 { "f-basis", required_argument, NULL, OPT_FaceBasis },
		 { "e-ogb", required_argument, NULL, OPT_OrthogonalizeElemBases },
		 { "f-ogb", required_argument, NULL, OPT_OrthogonalizeFaceBases },
		 { "no-static-cond", no_argument, NULL, OPT_NoStaticCondensation },
		 { "pen", required_argument, NULL, OPT_Penalization },
		 { "poly-space", required_argument, NULL, OPT_PolySpace },
		 // Solver
		 { "initial-guess", required_argument, NULL, OPT_InitialGuess },
		 { "tol", required_argument, NULL, OPT_Tolerance },
		 { "max-iter", required_argument, NULL, OPT_MaxIterations },
		 { "relax", required_argument, NULL, OPT_Relaxation },
		 { "block-size", required_argument, NULL, OPT_BlockSize },
		 // Multigrid
		 { "cycle", required_argument, NULL, OPT_MGCycle },
		 { "disable-hor", no_argument, NULL, OPT_DisableHigherOrderReconstruction },
		 { "disable-heterog-weight", no_argument, NULL, OPT_DisableHeterogeneousWeighting },
		 { "prolong", required_argument, NULL, OPT_MultigridProlongationCode },
		 { "p-prolong", required_argument, NULL, OPT_PProlongationCode },
		 { "p-restrict", required_argument, NULL, OPT_PRestrictionCode },
		 { "hp-config", required_argument, NULL, OPT_HPConfig },
		 { "coarsening-prolong", required_argument, NULL, OPT_CoarseningProlongationCode },
		 { "face-prolong", required_argument, NULL, OPT_FaceProlongationCode },
		 { "coarse-size", required_argument, NULL, OPT_CoarseMatrixSize },
		 { "smoothers", required_argument, NULL, OPT_Smoothers },
		 { "hp-cs", required_argument, NULL, OPT_HP_CS },
		 { "cs", required_argument, NULL, OPT_H_CS },
		 { "p-cs", required_argument, NULL, OPT_P_CS },
		 { "fcs", required_argument, NULL, OPT_FaceCoarseningStrategy },
		 { "num-meshes", required_argument, NULL, OPT_NumberOfMeshes },
		 { "bfc", required_argument, NULL, OPT_BoundaryFaceCollapsing },
		 { "coarse-solver", required_argument, NULL, OPT_CoarseSolver },
		 { "rcm", required_argument, NULL, OPT_ReEntrantCornerManagement },
		 { "coarsening-factor", required_argument, NULL, OPT_CoarseningFactor },
		 { "coarse-n", required_argument, NULL, OPT_CoarseN },
		 { "subtri", required_argument, NULL, OPT_ApproxL2ProjNSubtriangulations },
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
				args.Problem.TestCaseCode = optarg;
				break;
			case OPT_Source:
				args.Problem.SourceCode = optarg;
				break;
			case OPT_BoundaryConditions:
				args.Problem.BCCode = optarg;
				break;
			case OPT_HeterogeneityRatio: 
				args.Problem.HeterogeneityRatio = atof(optarg);
				if (args.Problem.HeterogeneityRatio <= 0)
					argument_error("heterogeneity ratio must be > 0. Check -heterog argument.");
				break;
			case OPT_AnisotropyRatio:
				args.Problem.AnisotropyRatio = atof(optarg);
				if (args.Problem.AnisotropyRatio <= 0)
					argument_error("anisotropy ratio must be > 0. Check -aniso argument.");
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
			case OPT_ElemBasis: 
			{
				string basisCode = optarg;
				if (basisCode.compare("monomials") != 0 && basisCode.compare("legendre") != 0 && basisCode.compare("nlegendre") != 0 && basisCode.compare("bernstein") != 0 && basisCode.compare("hemker") != 0)
					argument_error("unknown polynomial basis '" + basisCode + "'. Check -e-basis argument.");
				args.Discretization.ElemBasisCode = basisCode;
				break;
			}
			case OPT_FaceBasis:
			{
				string basisCode = optarg;
				if (basisCode.compare("monomials") != 0 && basisCode.compare("legendre") != 0 && basisCode.compare("nlegendre") != 0 && basisCode.compare("bernstein") != 0 && basisCode.compare("hemker") != 0)
					argument_error("unknown polynomial basis '" + basisCode + "'. Check -f-basis argument.");
				args.Discretization.FaceBasisCode = basisCode;
				break;
			}
			case OPT_OrthogonalizeElemBases:
			{
				int i = atoi(optarg);
				if (i < 0 || i > 4)
					argument_error("check -e-ogb argument. Accepted values: 0, 1, 2, 3, 4.");
				else
					args.Discretization.OrthogonalizeElemBasesCode = i;
				break;
			}
			case OPT_OrthogonalizeFaceBases:
			{
				int i = atoi(optarg);
				if (i < 0 || i > 4)
					argument_error("check -f-ogb argument. Accepted values: 0, 1, 2, 3, 4.");
				else
					args.Discretization.OrthogonalizeFaceBasesCode = i;
				break;
			}
			case 'p': 
				args.Discretization.PolyDegree = atoi(optarg);
				break;
			case OPT_HHO_K:
				args.Discretization.PolyDegree = atoi(optarg) + 1;
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
			case OPT_BlockSize:
				args.Solver.BlockSize = atoi(optarg);
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
			case OPT_DisableHigherOrderReconstruction:
				args.Solver.MG.UseHigherOrderReconstruction = false;
				break;
			case OPT_DisableHeterogeneousWeighting:
				args.Solver.MG.UseHeterogeneousWeighting = false;
				break;
			case OPT_MultigridProlongationCode:
			{
				int prolongationCode = atoi(optarg);
				if (prolongationCode < 1 || prolongationCode > 9)
					argument_error("unknown prolongation code. Check -prolong argument.");
				args.Solver.MG.ProlongationCode = prolongationCode;
				break;
			}
			case OPT_PProlongationCode:
			{
				int prolongationCode = atoi(optarg);
				if (prolongationCode < 1 || prolongationCode > 2)
					argument_error("unknown p-prolongation code. Check -p-prolong argument.");
				args.Solver.MG.GMG_P_Prolong = static_cast<GMG_P_Prolongation>(prolongationCode);
				break;
			}
			case OPT_PRestrictionCode:
			{
				int restrictionCode = atoi(optarg);
				if (restrictionCode < 1 || restrictionCode > 2)
					argument_error("unknown p-restriction code. Check -p-restrict argument.");
				args.Solver.MG.GMG_P_Restrict = static_cast<GMG_P_Restriction>(restrictionCode);
				break;
			}
			case OPT_CoarseningProlongationCode:
			{
				int prolongationCode = atoi(optarg);
				if (prolongationCode < 1 || prolongationCode > 9)
					argument_error("unknown prolongation code. Check -coarsening-prolong argument.");
				args.Solver.MG.CoarseningProlongationCode = prolongationCode;
				break;
			}
			case OPT_FaceProlongationCode:
			{
				int prolongationCode = atoi(optarg);
				if (prolongationCode < 1 || prolongationCode > 3)
					argument_error("unknown face prolongation code. Check -face-prolong argument.");
				args.Solver.MG.FaceProlongationCode = prolongationCode;
				break;
			}
			case OPT_HPConfig:
			{
				int code = atoi(optarg);
				if (code == 1)
					args.Solver.MG.HP_CS = HP_CoarsStgy::H_only;
				else if (code == 2)
				{
					args.Solver.MG.HP_CS = HP_CoarsStgy::P_then_H;
					args.Solver.MG.P_CS = P_CoarsStgy::Minus2;
					args.Solver.MG.GMG_P_Prolong = GMG_P_Prolongation::Injection;
					args.Solver.MG.GMG_P_Restrict = GMG_P_Restriction::RemoveHigherOrders;
				}
				else if (code == 3)
				{
					args.Solver.MG.HP_CS = HP_CoarsStgy::P_then_H;
					args.Solver.MG.P_CS = P_CoarsStgy::Minus2;
					args.Solver.MG.GMG_P_Prolong = GMG_P_Prolongation::H_Prolongation;
					args.Solver.MG.GMG_P_Restrict = GMG_P_Restriction::P_Transpose;
				}
				else if (code == 4)
				{
					args.Solver.MG.HP_CS = HP_CoarsStgy::HP_then_H;
					args.Solver.MG.P_CS = P_CoarsStgy::Minus1;
				}
				else if (code == 5)
				{
					args.Solver.MG.HP_CS = HP_CoarsStgy::Alternate;
					args.Solver.MG.P_CS = P_CoarsStgy::Minus2;
				}
				else
					argument_error("unknown hp-configuration code. Check -hp-config argument.");
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
				defaultCoarseOperator = false;
				break;
			case OPT_Smoothers:
			{
				string s(optarg);
				auto pos = s.find(",");
				args.Solver.MG.PreSmootherCode = s.substr(0, s.find(","));
				args.Solver.MG.PostSmootherCode = s.substr(pos + 1);
				break;
			}
			case OPT_HP_CS:
			{
				string code = optarg;
				if (code.compare("h") == 0)
					args.Solver.MG.HP_CS = HP_CoarsStgy::H_only;
				else if (code.compare("p") == 0)
					args.Solver.MG.HP_CS = HP_CoarsStgy::P_only;
				else if (code.compare("p_h") == 0)
					args.Solver.MG.HP_CS = HP_CoarsStgy::P_then_H;
				else if (code.compare("h_p") == 0)
					args.Solver.MG.HP_CS = HP_CoarsStgy::H_then_P;
				else if (code.compare("p_hp") == 0)
					args.Solver.MG.HP_CS = HP_CoarsStgy::P_then_HP;
				else if (code.compare("hp_h") == 0)
					args.Solver.MG.HP_CS = HP_CoarsStgy::HP_then_H;
				else if (code.compare("hp_p") == 0)
					args.Solver.MG.HP_CS = HP_CoarsStgy::HP_then_P;
				else if (code.compare("alt") == 0)
					args.Solver.MG.HP_CS = HP_CoarsStgy::Alternate;
				else
					argument_error("unknown hp-coarsening strategy code '" + code + "'. Check -hp-cs argument.");
				break;
			}
			case OPT_H_CS:
			{
				string coarseningStgyCode = optarg;
				if (coarseningStgyCode.compare("s") == 0)
					args.Solver.MG.H_CS = H_CoarsStgy::StandardCoarsening;
				else if (coarseningStgyCode.compare("r") == 0)
					args.Solver.MG.H_CS = H_CoarsStgy::GMSHSplittingRefinement;
				else if (coarseningStgyCode.compare("b") == 0)
					args.Solver.MG.H_CS = H_CoarsStgy::BeyRefinement;
				else if (coarseningStgyCode.compare("m") == 0)
					args.Solver.MG.H_CS = H_CoarsStgy::IndependentRemeshing;
				else if (coarseningStgyCode.compare("n") == 0)
					args.Solver.MG.H_CS = H_CoarsStgy::AgglomerationCoarseningByFaceNeighbours;
				else if (coarseningStgyCode.compare("mn") == 0)
					args.Solver.MG.H_CS = H_CoarsStgy::MultipleAgglomerationCoarseningByFaceNeighbours;
				else if (coarseningStgyCode.compare("dpa") == 0)
					args.Solver.MG.H_CS = H_CoarsStgy::DoublePairwiseAggregation;
				else if (coarseningStgyCode.compare("mpa") == 0)
					args.Solver.MG.H_CS = H_CoarsStgy::MultiplePairwiseAggregation;
				else if (coarseningStgyCode.compare("vr") == 0)
					args.Solver.MG.H_CS = H_CoarsStgy::AgglomerationCoarseningByVertexRemoval;
				else if (coarseningStgyCode.compare("mcf") == 0)
					args.Solver.MG.H_CS = H_CoarsStgy::AgglomerationCoarseningByMostCoplanarFaces;
				else if (coarseningStgyCode.compare("cc") == 0)
					args.Solver.MG.H_CS = H_CoarsStgy::AgglomerationCoarseningByClosestCenter;
				else if (coarseningStgyCode.compare("clf") == 0)
					args.Solver.MG.H_CS = H_CoarsStgy::AgglomerationCoarseningByClosestFace;
				else if (coarseningStgyCode.compare("li") == 0)
					args.Solver.MG.H_CS = H_CoarsStgy::AgglomerationCoarseningByLargestInterface;
				else if (coarseningStgyCode.compare("seed") == 0)
					args.Solver.MG.H_CS = H_CoarsStgy::AgglomerationCoarseningBySeedPoints;
				else if (coarseningStgyCode.compare("vn") == 0)
					args.Solver.MG.H_CS = H_CoarsStgy::AgglomerationCoarseningByVertexNeighbours;
				else if (coarseningStgyCode.compare("f") == 0)
					args.Solver.MG.H_CS = H_CoarsStgy::FaceCoarsening;
				else
					argument_error("unknown coarsening strategy code '" + coarseningStgyCode + "'. Check -cs argument.");
				break;
			}
			case OPT_P_CS:
			{
				string code = optarg;
				if (code.compare("-1") == 0)
					args.Solver.MG.P_CS = P_CoarsStgy::Minus1;
				else if (code.compare("-2") == 0)
					args.Solver.MG.P_CS = P_CoarsStgy::Minus2;
				else if (code.compare("/2") == 0)
					args.Solver.MG.P_CS = P_CoarsStgy::DivideBy2;
				else if (code.compare("=0") == 0)
					args.Solver.MG.P_CS = P_CoarsStgy::DirectToLow;
				else
					argument_error("unknown p-coarsening strategy code '" + code + "'. Check -p-cs argument.");
				break;
			}
			case OPT_FaceCoarseningStrategy:
			{
				string fcsCode = optarg;
				if (fcsCode.compare("c") == 0)
					args.Solver.MG.FaceCoarseningStgy = FaceCoarseningStrategy::InterfaceCollapsing;
				else if (fcsCode.compare("n") == 0)
					args.Solver.MG.FaceCoarseningStgy = FaceCoarseningStrategy::None;
				else if (fcsCode.compare("i") == 0)
					args.Solver.MG.FaceCoarseningStgy = FaceCoarseningStrategy::InterfaceCollapsingAndTryAggregInteriorToInterfaces;
				else
					argument_error("unknown face coarsening strategy code '" + fcsCode + "'. Check -fcs argument.");
				break;
			}
			case OPT_NumberOfMeshes:
				args.Solver.MG.NumberOfMeshes = atoi(optarg);
				break;
			case OPT_CoarseSolver:
			{
				defaultCoarseSolver = false;
				string code = optarg;
				args.Solver.MG.CoarseSolverCode = code;
				break;
			}
			case OPT_BoundaryFaceCollapsing:
			{
				string code = optarg;
				if (code.compare("d") == 0)
					args.Solver.MG.BoundaryFaceCollapsing = FaceCollapsing::Disabled;
				else if (code.compare("c") == 0)
					args.Solver.MG.BoundaryFaceCollapsing = FaceCollapsing::OnlyCollinear;
				else if (code.compare("p") == 0)
					args.Solver.MG.BoundaryFaceCollapsing = FaceCollapsing::ByPairs;
				else if (code.compare("m") == 0)
					args.Solver.MG.BoundaryFaceCollapsing = FaceCollapsing::Max;
				else
					argument_error("unknown boundary face collapsing code '" + code + "'. Check -bfc argument.");
				break;
			}
			case OPT_ReEntrantCornerManagement:
			{
				string code = optarg;
				if (code.compare("d") == 0)
					args.Solver.MG.ReEntrantCornerManagement = ReEntrantCornerMgmt::Disabled;
				else if (code.compare("f") == 0)
					args.Solver.MG.ReEntrantCornerManagement = ReEntrantCornerMgmt::AgglomerateFirst;
				else
					argument_error("unknown re-entrant corner management code '" + code + "'. Check -rcm argument.");
				break;
			}
			case OPT_CoarseningFactor:
				args.Solver.MG.CoarseningFactor = atof(optarg);
				break;
			case OPT_CoarseN:
				args.Solver.MG.CoarseN = stoul(optarg, nullptr, 0);
				break;
			case OPT_ApproxL2ProjNSubtriangulations:
				args.Solver.MG.NSubtriangulationsForApproxL2Proj = atoi(optarg);
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
					else if (code.compare("errgmsh") == 0)
						args.Actions.ExportErrorToGMSH = true;
					else if (code.compare("sourcegmsh") == 0)
						args.Actions.ExportSourceToGMSH = true;
					else if (code.compare("mg") == 0)
						args.Actions.ExportMultigridComponents = true;
					else if (code.compare("mgit") == 0)
						args.Actions.ExportMultigridIterationVectors = true;
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
	//                 Problem                  //
	//------------------------------------------//

	// dimension
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
			args.Problem.TestCaseCode = args.Problem.GeoCode;
		else
			args.Problem.TestCaseCode = FileSystem::FileNameWithoutExtension(args.Problem.GeoCode);
	}

	// Heterogeneity
	if (args.Problem.GeoCode.compare("square") == 0 && args.Problem.HeterogeneityRatio != 1)
		Utils::Warning("The geometry 'square' has only one physical part: -heterog argument is ignored. Use 'square4quadrants' instead to run an heterogeneous problem.");
	

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
	//              Discretization              //
	//------------------------------------------//

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

	// Elem polynomial bases
	if (args.Discretization.ElemBasisCode.empty())
	{
		if (args.Discretization.OrthogonalizeElemBasesCode == -1)
			args.Discretization.OrthogonalizeElemBasesCode = args.Discretization.MeshCode.compare("cart") == 0 ? 0 : 1;
		args.Discretization.ElemBasisCode = args.Discretization.OrthogonalizeElemBasesCode > 0 ? "monomials" : "legendre";
	}
	else if (args.Discretization.OrthogonalizeElemBasesCode == -1)
	{
		if (args.Discretization.ElemBasisCode.compare("legendre") == 0 && args.Discretization.MeshCode.compare("cart") == 0)
			args.Discretization.OrthogonalizeElemBasesCode = 0;
		else
			args.Discretization.OrthogonalizeElemBasesCode = 1;
	}

	// Face polynomial bases
	if (args.Discretization.FaceBasisCode.empty())
	{
		if (args.Discretization.OrthogonalizeFaceBasesCode == -1)
		{
			if (args.Problem.Dimension <= 2 || args.Discretization.MeshCode.compare("cart") == 0)
				args.Discretization.OrthogonalizeFaceBasesCode = 0;
			else
				args.Discretization.OrthogonalizeFaceBasesCode = 1;
		}
		args.Discretization.FaceBasisCode = args.Discretization.OrthogonalizeFaceBasesCode > 0 ? "monomials" : "legendre";
	}
	else if (args.Discretization.OrthogonalizeFaceBasesCode == -1)
	{
		if (args.Discretization.FaceBasisCode.compare("legendre") == 0 && args.Discretization.MeshCode.compare("cart") == 0)
			args.Discretization.OrthogonalizeFaceBasesCode = 0;
		else
			args.Discretization.OrthogonalizeFaceBasesCode = 1;
	}

	//------------------------------------------//
	//                  Solver                  //
	//------------------------------------------//


	if (args.Solver.SolverCode.compare("default") == 0)
	{
		if (args.Discretization.Method.compare("hho") == 0 && args.Discretization.StaticCondensation && args.Problem.Dimension > 1)
		{
			args.Solver.SolverCode = "mg";
			if ((args.Solver.MG.GMG_H_Prolong == GMG_H_Prolongation::Wildey || args.Solver.MG.GMG_H_Prolong == GMG_H_Prolongation::FaceInject) && !args.Solver.MG.UseGalerkinOperator)
			{
				Utils::Warning("The multigrid with prolongation code " + to_string((unsigned)args.Solver.MG.GMG_H_Prolong) + " requires the Galerkin operator. Option -g 0 ignored.");
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

	if (args.Solver.SolverCode.compare("mg") == 0 || args.Solver.SolverCode.compare("cgmg") == 0 || args.Solver.SolverCode.compare("fcgmg") == 0 || args.Solver.SolverCode.compare("p_mg") == 0)
	{
		if (args.Discretization.Method.compare("dg") == 0)
			argument_error("Multigrid only applicable on HHO discretization.");

		if (!args.Discretization.StaticCondensation)
			argument_error("Multigrid only applicable if the static condensation is enabled.");

		args.Solver.MG.GMG_H_Prolong = static_cast<GMG_H_Prolongation>(args.Solver.MG.ProlongationCode);

		if (args.Solver.MG.GMG_H_Prolong == GMG_H_Prolongation::Wildey && !args.Solver.MG.UseGalerkinOperator)
			argument_error("To use the prolongationCode " + to_string((unsigned)GMG_H_Prolongation::Wildey) + ", you must also use the Galerkin operator. To do so, add option -g 1.");

		if (args.Solver.MG.H_CS == H_CoarsStgy::FaceCoarsening && !args.Solver.MG.UseGalerkinOperator)
			argument_error("To use the face coarsening, you must also use the Galerkin operator. To do so, add option -g 1.");

		if (args.Solver.MG.H_CS == H_CoarsStgy::IndependentRemeshing && 
			Utils::RequiresNestedHierarchy(args.Solver.MG.GMG_H_Prolong) &&
			args.Solver.MG.GMG_H_Prolong != GMG_H_Prolongation::Default)
			argument_error("The coarsening by independent remeshing is only applicable with the non-nested versions of the multigrid (-prolong " + to_string((unsigned)GMG_H_Prolongation::CellInterp_ExactL2proj_Trace) + ", " + to_string((unsigned)GMG_H_Prolongation::CellInterp_ApproxL2proj_Trace) + " or " + to_string((unsigned)GMG_H_Prolongation::CellInterp_FinerApproxL2proj_Trace) + ").");

		if (args.Solver.SolverCode.compare("p_mg") == 0 && defaultCoarseSolver)
			args.Solver.MG.CoarseSolverCode = "mg";

		if (args.Solver.MG.H_CS == H_CoarsStgy::None)
		{
			if (args.Problem.Dimension < 3)
			{
				if (args.Solver.MG.GMG_H_Prolong == GMG_H_Prolongation::Default)
				{
					if (args.Discretization.Mesher.compare("inhouse") == 0)
						args.Solver.MG.H_CS = H_CoarsStgy::StandardCoarsening;
					else
						args.Solver.MG.H_CS = H_CoarsStgy::IndependentRemeshing;
				}
				else if (Utils::RequiresNestedHierarchy(args.Solver.MG.GMG_H_Prolong))
					args.Solver.MG.H_CS = H_CoarsStgy::GMSHSplittingRefinement;
				else
					args.Solver.MG.H_CS = H_CoarsStgy::IndependentRemeshing;
			}
			else
			{
				if (args.Discretization.Mesher.compare("inhouse") == 0 && args.Discretization.MeshCode.compare("tetra") == 0)
					args.Solver.MG.H_CS = H_CoarsStgy::BeyRefinement;
				else if (args.Discretization.Mesher.compare("inhouse") == 0 && args.Discretization.MeshCode.compare("cart") == 0)
					args.Solver.MG.H_CS = H_CoarsStgy::StandardCoarsening;
				else
					args.Solver.MG.H_CS = H_CoarsStgy::IndependentRemeshing;
			}
		}

		if (args.Solver.MG.GMG_H_Prolong == GMG_H_Prolongation::Default)
		{
			if (Utils::BuildsNestedMeshHierarchy(args.Solver.MG.H_CS))
				args.Solver.MG.GMG_H_Prolong = GMG_H_Prolongation::CellInterp_Trace;
			else
				args.Solver.MG.GMG_H_Prolong = GMG_H_Prolongation::CellInterp_FinerApproxL2proj_Trace;
		}

		if (defaultCycle)
		{
			args.Solver.MG.PreSmoothingIterations = 0;
			args.Solver.MG.PostSmoothingIterations = args.Problem.Dimension < 3 ? 3 : 6;
		}
	}

	if (args.Solver.SolverCode.compare("uamg") == 0 || args.Solver.SolverCode.compare("fcguamg") == 0)
	{
		if (args.Discretization.Method.compare("dg") == 0)
			argument_error("Multigrid only applicable on HHO discretization.");

		if (args.Solver.MG.ProlongationCode == 0)
			args.Solver.MG.UAMGMultigridProlong = UAMGProlongation::ChainedCoarseningProlongations;
		else
			args.Solver.MG.UAMGMultigridProlong = static_cast<UAMGProlongation>(args.Solver.MG.ProlongationCode);
		args.Solver.MG.UAMGFaceProlong = args.Solver.MG.FaceProlongationCode == 0 ? UAMGFaceProlongation::BoundaryAggregatesInteriorAverage : static_cast<UAMGFaceProlongation>(args.Solver.MG.FaceProlongationCode);
		args.Solver.MG.UAMGCoarseningProlong = args.Solver.MG.CoarseningProlongationCode == 0 ? UAMGProlongation::ReconstructSmoothedTraceOrInject : static_cast<UAMGProlongation>(args.Solver.MG.CoarseningProlongationCode);

		if (args.Solver.MG.H_CS == H_CoarsStgy::None)
			args.Solver.MG.H_CS = H_CoarsStgy::MultiplePairwiseAggregation;

		if ((args.Solver.MG.H_CS == H_CoarsStgy::MultiplePairwiseAggregation || 
			 args.Solver.MG.H_CS == H_CoarsStgy::MultipleAgglomerationCoarseningByFaceNeighbours)
			&& args.Solver.MG.CoarseningFactor == 0)
			args.Solver.MG.CoarseningFactor = 3.8;

		if (defaultCoarseOperator)
			args.Solver.MG.UseGalerkinOperator = true;

		if (defaultCycle)
			args.Solver.MG.CycleLetter = 'K';
	}

	if (args.Solver.SolverCode.compare("aggregamg") == 0 || args.Solver.SolverCode.compare("fcgaggregamg") == 0)
	{
		if (!defaultCoarseOperator && !args.Solver.MG.UseGalerkinOperator)
			Utils::Warning("AggregAMG uses the Galerkin operator. Argument -g 0 ignored.");

		if (args.Solver.MG.H_CS == H_CoarsStgy::None)
			args.Solver.MG.H_CS = H_CoarsStgy::DoublePairwiseAggregation;

		if (defaultCycle)
			args.Solver.MG.CycleLetter = 'K';
	}


	if (args.Solver.MG.H_CS == H_CoarsStgy::MultiplePairwiseAggregation && args.Solver.MG.CoarseningFactor == 0)
		args.Solver.MG.CoarseningFactor = 3.5;
	else if (args.Solver.MG.H_CS == H_CoarsStgy::IndependentRemeshing && args.Solver.MG.CoarseningFactor == 0)
		args.Solver.MG.CoarseningFactor = 2;
	
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
