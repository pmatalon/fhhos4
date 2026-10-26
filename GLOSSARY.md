# fhhos4

Solvers for hybrid discretizations of diffusion problems, centred on algebraic multigrid: a research program that
reproduces the papers, and a library that gives the multigrid to other codes.

## Discretization

**HHO** (Hybrid High-Order):
A discretization whose unknowns are polynomials on the cells and on the faces of the mesh.

**Cell unknowns**, **face unknowns**:
The coefficients of the polynomial of a cell, of a face, in its basis.
_Avoid_: element DoFs, facet DoFs, skeleton unknowns

**Face degree** (k), **cell degree** (l):
The polynomial degrees of the face and cell unknowns.
_Avoid_: order

**Polynomial space**:
The space the bases span: P (total degree) or Q (largest partial degree).

**Hierarchical basis**:
A basis whose first functions span the polynomials of lower degree, the constant first.

**Interpolation of 1**:
The coefficient c of each cell and each face such that the constant function 1 = c phi_0, phi_0 being the first basis
function.
_Avoid_: scaling, constant value

**Blocks**:
The four parts of the uncondensed system: `A_TT` (cell-cell, block-diagonal), `A_TF` (cell-face), `A_FT` (face-cell),
`A_FF` (face-face).

**Uncondensed system**:
The linear system on the cell and face unknowns, given by its blocks.
_Avoid_: interface 1, full system

**Static condensation**:
The elimination of the cell unknowns, cell by cell, which leaves a system on the face unknowns.

**Condensed system**:
The system on the face unknowns that static condensation gives (the Schur complement).
_Avoid_: interface 2, global matrix, Schur system

**Boundary face**:
A face of a single cell that stays in the system (Neumann, Robin conditions). Dirichlet faces are not in the system.

**Dirichlet lifting**:
Moving the known Dirichlet values to the right-hand side; always the caller's job.

**Full Neumann problem**:
A pure Neumann problem: its system is singular, its solution defined up to a constant.

## Multigrid

**U-AMG**:
The uncondensed algebraic multigrid: an aggregation multigrid on the condensed system, built from the blocks of the
uncondensed system, whose coarse levels are again hybrid (cells and faces).
_Avoid_: UAMG, fhhos4_AMG

**Level**:
One operator of the multigrid hierarchy, with its smoother and its transfers to the next.

**p-level**:
A level on the same mesh with a lower face degree.

**h-level**:
A level on a coarser mesh, whose cells are aggregates.

**Aggregate**:
A set of cells of a level that becomes one cell of the next h-level.

**Coarsening pass**:
One pairwise aggregation of cells; an h-level results from several passes.
_Avoid_: round, sweep

**Coupling**:
The strength of the connection of two cells through a face, measured on the blocks, which drives the aggregation.

**Coarse face**:
A face of an h-level: the fine faces between two aggregates, collapsed into one (interface collapsing).

**Removed face**:
A fine face that is not part of any coarse face: interior to an aggregate, or a boundary face.

**Coarsening prolongation**:
The map from the coarse faces of a pass to its fine faces.

**Multigrid prolongation**:
The map from a level to the next finer one: the product of the coarsening prolongations of its passes.

**K-cycle**:
A multigrid cycle with Krylov iterations on the coarse levels, which makes the multigrid a non-linear preconditioner.

**Smoother**:
The iterative method applied on each level before and after the coarse correction.

**Downwind order**:
An order of the faces that follows the flow, used by the Gauss-Seidel sweeps of a non-symmetric system.

**Mixed precision**:
The multigrid levels stored and applied in single precision, the outer Krylov method in double.

## Parallelism

**Hybrid smoother**:
Gauss-Seidel within the rows of each thread, Jacobi between the threads.

**Owner**:
The process that holds a face's row in a distributed system: the smaller rank among its cells' processes.

**Ghost face**:
A face that a process uses without owning it.

## Project

**Program**:
The executable `fhhos4`: the research test bed, which reproduces the papers.

**Library**:
`libfhhos4`: U-AMG for other codes.
_Avoid_: fhhos4_AMG

**Client code**:
A code that calls the library: HArDCore3D, Disk++, NGSolve (through the adapter).
_Avoid_: binding

**Adapter**:
The NGSolve operator `UAMG(a)` that wraps the library, outside the repository.

**Papers' iteration counts**:
The iteration counts of `reproducibility/*.md`, which the tests check and the default options of the program must keep.
