# The discretization inferred from the stored pattern of the blocks

`Dimension`, `FaceDegree`, `CellDegree` and `PolynomialSpace`, when they are not set, are inferred by the setups from the
blocks (2026-11-11 to 2026-11-14, the user: "a user trying fhhos4 should not have to set every parameter of the
discretization first"; the use case is `-pc_type fhhos4` with almost no option, but every interface has it: C++,
Python, PETSc). The documentation still tells the users to set them.

- **The DoFs per cell and per face come from the stored entries of `A_TT` and `A_TF`**, not from their values (exact
  zeros inside a block are common: Legendre bases on a square cell), nor from the interpolations of 1 ("unreliable";
  PETSc gives the interpolation of 1 on all the DoFs, which needs the block sizes first). The zeros inside the blocks
  of each cell must be stored, as the clients store them (HArDCore3D, Disk++, PETSc by default, NGSolve).
- **The readings are the discretizations the library takes**: 2D or 3D, P or Q, `CellDegree` k - 1, k or k + 1. Hence
  the library refuses any other `CellDegree`, on every path: a discretization outside the readings could otherwise be
  read silently as one inside.
- **Any ambiguity is an error**: readings that give different solvers (the block sizes of the p-levels of degree >= 1,
  or k = 0 against k >= 1), whose message lists them and names the parameter that resolves them. Readings that give the
  same solver are not one (k = 0: the dimension doesn't matter): the parameters they don't determine take the value of
  any reading. Up to k = 10, one parameter always resolves an ambiguity.
- **`VectorComponents` is never inferred**: not from the sizes (39 of 112 common configurations ambiguous), nor from the
  number of near-kernel vectors (elasticity: 3 in 2D, 6 in 3D). It is required when `Dimension` or a degree is inferred.
  With the three set, it is 1 if not set, and the space alone comes from the stored entries too: Q if they show its
  DoFs, P otherwise (since 2026-11-14: the near null space is given on all the DoFs, and no longer gives the number of
  cells and faces; before, from the rows of `A_TT` per coefficient of the interpolation of 1, without stored zeros).
  A wrong space is then refused by the checks of the sizes or of the near null space.

An inferred value is not a default (ADR-0005): it is determined by the blocks, or refused.

## Considered Options

- `PolynomialSpace::P` kept as the default: a Q discretization with nothing set could be read silently as a P one (2D Q
  k=2 l=1 has the DoFs of 3D P k=1 l=1, whose solvers differ with `PCoarsening::Minus1`). Hence
  `PolynomialSpace::Inferred`, the default.
- A public function to infer the discretization (the user, 2026-11-14: not for the user to call): the setups infer; the
  PETSc preconditioner, which needs the block sizes before it merges the blocks, calls the internal functions itself.
- The analysis when every parameter is set: none, there is nothing to infer.

## Consequences

- The setup writes the values found into the members of `Solver` (printed at Verbosity >= 1): a later setup of the same
  solver reads them, with no inference (PETSc's `PCSetUp()` on new values of the same operator); on another
  discretization, the checks of the sizes refuse it.
- On several processes, each process analyzes its own cells and faces, then one gather; a process with no cell has no
  opinion. The DoFs per face need a face between two cells of a process: without one (one cell per process), an error.
- The PETSc preconditioner, as the setups since 2026-11-14, has no number of cells and faces (the near null space is
  on all the DoFs): with the dimension and the degrees given, it reads the space from the stored
  entries too, unless P and Q have the same DoFs and solver (e.g. k = 0), where the space doesn't matter; P if they
  don't give it (as the setups, which never refuse the dimension and the degrees without the space).
- The PETSc preconditioner includes an internal header (`src/Solver/DiscretizationInference.h`, header-only, compiled
  into it): moving it into PETSc (ADR-0012) means installing that header with the library's.
