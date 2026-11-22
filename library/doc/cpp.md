[← The fhhos4 library](../../README_library.md)

# The library in C++

- [Installation](#installation)
  - [Install from source](#install-from-source)
  - [Conda and Spack packages](#conda-and-spack-packages)
  - [Link](#link)
- [How to use](#how-to-use)
  - [Symmetric systems](#symmetric-systems)
  - [Non-symmetric systems](#non-symmetric-systems)
  - [Use as a preconditioner](#use-as-a-preconditioner)
- Common to the three interfaces: [The system to give](system.md), [Several processes (MPI)](mpi.md),
  [Parameters](parameters.md)

## Installation

### Install from source

For use in shared memory only, you'll need
- a C++17 compiler,
- CMake (>= 3.20),
- Ninja or Make,
- [Eigen](https://eigen.tuxfamily.org) 3.4 or 5.x.

For distributed memory, you'll also need
- MPI
- [MUMPS](https://mumps-solver.org)  (e.g. conda-forge's `mumps-mpi`)

```bash
git clone https://github.com/pmatalon/fhhos4.git
cd fhhos4
git checkout $(git describe --tags --abbrev=0)    # the latest release
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
cmake --install build --prefix <prefix>    # lib/, include/ and the CMake package
```

Between releases, these pages describe the development version (`master`); those of the release are in its clone, or
on GitHub from [the latest release](https://github.com/pmatalon/fhhos4/releases/latest).

To update to a later release: `git fetch && git checkout $(git describe --tags --abbrev=0 origin/master)` in the clone,
then build and install again.

Without installing, the library is `build/lib/libfhhos4.so`.

Build options (`cmake -D<OPTION>=<value>`):

| Option | Default | Effect |
|---|---|---|
| `ENABLE_OPENMP` | `ON` | Parallel (OpenMP) |
| `ENABLE_MPI` | `OFF` | The solver on several processes, see [Several processes](mpi.md). Needs MPI and the MPI version of [MUMPS](https://mumps-solver.org); your code is then compiled and linked with MPI too, and the Python module needs mpi4py. In 1.5.0 and before, also `-DENABLE_TESTS=OFF`, unless GoogleTest is installed |
| `ENABLE_PYTHON` | `OFF` | The Python module, see [Python](python.md) |
| `ENABLE_PETSC_PC` | `OFF` | The preconditioner for PETSc, `libfhhos4_petsc`, see [With PETSc](petsc.md). Needs PETSc (found by `pkg-config`) |
| `ENABLE_MUMPS` | `OFF` | The coarse solver `CoarseSolver::MUMPS` on one process; always on with `ENABLE_MPI`. Without MPI, needs MUMPS's sequential version, e.g. conda-forge's `mumps-seq` |
| `MUMPS_INCLUDE_DIR`, `MUMPS_LIBRARY` | found | MUMPS is found from its header `dmumps_c.h` and its library `libdmumps`, or given by `-DMUMPS_INCLUDE_DIR=<dir> -DMUMPS_LIBRARY=<path>/libdmumps.so` |
| `ENABLE_NATIVE_ARCH` | `OFF` | `-march=native`: optimized for the CPU of the build machine, may not run on others |
| `ENABLE_64BIT_INDEX` | `OFF` | 64-bit indices (`int64_t`) in `fhhos4::SparseMatrix` (`int` otherwise) |

### Conda and Spack packages

Coming soon. Currently waiting for review by the Conda and Spack maintainers.

### Link

```cmake
find_package(fhhos4 REQUIRED)
target_link_libraries(my_code PRIVATE fhhos4::fhhos4)
```

and `-DCMAKE_PREFIX_PATH=<prefix>` (installation), or `-Dfhhos4_DIR=<fhhos4>/build` (build tree). CMake puts the
directory of the library in the RPATH of the executable. A library compiled with MPI brings MPI along
(`MPI::MPI_CXX`). The preconditioner for PETSc, if built: `fhhos4::petsc`.

## How to use

### Symmetric systems

```cpp
#include <fhhos4>

fhhos4::Solver solver;
solver.VectorComponents = 1; // components of the unknown: 1 if scalar, e.g. 3 for a vector in 3D (see Inputs)
// Optional, inferred from the blocks if not set (see Inputs):
solver.Dimension = 3;  // space dimension
solver.FaceDegree = k; // polynomial degree in the faces
solver.CellDegree = k; // polynomial degree in the cells
solver.PolynomialSpace = fhhos4::PolynomialSpace::P; // or Q: bases of Q^k instead of P^k

#if FHHOS4_MPI // 1 if the library is compiled with MPI (fhhos4_version.h)
solver.Comm = MPI_COMM_WORLD;
// The blocks: the rows of the cells and faces of this process, the columns of the faces numbered locally (see With MPI)
std::vector<std::int64_t> faceGhostGlobalNumbers = ...; // or a std::span, an Eigen vector, GhostGlobalNumbers(n, p)
solver.SetupSym(A_TT, A_TF, A_FF, cellNearNullSpace, faceNearNullSpace, faceGhostGlobalNumbers);
#else
solver.SetupSym(A_TT, A_TF, A_FF, cellNearNullSpace, faceNearNullSpace);
#endif

Eigen::VectorXd x_T(b_T.size()), x_F = Eigen::VectorXd::Zero(b_F.size()); // with MPI: of the process
fhhos4::Result result = solver.Solve(b_T, b_F, x_T, x_F); // flexible conjugate gradient preconditioned by the AMG
```

> [!TIP]
> `VectorComponents` is the one parameter to set: in the vast majority of common cases, the setup infers the
> dimension, the degrees and the polynomial space from the blocks (see [Inputs](system.md#inputs) for its assumptions, and the
> rare ambiguities, where an error names the parameter to set).

A complete example, an HHO discretization assembled into the blocks, on one process or several with MPI:
[`examples/fhhos4_solver.cpp`](../../examples/fhhos4_solver.cpp).

The blocks are of type `fhhos4::SparseMatrix`, which is an Eigen sparse matrices (row-major, compressed).
It can also an Eigen map of your own arrays (see [Inputs](system.md#inputs)).

### Non-symmetric systems

What changes for a non-symmetric matrix (e.g. convection-diffusion):

- The block $A_{FT}$ too, given to `SetupNonSym()`:
  ```cpp
  solver.SetupNonSym(A_TT, A_TF, A_FT, A_FF, cellNearNullSpace, faceNearNullSpace);
  ```
- The same algorithm runs, with GCR instead of FCG, and Gauss-Seidel smoothers that follow the direction of the flow
  (downwind order). The default smoother, `L1HybridBlockGaussSeidel`, converges whatever the numbers of threads and
  processes. For convection-dominated problems, `HybridBlockGaussSeidel` usually takes fewer iterations; keep the
  default if the iterations stall, e.g. in 3D at degree 2 or more, or on several processes:
  ```cpp
  solver.PreSmoother = fhhos4::Smoother::HybridBlockGaussSeidel;
  ```
- With MPI: see [Several processes](mpi.md#non-symmetric-systems).

### Use as a preconditioner

The library includes everything to use it as a solver, Krylov method included. 
However, it can also precondition an external Krylov
method working on the condensed system:

```cpp
solver.Krylov = fhhos4::Krylov::None;                 // disable Krylov method in fhhos4
solver.MaxIterations = 1;                             // one cycle per application
solver.SetupSym(A_TT, A_TF, A_FF, cellNearNullSpace, faceNearNullSpace);
Eigen::VectorXd b(b_F.size());
solver.CondenseRHS(b_T, b_F, b);                      // b: the condensed right-hand side

// Your Krylov method on A x_F = b
Eigen::VectorXd x_F = Eigen::VectorXd::Zero(b.size()), r = b, z(b.size()), p, Ap(b.size());
while (r.norm() > tolerance * b.norm())
{
    z.setZero();
    solver.Solve(r, z); // the preconditioner, z = M^{-1} r: one cycle
    p = ...;            // the new search direction, from z
    solver.ApplyCondensedMatrix(p, Ap); // the matrix-vector product, Ap = A p
    ...                 // the updates of x_F and r
}

Eigen::VectorXd x_T(b_T.size());
solver.SolveCellUnknowns(b_T, x_F, x_T); // use fhhos4 to solve the cell unknowns at the end
```

On several processes, the same calls with the vectors of the faces of each process (the dot products of your Krylov
method summed over the processes). On one process, `CondensedMatrix()` gives the matrix itself (an Eigen map).

The default K-cycle is non-linear, which means that your Krylov method **must be flexible** (FCG, FGMRES, GCR).
V- and W-cycles work with any Krylov method (CG, GMRES), but are not recommended for performance reasons.
