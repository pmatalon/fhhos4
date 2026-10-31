# The Python module

Read before changing `library/python/`. Decisions: `docs/adr/0014-python-module-mirrors-the-cpp-api.md`, 0012. Its
first client is NGSolve's HDG (`L2 * FacetFESpace`, primal, SPD), through an adapter outside the repository (ADR-0012).

## What a client framework needs from the module

What NGSolve showed, likely true of other frameworks:
- Its arrays read in place as they are: uint64 row pointers and int32 columns (`a.mat.CSR()`), hence the index types
  of `PyCSRView`.
- Bases not ordered by degree (e.g. Dubiner's order): the caller sorts the DoFs of each cell and face by degree in the
  DoF lists of `ExtractSymmetricBlocks()`, otherwise the p-levels silently use the wrong low-order space
  (`PCoarsening::ToZero` is the fallback, slower from k = 3). phi_0 may differ from 1: the interpolations of 1 are
  1/phi_0.
- On several processes: a process may own nothing (NGSolve's rank 0 holds no element), and the framework's matrices
  are sub-assembled (a shared face's row holds the local contributions only): the caller merges them onto the owners.
  Time the collective calls after a barrier, otherwise a process waits inside the setup for the others' preparation.

## Code

`library/python/fhhos4_python.cpp`: the module `fhhos4` (nanobind, `NB_STATIC`), the public API with the C++ names, on
NumPy arrays read in place (`nb::try_cast` without conversion; the vectors `noconvert()`, float64 only):
- the matrices as `CSRView` (`PyCSRView`: the arrays of any index types, int32/int64, uint64 read as int64, and their
  Python owners), `CSRMatrix` (arrays exposed as read-only NumPy views that keep it alive), or any object with
  `indptr`, `indices`, `data`, `shape` and `format` "csr" (SciPy, duck-typed: SciPy is needed only by the test);
- `PySolver` derives from `fhhos4::Solver` and keeps the views the solves read (`A`, `A_FT`), replaced only after a
  successful setup: the checks of a setup run before it resets the previous one, which still reads them;
- `CondensedMatrix()`: a copy (`ToCSR()`); the setups, solves and helpers release the GIL.

Built with `-DENABLE_PYTHON=ON` (in the library section of `CMakeLists.txt`, also with `ENABLE_PROGRAM=OFF`) into
`<build>/python`, not installed. Test: `ctest -R Python` (`PythonModule`: `library/python/tests/test_module.py`, the
toy system of `UAMGLibraryTest.cpp` against SciPy's direct solver).

## With MPI (`FHHOS4_MPI`)

`Solver.Comm` from an mpi4py communicator (`py2f()`, `MPI_Comm_f2c`; refused if MPI is not initialized; the getter
returns the Python object), not set: one process, as in C++; `SetupSym(..., faceGhosts=None)`,
`SetupNonSym(..., faceGhosts=None, cellGhosts=None)` (int64 arrays read in place, `ToGhostNumbers()`);
`fhhos4.WithMPI`. Tests: `PythonModule` without `Comm` nor mpi4py (MPI not initialized: `CoarseSolver.MUMPS` refused);
`PythonModuleMPI_np2`/`_np4` (`test_module_mpi.py`, under `mpiexec`: the toy system distributed by strips, also with an
empty rank 0, symmetric or not, against SciPy; an input error on one process raised on all).

## Environments and build directories

- Conda env `fhhos4-python` (python 3.12, numpy, scipy, nanobind, cmake, ninja, cxx-compiler, eigen 5.0.1; mpich 5.0.1,
  mpi4py 4.1.2 and MUMPS 5.8.2, installed as in `docs/mpi.md`; NGSolve from pip, not on conda-forge), not the `fhhos4`
  env.
- `build-python` (`-DENABLE_PROGRAM=OFF -DENABLE_PYTHON=ON`, ~1 min); `build-python-mpi` (the same with
  `-DENABLE_MPI=ON -DENABLE_TESTS=OFF`: no GTest in that env; the GTest MPI tests run in `build-mpi-lib`, env `fhhos4`).
- For NGSolve on several processes: env `fhhos4-python-openmpi` (the same with openmpi 5.0.11 instead of mpich, and
  NGSolve 6.2.2608 from pip: NGSolve's runs on several processes use OpenMPI; MUMPS 5.7.3: `mumps-mpi=5.7.3
  openmpi=5.0.11 _openmp_mutex=*=*_gnu libblas=*=*blis liblapack=*=*netlib`, since its OpenMPI builds 5.8.x require
  OpenBLAS's OpenMP build, i.e. LLVM's OpenMP), build directory
  `build-python-openmpi` (the options of `build-python-mpi`; `ctest -R Python` passes with
  `OMPI_MCA_rmaps_base_oversubscribe=1`).
