# The fhhos4 multigrid solver as an external library

1. [The algorithm in a nutshell](#1-the-algorithm-in-a-nutshell)
2. [Applicability](#2-applicability)
3. [How to use](#3-how-to-use)
4. [License](#license)

# 1. The algorithm in a nutshell

The shared library `fhhos4` offers an algebraic multigrid method for the solution of linear systems 
arising from hybrid discretizations (HHO, HDG...) of diffusion and convection-diffusion problems, scalar or
vector-valued, of the form
```math
\begin{pmatrix} A_{TT} & A_{TF} \\ A_{FT} & A_{FF} \end{pmatrix}
\begin{pmatrix} x_T \\ x_F \end{pmatrix}
=
\begin{pmatrix} b_T \\ b_F \end{pmatrix}
```
While the smoothers work on the condensed system $A x_F = b$ (given by the elimination of the cell unknowns)
```math
A = A_{FF} - A_{FT} A_{TT}^{-1} A_{TF}, \qquad b = b_F - A_{FT} A_{TT}^{-1} b_T
```
the aggregation method used for the construction of coarse levels requires the blocks of the uncondensed matrix
(to recover the cell-to-cell connectivity amongst other things). Consequently, instead of giving only the condensed system like any black-box AMG method, 
the blocks of the uncondensed matrix must in our case be passed to the solver.

The multigrid hierarchy unfolds as follows:
First, the high polynomial orders are handled by p-multigrid levels, successively reducing the order to the lowest one.
Then, once the system is reduced to one unknown per cell and face, h-coarsening follows, based on the algorithm described in

> D. A. Di Pietro, F. Hülsemann, P. Matalon, P. Mycek, U. Rüde, *Algebraic multigrid preconditioner for statically
> condensed systems arising from lowest-order hybrid discretizations*, SIAM J. Sci. Comput., 2023.

# 2. Applicability

- [x] Linear systems:
   - [x] Symmetric positive definite (e.g. diffusion)
   - [x] Non-symmetric, indefinite systems (e.g. convection-diffusion, both diffusion- and convection-dominated)
   - [x] Singular, symmetric systems
      - [x] with translation kernel (e.g. pure Neumann diffusion)
      - [ ] user-defined kernel vectors (e.g. for elasticity problems)
- [x] Vector-valued unknowns
   - [x] Linear elasticity, with the rigid body modes given as the near null space
- [x] `P` and `Q` polynomial spaces in the cells and faces
- [x] Arbitrary polynomial degree:
   - [x] Hierarchical bases **ordered by degree**
   - [x] If not ordered by degree: the degree 0 must at least be first. Then OK up to degree 2.
   - [x] Helper to reorder the bases
   - [ ] Non-hierarchical bases not managed
   - [ ] p-adapdability not managed (fixed polynomial degrees across the mesh!)

# 3. How to use

Your code assembles the uncondensed system, i.e.
- the matrix blocks $A_{TT}$, $A_{TF}$, $A_{FF}$ for a symmetric system, and $A_{FT}$ for a non-symmetric one,
- the two blocks of the right-hand side $(b_T, b_F)$, without the DoFs of the Dirichlet faces (moved to the right-hand side).

That's it, no static condensation.
The library condenses the system internally and solves it for both the
cell and face unknowns $(x_T, x_F)$.

Then, by the interface of your code:

| Interface | Page |
|---|---|
| C++: your code links the library | [The library in C++](library/doc/cpp.md): installation, link, symmetric and non-symmetric systems, use as a preconditioner |
| Python: NumPy arrays and SciPy matrices | [Python](library/doc/python.md) |
| PETSc: `-pc_type fhhos4`, a preconditioner of its Krylov methods | [PETSc](library/doc/petsc.md) |

Common to the three:
- [The system to give](library/doc/system.md): the blocks and their format, the near null space, the bases, what the
  setup infers, and the helpers for other matrix formats.
- [Several processes (MPI)](library/doc/mpi.md): who owns what, the ghost faces, what changes in the solver.
- [Parameters](library/doc/parameters.md) and constraints.

# License

LGPL-3.0-or-later, like the rest of fhhos4 (see [README.md](README.md)). The library only depends on Eigen (MPL-2.0),
and, compiled with MPI or `ENABLE_MUMPS`, on MPI and MUMPS (CeCILL-C, compatible with the LGPL; MUMPS's own
dependencies depend on its build: conda-forge's `mumps-mpi` brings ParMETIS, whose license forbids commercial use), not
on the GPL dependencies of the program (GMSH, CGAL): codes under other licenses can link it, under the conditions
of section 4 of the LGPL (in short: keep the license notices, and let the user replace the shared library).
