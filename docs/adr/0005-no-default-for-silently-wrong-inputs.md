# No default for an input whose wrong value silently breaks the solver

The interpolations of 1 of the cells and faces, `Dimension`, `FaceDegree`, `CellDegree` and, with MPI, `Comm` have no
default: the caller must give them (2026-10-14, 2026-10-19, 2026-10-20). A wrong interpolation of 1 does not fail, it
slows U-AMG down or makes it diverge (L2-orthonormal bases with c = 1: Cube-tet n=16, k=0, 26 -> 44 iterations; graded
meshes fail), and it cannot be recovered from the blocks (ambiguity c*h, h discrete harmonic). "It's too risky to
leave the user the possibility to leave the solver with very wrong values." A compact form (one constant for all) only
if a measured gain justifies it.

Exception: `SingularMatrix` defaults to `false` (the user, 2026-10-27), although a singular matrix left with `false`
converges erratically or not at all. The library cannot know that a matrix is singular: a test on A z needs a threshold
and cannot tell a tiny reaction term from none. Most clients never solve a singular system, and the wrong value the
other way is caught: with `true`, the setup checks that the interpolation of 1 is in the kernel. Not `KernelDimension`
with vectors: the kernel is the near-kernel the solver already has (the interpolation of 1; PETSc's `has_cnst`), and a
kernel given apart would have to equal it; for elasticity, the rigid body modes, and a kernel that is a part of them
(sliding conditions) would need vectors then (ADR-0020).
