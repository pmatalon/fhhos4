# One algorithm for symmetric and non-symmetric systems

U-AMG runs the same code on symmetric and non-symmetric matrices; the symmetric case differs only by optimizations that
keep the method (Cholesky instead of LU, lower triangles, lower-triangle local Galerkin products) (2026-10-11: "same
code for both, except the symmetric optimizations"). Only a method that needs symmetry (FCG, the Cholesky coarse
solver) is refused on a non-symmetric system. To extend a solver to a new class of matrices, generalize the existing
algorithm and branch only where symmetry is exploited; measure and ask before making a variant case-specific.

## Considered Options

- The hierarchy built on the symmetric part of the matrix (first version, 2026-10-10): rejected, no symmetric part
  assembled.
- Couplings from the rows of the matrix, as AGMG's strong couplings: 20% more iterations on diffusion-dominated
  convection-diffusion at k = 0, worse without the downwind order. Chosen: the averaged traces of `A_TF` and `A_FT`.
