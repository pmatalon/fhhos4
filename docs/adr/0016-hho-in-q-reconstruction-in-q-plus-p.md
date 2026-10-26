# HHO in Q: the reconstruction in Q^l + P^{k+1}

HHO with Q^k unknowns (`-poly-space q`, 2026-10-23) reconstructs in Q^l + P^{k+1}, not Q^{k+1}: Q^{k+1} is
inconsistent with the Q^k faces (normal derivatives outside Q^k(F)), and the method did not converge for k >= 1. The
bases of Q are ordered by largest partial degree, hierarchical in Q, which the element's rule "cell basis = first
functions of the reconstruction basis" and the p-levels of U-AMG need. Scope: diffusion, Cartesian meshes, `-kc 0|-1`
(no consistent reconstruction for Q^{k+1} cells), no geometric p-multigrid. The library takes the space as one
parameter for cells and faces (`PolynomialSpace`).
