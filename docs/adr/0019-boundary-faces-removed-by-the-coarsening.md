# Boundary faces are removed by the coarsening

Faces of a single cell (Neumann, Robin) are removed by the coarsening like the faces interior to an aggregate, their
coarse values given by the prolongation: the default of the program and the library (2026-10-25). Keeping them,
collapsed into one coarse face per aggregate (program option `-uamg-bfc`), costs 10-17% more iterations on
HArDCore3D's tetrahedral meshes and gains at most 1-2 iterations elsewhere (Disk++). Measurements: `PERFORMANCE.md`.
