# fhhos4

fhhos4 (pronounce 'phosphore') stands for Fast HHO Solver for Research. It is an open-source code that solves diffusion (and
biharmonic) problems with the Hybrid High-Order (HHO) discretizations, with a focus on fast multigrid solvers.

It comes in two parts:

- **The program `fhhos4`**: a full solver for research (HHO and DG discretizations, high-order, unstructured meshes,
  GMSH binding, geometric and algebraic multigrid, many parameters, exports of the linear systems, shared-memory
  parallelism). The numerical experiments of related papers are reproduced following the instructions in the
  `reproducibility` folder.
  Install and use: [README_program.md](README_program.md).
- **The library `fhhos4`**: the algebraic multigrid for hybrid discretizations (U-AMG) for your own HHO/HDG code.
  Only needs Eigen. Install, interfaces and usage: [README_library.md](README_library.md).

# License

Copyright (C) 2018-2022 CERFACS\
Copyright (C) 2023-2026 Pierre Matalon

fhhos4 is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
version (`LGPL-3.0-or-later`): see [COPYING.LESSER](COPYING.LESSER), which supplements the GNU General Public License
version 3 in [COPYING](COPYING).

- The library `fhhos4` only depends on Eigen (MPL-2.0): codes under other licenses can link it (conditions in
  section 4 of the LGPL). Modified versions of the library must be distributed under the LGPL.
- The program `fhhos4` links GMSH (GPL-2.0-or-later) and CGAL (GPL-3.0-or-later for the packages it uses): built with
  them (the default), it is distributed under the GNU GPL version 3 or later.
- AGMG (`-DENABLE_AGMG=ON`) has its own license, not compatible with the GPL: do not distribute a program built with it.
- The quadrature rules of Dunavant and Keast (`src/QuadratureRules/`) are code by John Burkardt, distributed under the
  GNU LGPL.
