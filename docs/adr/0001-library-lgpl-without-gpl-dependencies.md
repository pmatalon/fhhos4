# The library is LGPL and depends on Eigen only

fhhos4 is LGPL-3.0-or-later (2026-10-09; rights cleared with CERFACS, holder for 2018-2022, and Paul Mycek, author of
the Hemker basis still in the code) so that codes under other licenses can link the library: HArDCore3D (LGPL),
Disk++ (MPL-2.0), NGSolve, DOLFINx. The program's GPL dependencies (GMSH, CGAL) must never reach the library: its
path includes no mesh or discretization header, and it builds with `ENABLE_PROGRAM=OFF` from Eigen alone.
