# Packaging and release

Read before a release or a change of the conda recipe, the release workflow or the CMake package. Version rules:
`CLAUDE.md` "Version". Open work: `.scratch/release/`.

## Release

Change the version in `project()`, commit, then push an annotated tag `vX.Y.Z` (the tags before 1.3.0 are
`Release_X.Y.Z`). `.github/workflows/release.yml` checks the tag against the version, builds, runs `ctest`, the
installation and the library alone with 64-bit indices (`SMALL_INDEX=OFF`) and the conda test client, and in a parallel
job the conda package with its tests, then creates the GitHub release "Release X.Y.Z" (notes: the message of the tag,
which must not mention HArDCore3D). It pins GMSH, Eigen, CGAL, the compiler and rattler-build to the versions the tests
were validated with: update them with the local environment (envs `fhhos4` and `rattler`). After a change of
`release.yml` or of its pins, run it by hand on master (Actions > Release > Run workflow: no release) before pushing
the tag; otherwise push the tag directly. After the release, conda-forge's bot updates the feedstock; the Spack package
needs a PR to spack-packages (`spack checksum --add-to-package fhhos4 X.Y.Z`, then `spack install --test=root` and
`spack test run`).

## Conda package

The library alone (`conda/recipe/recipe.yaml`, rattler-build v1 recipe; `rattler-build build -r conda/recipe -c
conda-forge --experimental`, ~1 min; rattler-build is in the conda env `rattler`, not `fhhos4`).
- `ENABLE_PROGRAM=OFF`, LGPL (`license_file`: `COPYING.LESSER` and `COPYING`, the LGPL needs both). The program is not
  packaged: the name `fhhos4` is the library's, the clients install it.
- Run dependencies: Eigen >= 3.4 (ADR-0013) and the compiler's runtime (`libstdcxx`, `libgcc`; the overlinking warning
  on `libgomp` is normal on conda-forge). `run_exports` pinned to major.minor, the SOVERSION.
- Linux only: not tested elsewhere (macOS would need llvm-openmp, Windows MSVC untested).
- Its tests: `conda/recipe/test`, a CMake client (version, a solve on a lowest-order hybrid system), compiled as is,
  with `-march=native` (a client with other alignment settings than conda-forge's baseline) and with Eigen 3.4.
- The source is the clone (`path: ../..`, without the gitignored files). The version is read from `project()` with
  `load_from_file` (rattler-build `--experimental`). `variants.yaml`: local builds only.

## conda-forge

staged-recipes PR #35238 (2026-10-17), then the feedstock `conda-forge/fhhos4-feedstock` holds its own recipe: this one
with the version written out (no experimental feature there) and the archive of the release tag as source (url and
sha256; its bot updates them at each release). Both take the test client from the sources (`files: source:
conda/recipe/test/`): its changes reach conda-forge with the release; any other change of `recipe.yaml` goes to the
feedstock too.

## Spack

The library alone too, in Spack's builtin repository (`spack/spack-packages`,
`repos/spack_repo/builtin/packages/fhhos4/package.py`; PR 7139), not in this repository. `CMakePackage`,
`ENABLE_PROGRAM=OFF`, `ENABLE_NATIVE_ARCH=OFF` (Spack sets the target's flags), variant `openmp` (default on),
`eigen@3.4:`. No variant for the index type: `SMALL_INDEX=OFF` does not compile in 1.4.0 (fixed after it). Its
stand-alone test (`spack test run fhhos4`) builds and runs the conda test client `conda/recipe/test`, cached from the
sources: moving it changes the package too (per version).
