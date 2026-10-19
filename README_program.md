# fhhos4: the program

Solver for diffusion problems with HHO discretizations, and its algebraic and geometric multigrid methods. Looking for
the multigrid as a library for your own code? See [README_library.md](README_library.md).

## Install

This program depends on the following third-party software packages:
- Eigen
- CGAL (version 5.0 or later)
- GMSH (version 4.9.5 or later)

and optionally
- AGMG (version 3.3.5 or later).

Eigen and CGAL are provided through a conda environment (`conda/environment.yml`). GMSH is also available in that environment, but can alternatively be built from source if you need a specific version.

### 1. Create and activate the conda environment

```bash
conda env create -f conda/environment.yml
conda activate fhhos4
```

Keep this environment activated whenever you configure or build the project, so that CMake can find Eigen and CGAL automatically.

### 2. (Optional) Build GMSH from source

The `fhhos4` conda environment already provides GMSH. Skip this step unless you need a different version.

(See http://gmsh.info/#Download)

```bash
mkdir gmsh && cd gmsh      # create it wherever you want
wget gmsh.info/src/gmsh-4.9.5-source.tgz
tar zxvf gmsh-4.9.5-source.tgz
cd gmsh-4.9.5-source/
mkdir build && cd build
cmake -DENABLE_BUILD_DYNAMIC=1 -DENABLE_FLTK=0 .. # if issue with cgns, add option -DENABLE_CGNS=0 
make                       # and go get a coffee
```

### 3. (Optional) Install AGMG version 3.3.5 or later
AGMG must be compiled to get the .o files

```bash
mkdir agmg && cd agmg     # create it wherever you want
mv <path>/AGMG_3.3.5-aca.for.tar.gz .
tar -xvf AGMG_3.3.5-aca.for.tar.gz
cd AGMG_3.3.5-aca/Example_seq/
make                      # the .o files should be in the SRC/ directory
```

### 4. Use CMake to build fhhos4

With the `fhhos4` conda environment activated:

```bash
cd <path-to-fhhos4>
mkdir build && cd build
cmake -G Ninja -DCMAKE_BUILD_TYPE=Release ..
ninja
```

Ninja compiles at most 2 files in parallel by default (`-DMAX_COMPILE_JOBS=<n>` to change
it), as the compilation of `Program.cpp` takes ~4 GB of RAM. With Make (`cmake` without
`-G Ninja`), run `make -j2`.

Build options (`cmake -D<OPTION>=ON|OFF ..`):

| Option | Default | |
|---|---|---|
| `ENABLE_2D`, `ENABLE_3D` | ON | Dimensions to compile (`ENABLE_1D`: OFF) |
| `ENABLE_DG` | OFF | DG programs (`-discr dg`) |
| `ENABLE_FEM` | OFF | FEM programs (`-discr fem`) |
| `ENABLE_TESTS` | ON | Test suite (`fhhos4_tests`, run with `ctest`). Turning it off saves little time with Ninja, more with Make. |
| `ENABLE_GMSH` | ON | GMSH interface |
| `ENABLE_CGAL` | ON | Polygonal meshes (2D only) |
| `ENABLE_PROGRAM` | ON | The program `fhhos4` and its tests. OFF: the library alone, which only needs Eigen |
| `ENABLE_AGMG` | OFF | AGMG solver (`-s agmg`, see below) |
| `ENABLE_OPENMP` | ON | Shared-memory parallelism |
| `ENABLE_NATIVE_ARCH` | OFF | `-march=native`: optimized for the CPU of the build machine, may not run on others |
| `ENABLE_PYTHON` | OFF | The Python module of the library (see `README_library.md`); needs Python, nanobind, NumPy and SciPy |
| `SMALL_INDEX` | ON | `int` indices in the sparse matrices (`intmax_t` otherwise) |

If you built GMSH from source instead of using the conda package, pass its location explicitly:

```bash
cmake -DCMAKE_BUILD_TYPE=Release -DGMSH_API=<path>/gmsh/gmsh-4.9.5-source/api -DGMSH_LIB=<path>/gmsh/gmsh-4.9.5-source/build/libgmsh.so ..
```

And to enable AGMG:

```bash
cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_AGMG=ON -DAGMG_DIR=<path>/agmg/AGMG_3.3.5-aca/SRC ..
```

The build also produces the library `lib/libfhhos4.so` (see [README_library.md](README_library.md)).

### 5. (Optional) Install

```bash
cmake --install . --prefix <prefix>
```

installs the program (`bin/fhhos4`), the meshes (`share/fhhos4/meshes`), and the library with its headers and CMake
package (`lib/`, `include/`). The program also runs from the build directory. `--component library` or `--component
program` installs one part only.

### 6. Launch help command to view arguments and examples

```bash
./bin/fhhos4 -h
```

The program reads and writes:
- the meshes (`-geo <name>`): `data/meshes`, found through `share/fhhos4/meshes` next to the directory of the executable
  (in the build directory, a link to `data/meshes`; in the installation, a copy). The environment variable
  `FHHOS4_DATA_DIR` replaces `share/fhhos4`.
- the exported files: `./out` in the working directory, or `-o <directory>`.
- the cache of the GMSH meshes: `$XDG_CACHE_HOME/fhhos4`, by default `~/.cache/fhhos4` (`-no-cache` to bypass it).


## Troubleshooting

### Could NOT find CGAL / Eigen3

Make sure the `fhhos4` conda environment is activated before running `cmake` — Eigen and CGAL are resolved via that environment, not vendored in the repository.

### Could NOT find GMP or MPFR

CGAL depends on GMP and MPFR; the conda `cgal` package normally pulls these in automatically. If CMake still reports them missing, install them explicitly into the environment:

```bash
conda install -n fhhos4 gmp mpfr
```
