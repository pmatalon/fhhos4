# fhhos4

fhhos4 (pronounce 'phosphore') stands for Fast HHO Solver for Research. It is an open-source code that solves diffusion problems with the Hybrid High-Order (HHO) discretizations. 

Features:
- solves diffusion problems with piecewise constant tensors
- HHO and DG discretizations
- high-order solutions
- handles unstructured meshes
- binded to the GMSH mesher
- implements geometric and algebraic multigrid methods
- exports linear systems
- shared-memory parallelism
- many parameters (various choices of polynomial bases, multigrid cycles, smoothers, coarsening strategies...)

The numerical experiments of related papers can be reproduced following the instructions in the `reproducibility` folder.

# Install

This program depends on the following third-party software packages:
- Eigen
- CGAL (version 5.0 or later)
- GMSH (version 4.9.5 or later)

and optionally
- AGMG (version 3.3.5 or later).

Eigen and CGAL are provided through a conda environment (`conda/environment.yml`). GMSH is also available in that environment, but can alternatively be built from source if you need a specific version.

## 1. Create and activate the conda environment

```bash
conda env create -f conda/environment.yml
conda activate fhhos4
```

Keep this environment activated whenever you configure or build the project, so that CMake can find Eigen and CGAL automatically.

## 2. (Optional) Build GMSH from source

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

## 3. (Optional) Install AGMG version 3.3.5 or later
AGMG must be compiled to get the .o files

```bash
mkdir agmg && cd agmg     # create it wherever you want
mv <path>/AGMG_3.3.5-aca.for.tar.gz .
tar -xvf AGMG_3.3.5-aca.for.tar.gz
cd AGMG_3.3.5-aca/Example_seq/
make                      # the .o files should be in the SRC/ directory
```

## 4. Use CMake to build fhhos4

With the `fhhos4` conda environment activated:

```bash
cd <path-to-fhhos4>
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

If you built GMSH from source instead of using the conda package, pass its location explicitly:

```bash
cmake -DCMAKE_BUILD_TYPE=Release -DGMSH_API=<path>/gmsh/gmsh-4.9.5-source/api -DGMSH_LIB=<path>/gmsh/gmsh-4.9.5-source/build/libgmsh.so ..
```

And to enable AGMG:

```bash
cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_AGMG=ON -DAGMG_DIR=<path>/agmg/AGMG_3.3.5-aca/SRC ..
```

## 5. Launch help command to view arguments and examples

```bash
./bin/fhhos4 -h
```


# Troubleshooting

### Could NOT find CGAL / Eigen3

Make sure the `fhhos4` conda environment is activated before running `cmake` — Eigen and CGAL are resolved via that environment, not vendored in the repository.

### Could NOT find GMP or MPFR

CGAL depends on GMP and MPFR; the conda `cgal` package normally pulls these in automatically. If CMake still reports them missing, install them explicitly into the environment:

```bash
conda install -n fhhos4 gmp mpfr
```
