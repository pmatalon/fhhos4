#### 1. Install CGAL (see https://doc.cgal.org/latest/Manual/usage.html)

     ## install with the package manager

> sudo apt-get install libcgal-dev

     ## or install with Spack (Spack installation steps here: https://spack-tutorial.readthedocs.io/en/latest/tutorial_basics.html)

> spack install cgal
> spack load cgal 

#### 2. Install GMSH (see http://gmsh.info/#Download)

> mkdir gmsh       # create it wherever you want
> cd gmsh
> wget gmsh.info/src/gmsh-4.6.0-source.tgz
> tar zxvf gmsh-4.6.0-source.tgz
> gmsh-4.6.0-source
> mkdir build
> cd build
> cmake -DENABLE_BUILD_DYNAMIC=1 .. # if issue with cgns, add option -DENABLE_CGNS=0 
> make

#### 3. Use CMake to build the makefile

> cd dghho
> cmake -DCMAKE_BUILD_TYPE=Release .    # add -G "Unix Makefiles" if necessary
> make

#### 4. Launch help command to view arguments and examples

> ./bin/Release/dghho -h