###############################################################################################
#### 1. Install CGAL version 5.0 or later (see https://doc.cgal.org/latest/Manual/usage.html)

     ## install with the package manager

> sudo apt-get install libcgal-dev

     ## or install with Spack (Spack installation steps here: https://spack-tutorial.readthedocs.io/en/latest/tutorial_basics.html)

> spack install cgal
> spack load cgal 

###############################################################################################
#### 2. Install GMSH version 4.6 or later (see http://gmsh.info/#Download)

> mkdir gmsh && cd gmsh      # create it wherever you want
> wget gmsh.info/src/gmsh-4.6.0-source.tgz
> tar zxvf gmsh-4.6.0-source.tgz
> cd gmsh-4.6.0-source/
> mkdir build && cd build
> cmake -DENABLE_BUILD_DYNAMIC=1 .. # if issue with cgns, add option -DENABLE_CGNS=0 
> make

###############################################################################################
#### 3. Use CMake to build dghho

> cd <path/to/dghho>
> mkdir build && cd build
> cmake -DCMAKE_BUILD_TYPE=Release ..    # add -G "Unix Makefiles" if necessary
> make

###############################################################################################
#### 4. Launch help command to view arguments and examples

> ./bin/dghho -h