This program depends on the following third-party softwares:
- Eigen
- CGAL (version 5.0 or later)
- GMSH (version 4.6 or later)

Eigen and CGAL are header-only libraries, they are shipped with this code and shall be compiled simultaneously with the program.
So without the need for a specific version of those libraries, you have nothing to do. You can then skip step 1.
Only GMSH is to be installed. To do so, follow step 2.
Finally, build the program following step 3.


###############################################################################################
#### 1. Install CGAL version 5.0 or later (see https://doc.cgal.org/latest/Manual/usage.html)

     ## install with the package manager

> sudo apt-get install libcgal-dev

     ## or install with Spack (Spack installation steps here: https://spack-tutorial.readthedocs.io/en/latest/tutorial_basics.html)

> spack install cgal
> spack load cgal 

     ## or install with the sources

> mkdir cgal && cd cgal      # create it wherever you want
> wget https://github.com/CGAL/cgal/releases/download/v5.1/CGAL-5.1.tar.xz
> tar xf CGAL-5.1.tar.xz     # the path to CGAL must then be added to the cmake command via -DCGAL_DIR=<path>/cgal/CGAL-5.1

###############################################################################################
#### 2. Install GMSH version 4.6 or later (see http://gmsh.info/#Download)

> mkdir gmsh && cd gmsh      # create it wherever you want
> wget gmsh.info/src/gmsh-4.6.0-source.tgz
> tar zxvf gmsh-4.6.0-source.tgz
> cd gmsh-4.6.0-source/
> mkdir build && cd build
> cmake -DENABLE_BUILD_DYNAMIC=1 ..       # if issue with cgns, add option -DENABLE_CGNS=0 
> make

###############################################################################################
#### 3. Use CMake to build dghho

> cd <path-to-dghho>
> mkdir build && cd build
> cmake -DCMAKE_BUILD_TYPE=Release -DGMSH_API=<path-to-gmsh>/api -DGMSH_LIB=<path-to-gmsh>/build/libgmsh.so ..
> make

###############################################################################################
#### 4. Launch help command to view arguments and examples

> ./bin/dghho -h