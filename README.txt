This program depends on the following third-party softwares:
- Eigen
- CGAL (version 5.0 or later)
- GMSH (version 4.9.5 or later)
Optional:
- AGMG (version 3.3.5 or later)

Eigen and CGAL are header-only libraries, they are shipped with this code and shall be compiled simultaneously with the program.
So without the need for a specific version of those libraries, you have nothing to do. You can then skip step 1.
Only GMSH is to be installed. To do so, follow step 2.
Finally, build the program following step 4.


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
#### 2. Install GMSH version 4.9.5 or later (see http://gmsh.info/#Download)

> mkdir gmsh && cd gmsh      # create it wherever you want
> wget gmsh.info/src/gmsh-4.9.5-source.tgz
> tar zxvf gmsh-4.9.5-source.tgz
> cd gmsh-4.9.5-source/
> mkdir build && cd build
> cmake -DENABLE_BUILD_DYNAMIC=1 -DENABLE_FLTK=0 ..       # if issue with cgns, add option -DENABLE_CGNS=0 
> make                       # and go get a coffee

###############################################################################################
#### 3. (Optional) Install AGMG version 3.3.5 or later
#### AGMG must be compiled to get the .o files

> mkdir agmg && cd agmg     # create it wherever you want
> mv <path>/AGMG_3.3.5-aca.for.tar.gz .
> tar -xvf AGMG_3.3.5-aca.for.tar.gz
> cd AGMG_3.3.5-aca/Example_seq/
> make                      # the .o files should be in the SRC/ directory

###############################################################################################
#### 4. Use CMake to build fhhos4

> cd <path-to-fhhos4>
> mkdir build && cd build
> cmake -DCMAKE_BUILD_TYPE=Release -DGMSH_API=<path>/gmsh/gmsh-9.5.0-source/api -DGMSH_LIB=<path>/gmsh/gmsh-4.9.5-source/build/libgmsh.so -DENABLE_AGMG=ON -DAGMG_DIR=<path>/agmg/AGMG_3.3.5-aca/SRC ..
> make

###############################################################################################
#### 4. Launch help command to view arguments and examples

> ./bin/fhhos4 -h





###############################################################################################
#### TROUBLESHOOTING

## 1

CMake Error at /usr/share/cmake-3.10/Modules/FindPackageHandleStandardArgs.cmake:137 (message):
  Could NOT find GMP (missing: GMP_LIBRARIES GMP_INCLUDE_DIR)
Call Stack (most recent call first):
  /usr/share/cmake-3.10/Modules/FindPackageHandleStandardArgs.cmake:378 (_FPHSA_FAILURE_MESSAGE)
  dependencies/cgal/CGAL-5.1/cmake/modules/FindGMP.cmake:53 (find_package_handle_standard_args)
  dependencies/cgal/CGAL-5.1/cmake/modules/CGAL_SetupGMP.cmake:24 (find_package)
  dependencies/cgal/CGAL-5.1/cmake/modules/CGAL_SetupCGALDependencies.cmake:41 (include)
  dependencies/cgal/CGAL-5.1/lib/cmake/CGAL/CGALConfig.cmake:128 (include)
  dependencies/cgal/CGAL-5.1/CGALConfig.cmake:6 (include)
  CMakeLists.txt:64 (find_package)

> sudo apt-get install libgmp-dev

## 2

CMake Error at /usr/share/cmake-3.10/Modules/FindPackageHandleStandardArgs.cmake:137 (message):
  Could NOT find MPFR (missing: MPFR_LIBRARIES MPFR_INCLUDE_DIR)
Call Stack (most recent call first):
  /usr/share/cmake-3.10/Modules/FindPackageHandleStandardArgs.cmake:378 (_FPHSA_FAILURE_MESSAGE)
  dependencies/cgal/CGAL-5.1/cmake/modules/FindMPFR.cmake:52 (find_package_handle_standard_args)
  dependencies/cgal/CGAL-5.1/cmake/modules/CGAL_SetupGMP.cmake:25 (find_package)
  dependencies/cgal/CGAL-5.1/cmake/modules/CGAL_SetupCGALDependencies.cmake:41 (include)
  dependencies/cgal/CGAL-5.1/lib/cmake/CGAL/CGALConfig.cmake:128 (include)
  dependencies/cgal/CGAL-5.1/CGALConfig.cmake:6 (include)
  CMakeLists.txt:64 (find_package)

> sudo apt-get install libmpfr-dev

and probably also needs

> sudo apt-get install libboost-all-dev

