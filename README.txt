#### Install CGAL and its dependencies (see https://doc.cgal.org/latest/Manual/usage.html)

## install with the package manager

>   sudo apt-get install libcgal-dev

## or install with Spack (Spack installation steps here: https://spack-tutorial.readthedocs.io/en/latest/tutorial_basics.html)

> spack install cgal
> module load boost-headers
> module load cgal

#### Use CMake to build the makefile

>	cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release

#### Compile

>	make

#### Launch help command to view arguments and examples

>	./bin/Release/dghho -h