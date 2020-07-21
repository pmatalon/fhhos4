#### Install CGAL and its dependencies

>   sudo apt-get install libcgal-dev

#### Use CMake to build the makefile

>	cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DENABLE_GMSH=ON -DENABLE_AGMG=OFF

#### Compile

>	make

#### Launch program

>	./bin/Release/dghho -d 2 -discr hho -mesh cart -n 8 -p 2

#### Help

>	./bin/dghho -h