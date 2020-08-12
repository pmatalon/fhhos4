#### Install CGAL and its dependencies

>   sudo apt-get install libcgal-dev

#### Use CMake to build the makefile

>	cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DENABLE_GMSH=ON -DENABLE_AGMG=OFF

#### Compile

>	make

#### Launch help command to view arguments and examples

>	./bin/Release/dghho -h