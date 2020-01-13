#### Use CMake to build the makefile

>	cmake -G "Unix Makefiles" -DENABLE_GMSH=ON -DENABLE_AGMG=OFF

#### Compile

>	make

#### Launch program

>	./bin/dghho -d 2 -discr hho -a s -n 8 -p 1

#### Help

>	./bin/dghho -h