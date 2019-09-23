#!/bin/bash
./bin/dghho -d 2 -a s -s pcgmg -n 32   -p 1 > log/Poisson2D_pcgmg_n32_p1.txt
./bin/dghho -d 2 -a s -s pcgmg -n 64   -p 1 > log/Poisson2D_pcgmg_n64_p1.txt
./bin/dghho -d 2 -a s -s pcgmg -n 128  -p 1 > log/Poisson2D_pcgmg_n128_p1.txt
./bin/dghho -d 2 -a s -s pcgmg -n 256  -p 1 > log/Poisson2D_pcgmg_n256_p1.txt
./bin/dghho -d 2 -a s -s pcgmg -n 512  -p 1 > log/Poisson2D_pcgmg_n512_p1.txt
./bin/dghho -d 2 -a s -s pcgmg -n 1024 -p 1 > log/Poisson2D_pcgmg_n1024_p1.txt