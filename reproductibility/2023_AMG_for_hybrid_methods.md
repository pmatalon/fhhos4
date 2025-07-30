# Algebraic multigrid preconditioner for statically condensed systems arising from lowest-order hybrid discretizations
D. A. Di Pietro, F. Hülsemann, P. Matalon, P. Mycek, U. Rüde, *SIAM Journal on Scientific Computing*, 2023


The experiments can be reproduced with release 1.0.
### Figure 4.2 and Table 4.2

Execute the following command line
```bash
> fhhos4 [[test_case]] [[solver]]
```
where ``[[test_case]]`` corresponds to
```bash
# Cube-cart
-geo cube -mesh cart -mesher inhouse -k 0 -n 128
# Cube-tet
-geo cube -mesh tetra -k 0 -n 64
# Complex-tet
-geo platewith4holes -k 0 -n 16
# Heterog1e8
-geo square4quadrants -k 0 -n 2048 -heterog 1e8
# Cube-cart-aniso100
-geo cube -mesh cart -mesher inhouse -k 0 -n 128 -aniso 100
# Cube-tet-aniso100
-geo cube -mesh tetra -k 0 -n 64 -aniso 20
```
and ``[[solver]]`` corresponds to
```bash
# U-AMG
-s fcguamg
# C-AMG
-s fcgaggregamg
# AGMG
-s agmg
```

### Figure 4.4
```bash
> fhhos4 -geo cube -mesh tetra -k 0 -n {16|32|64|128} [[solver]]
```
### Section 4.3.3 (Table 4.3 of the submitted version)
Details of the adaptive multiple coarsening strategy of U-AMG for the test
case Cube-tet.
```bash
> fhhos4 -geo cube -mesh tetra -k 0 -n 64 -s fcguamg
```
### Section 4.3.3 (Table 4.4 of the submitted version)
Details of the fixed double coarsening strategy of U-AMG for the test case
Cube-tet.
Important parameter: ``-cs dpa``.
```bash
> fhhos4 -geo cube -mesh tetra -k 0 -n 64 -s fcguamg -cs dpa
```
### Table 4.3
```bash
# U-AMG (multiple coarsening)
> fhhos4 -geo cube -mesh tetra -k 0 -n 64 -s fcguamg -cs mpa
# U-AMG (double coarsening)
> fhhos4 -geo cube -mesh tetra -k 0 -n 64 -s fcguamg -cs dpa
# C-AMG
> fhhos4 -geo cube -mesh tetra -k 0 -n 64 -s fcgaggregamg
# AGMG
> fhhos4 -geo cube -mesh tetra -k 0 -n 64 -s agmg
```
### Table 4.4
Important parameter: ``-coarsening-prolong {3|4|5|6}``.
```bash
# Cube-tet
> fhhos4 -geo cube -mesh tetra -k 0 -n 64 -s fcguamg -coarsening-prolong 6
> fhhos4 -geo cube -mesh tetra -k 0 -n 64 -s fcguamg -coarsening-prolong 4
> fhhos4 -geo cube -mesh tetra -k 0 -n 64 -s fcguamg -coarsening-prolong 5
> fhhos4 -geo cube -mesh tetra -k 0 -n 64 -s fcguamg -coarsening-prolong 3
# Cube-cart-aniso100
> fhhos4 -geo cube -mesh cart -mesher inhouse -k 0 -n 128 -aniso 100 -s fcguamg -coarsening-prolong 6
> fhhos4 -geo cube -mesh cart -mesher inhouse -k 0 -n 128 -aniso 100 -s fcguamg -coarsening-prolong 4
> fhhos4 -geo cube -mesh cart -mesher inhouse -k 0 -n 128 -aniso 100 -s fcguamg -coarsening-prolong 5
> fhhos4 -geo cube -mesh cart -mesher inhouse -k 0 -n 128 -aniso 100 -s fcguamg -coarsening-prolong 3
```
### Figure 4.4
Same commands as for Table 4.4, also with the solver ``-s fcgaggregamg``.
