
# An $h$-multigrid method for Hybrid High-Order discretizations
D. A. Di Pietro, F. Hülsemann, P. Matalon, P. Mycek, U. Rüde, D. Ruiz, *SIAM Journal on Scientific Computing*, 2021

The experiments can be reproduced with release 1.0.
### Figure 4.1
```bash
# 2D cart
> fhhos4 -geo square -mesh cart -mesher inhouse -s mg -cycle V,1,1 -k 0 -n {32|64|128|256|512|1024}
> fhhos4 -geo square -mesh cart -mesher inhouse -s mg -cycle V,1,1 -k 1 -n {32|64|128|256|512|1024}
> fhhos4 -geo square -mesh cart -mesher inhouse -s mg -cycle V,1,1 -k 2 -n {32|64|128|256|512|1024}
> fhhos4 -geo square -mesh cart -mesher inhouse -s mg -cycle V,1,1 -k 3 -n {32|64|128|256|512|1024}
# 2D tri
> fhhos4 -geo square -mesh stri -mesher inhouse -s mg -cycle V,1,1 -k 0 -n {32|64|128|256|512|1024|2048}
> fhhos4 -geo square -mesh stri -mesher inhouse -s mg -cycle V,1,1 -k 1 -n {32|64|128|256|512|1024}
> fhhos4 -geo square -mesh stri -mesher inhouse -s mg -cycle V,1,1 -k 2 -n {32|64|128|256|512}
> fhhos4 -geo square -mesh stri -mesher inhouse -s mg -cycle V,1,1 -k 3 -n {32|64|128|256|512}
# 3D cart
> fhhos4 -geo cube -mesh cart -mesher inhouse -s mg -cycle V,1,1 -k 0 -n {16|32|64|128}
> fhhos4 -geo cube -mesh cart -mesher inhouse -s mg -cycle V,1,1 -k 1 -n {16|32|64|128}
> fhhos4 -geo cube -mesh cart -mesher inhouse -s mg -cycle V,1,1 -k 2 -n {8|16|32|64}
> fhhos4 -geo cube -mesh cart -mesher inhouse -s mg -cycle V,1,1 -k 3 -n {8|16|32|64}
# 3D tetra
> fhhos4 -geo cube -mesh stetra -mesher inhouse -s mg -cycle V,2,2 -k 0 -n 16 # diverging
> fhhos4 -geo cube -mesh stetra -mesher inhouse -s mg -cycle V,2,2 -k 1 -n {8|16|32|64}
> fhhos4 -geo cube -mesh stetra -mesher inhouse -s mg -cycle V,2,2 -k 2 -n {8|16|32}
> fhhos4 -geo cube -mesh stetra -mesher inhouse -s mg -cycle V,2,2 -k 3 -n {8|16|32}
```
### Figure 4.2
Replace the \* characters with numerical values in
```bash
> fhhos4 -geo square -mesh stri -mesher inhouse -k 1 -n 512 -s mg -cycle V,*,* 
```
### Figure 4.3
Replace the \* characters with numerical values in
```bash
> fhhos4 -geo square -mesh stetra -mesher inhouse -k 1 -n 32 -s mg -cycle V,*,* 
```
### Figure 4.4
Important parameter: ``-smoothers bj23,bj23``.

Replace the \* characters with numerical values in
```bash
> fhhos4 -geo square -mesh stri -mesher inhouse -k 1 -n 512 -s mg -smoothers bj23,bj23 -cycle V,*,* 
```
### Figure 4.5
```bash
# 2D cart
> fhhos4 -geo square -mesh cart -mesher inhouse -s mg -cycle V,0,3 -k 0 -n {32|64|128|256|512|1024}
> fhhos4 -geo square -mesh cart -mesher inhouse -s mg -cycle V,0,3 -k 1 -n {32|64|128|256|512|1024}
> fhhos4 -geo square -mesh cart -mesher inhouse -s mg -cycle V,0,3 -k 2 -n {32|64|128|256|512|1024}
> fhhos4 -geo square -mesh cart -mesher inhouse -s mg -cycle V,0,3 -k 3 -n {32|64|128|256|512|1024}
# 2D tri
> fhhos4 -geo square -mesh stri -mesher inhouse -s mg -cycle V,0,3 -k 0 -n {32|64|128|256|512|1024|2048}
> fhhos4 -geo square -mesh stri -mesher inhouse -s mg -cycle V,0,3 -k 1 -n {32|64|128|256|512|1024}
> fhhos4 -geo square -mesh stri -mesher inhouse -s mg -cycle V,0,3 -k 2 -n {32|64|128|256|512}
> fhhos4 -geo square -mesh stri -mesher inhouse -s mg -cycle V,0,3 -k 3 -n {32|64|128|256|512}
# 3D cart
> fhhos4 -geo cube -mesh cart -mesher inhouse -s mg -cycle V,0,6 -k 0 -n {16|32|64|128}
> fhhos4 -geo cube -mesh cart -mesher inhouse -s mg -cycle V,0,6 -k 1 -n {16|32|64|128}
> fhhos4 -geo cube -mesh cart -mesher inhouse -s mg -cycle V,0,6 -k 2 -n {8|16|32|64}
> fhhos4 -geo cube -mesh cart -mesher inhouse -s mg -cycle V,0,6 -k 3 -n {8|16|32|64}
# 3D tetra
> fhhos4 -geo cube -mesh stetra -mesher inhouse -s mg -cycle V,0,6 -k 0 -n 16 # diverging
> fhhos4 -geo cube -mesh stetra -mesher inhouse -s mg -cycle V,0,6 -k 1 -n {8|16|32|64}
> fhhos4 -geo cube -mesh stetra -mesher inhouse -s mg -cycle V,0,6 -k 2 -n {8|16|32}
> fhhos4 -geo cube -mesh stetra -mesher inhouse -s mg -cycle V,0,6 -k 3 -n {8|16|32}
```
### Figure 4.7
```bash
> fhhos4 -geo square4quadrants -tc kellogg -mesh cart -mesher inhouse -s mg -cycle V,1,1 -k 0 -n {32|64|128|256|512|1024}
> fhhos4 -geo square4quadrants -tc kellogg -mesh cart -mesher inhouse -s mg -cycle V,1,1 -k 1 -n {32|64|128|256|512|1024}
> fhhos4 -geo square4quadrants -tc kellogg -mesh cart -mesher inhouse -s mg -cycle V,1,1 -k 2 -n {32|64|128|256|512}
> fhhos4 -geo square4quadrants -tc kellogg -mesh cart -mesher inhouse -s mg -cycle V,1,1 -k 3 -n {32|64|128|256|512}
```
### Figure 4.8
Important parameters: ``-prolong {1|2} [-disable-hor]``.
#### (a) V(0,3)
```bash
# cell k and injection
> fhhos4 -geo square -mesh stri -mesher inhouse -k 1 -s mg -cycle V,0,3 -prolong 2 -disable-hor -n {32|64|128|256|512|1024}
# cell k and average
> fhhos4 -geo square -mesh stri -mesher inhouse -k 1 -s mg -cycle V,0,3 -prolong 1 -disable-hor -n {32|64|128|256|512|1024}
# cell k+1 and injection
> fhhos4 -geo square -mesh stri -mesher inhouse -k 1 -s mg -cycle V,0,3 -prolong 2 -n {32|64|128|256|512|1024}
# cell k+1 and average (final algo)
> fhhos4 -geo square -mesh stri -mesher inhouse -k 1 -s mg -cycle V,0,3 -prolong 1 -n {32|64|128|256|512|1024}
```
#### (b) V(1,2)
```bash
# cell k and injection
> fhhos4 -geo square -mesh stri -mesher inhouse -k 1 -s mg -cycle V,1,2 -prolong 2 -disable-hor -n {32|64|128|256|512|1024}
# cell k and average
> fhhos4 -geo square -mesh stri -mesher inhouse -k 1 -s mg -cycle V,1,2 -prolong 1 -disable-hor -n {32|64|128|256|512|1024}
# cell k+1 and injection
> fhhos4 -geo square -mesh stri -mesher inhouse -k 1 -s mg -cycle V,1,2 -prolong 2 -n {32|64|128|256|512|1024}
# cell k+1 and average (final algo)
> fhhos4 -geo square -mesh stri -mesher inhouse -k 1 -s mg -cycle V,1,2 -prolong 1 -n {32|64|128|256|512|1024}
```
### Figure 4.9
Modify the value of the parameter ``-heterog`` from 1e0 to 1e8:
#### (a) Heterogeneous weighting
```bash
> fhhos4 -geo square4quadrants -mesh cart -mesher inhouse -n 64 -s mg -g 1 -cycle V,0,3 -k 0 -heterog {1e0-1e8}
> fhhos4 -geo square4quadrants -mesh cart -mesher inhouse -n 64 -s mg -g 1 -cycle V,0,3 -k 1 -heterog {1e0-1e8}
> fhhos4 -geo square4quadrants -mesh cart -mesher inhouse -n 64 -s mg -g 1 -cycle V,0,3 -k 2 -heterog {1e0-1e8}
> fhhos4 -geo square4quadrants -mesh cart -mesher inhouse -n 64 -s mg -g 1 -cycle V,0,3 -k 3 -heterog {1e0-1e8}
```
#### (b) Homogeneous weighting
Important parameter: ``-disable-heterog-weight``.
```bash
> fhhos4 -geo square4quadrants -mesh cart -mesher inhouse -n 64 -s mg -g 1 -cycle V,0,3 -disable-heterog-weight -k 0 -heterog {1e0-1e8}
> fhhos4 -geo square4quadrants -mesh cart -mesher inhouse -n 64 -s mg -g 1 -cycle V,0,3 -disable-heterog-weight -k 1 -heterog {1e0-1e8}
> fhhos4 -geo square4quadrants -mesh cart -mesher inhouse -n 64 -s mg -g 1 -cycle V,0,3 -disable-heterog-weight -k 2 -heterog {1e0-1e8}
> fhhos4 -geo square4quadrants -mesh cart -mesher inhouse -n 64 -s mg -g 1 -cycle V,0,3 -disable-heterog-weight -k 3 -heterog {1e0-1e8}
```
### Figure 4.10
Important parameter: ``-fcs {0|1}``.
#### (a) With face coarsening
```bash
> fhhos4 -geo square -mesh cart -mesher inhouse -s mg -cycle V,0,3 -cs s -fcs 1 -k 0 -n {32|64|128|256|512}
> fhhos4 -geo square -mesh cart -mesher inhouse -s mg -cycle V,0,3 -cs s -fcs 1 -k 1 -n {32|64|128|256|512}
> fhhos4 -geo square -mesh cart -mesher inhouse -s mg -cycle V,0,3 -cs s -fcs 1 -k 2 -n {32|64|128|256|512}
> fhhos4 -geo square -mesh cart -mesher inhouse -s mg -cycle V,0,3 -cs s -fcs 1 -k 3 -n {32|64|128|256|512}
```
#### (b) Without face coarsening
```bash
> fhhos4 -geo square -mesh cart -mesher inhouse -s mg -cycle V,0,3 -cs s -fcs 0 -k 0 -n {32|64|128|256|512}
> fhhos4 -geo square -mesh cart -mesher inhouse -s mg -cycle V,0,3 -cs s -fcs 0 -k 1 -n {32|64|128|256|512}
> fhhos4 -geo square -mesh cart -mesher inhouse -s mg -cycle V,0,3 -cs s -fcs 0 -k 2 -n {32|64|128|256|512}
> fhhos4 -geo square -mesh cart -mesher inhouse -s mg -cycle V,0,3 -cs s -fcs 0 -k 3 -n {32|64|128|256|512}
```
### Figure 4.11
Important parameter: ``-cs b``.
```bash
# Custom Bey's refinement, V(0,6)
> fhhos4 -geo cube -mesh tetra -mesher inhouse -k 1 -s mg -cs b -cycle V,0,6 -n {8|16|32|64}
# Custom Bey's refinement, V(0,8)
> fhhos4 -geo cube -mesh tetra -mesher inhouse -k 1 -s mg -cs b -cycle V,0,8 -n {8|16|32|64}
# Custom Bey's refinement, V(0,10)
> fhhos4 -geo cube -mesh tetra -mesher inhouse -k 1 -s mg -cs b -cycle V,0,10 -n {8|16|32|64}
# Cartesian tet. refinement, V(0,6)
> fhhos4 -geo cube -mesh stetra -mesher inhouse -k 1 -s mg -cs s -cycle V,0,6 -n {8|16|32|64}
```
### Figure 4.12
```bash
> fhhos4 -geo platewith4holes -s mg -cs b -cycle V,0,10 -k 0 -n {8|16|32} # diverges at n=32
> fhhos4 -geo platewith4holes -s mg -cs b -cycle V,0,10 -k 1 -n {8|16|32}
> fhhos4 -geo platewith4holes -s mg -cs b -cycle V,0,10 -k 2 -n {8|16|32}
> fhhos4 -geo platewith4holes -s mg -cs b -cycle V,0,10 -k 3 -n {8|16|32}
```
### Figure 4.15
Important parameter: ``-cs r``.
```bash
> fhhos4 -geo squarecircle -s mg -cs r -cycle V,0,3 -k 0 -n {32|64|128|256|512|1024}
> fhhos4 -geo squarecircle -s mg -cs r -cycle V,0,3 -k 1 -n {32|64|128|256|512|1024}
> fhhos4 -geo squarecircle -s mg -cs r -cycle V,0,3 -k 2 -n {32|64|128|256|512}
> fhhos4 -geo squarecircle -s mg -cs r -cycle V,0,3 -k 3 -n {32|64|128|256|512}
```
