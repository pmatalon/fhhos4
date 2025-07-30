# High-order multigrid strategies for HHO discretizations of elliptic equations
D. A. Di Pietro, P. Matalon, P. Mycek, U. Rüde, *Numerical Linear Algebra with Applications*, 2022


The experiments can be reproduced with release 1.0.

In the following command lines, the strategies are configured by the parameter ``-hp-config {1|2|3|4}``:
```bash
-hp-config 1 # h-only
-hp-config 2 # p-h
-hp-config 3 # p-h*
-hp-config 4 # hp-h
```
### Section 3.4.1 (effect of basis normalization)
Execute the following command to show that the multigrid method diverges when orthonormal bases are used with local refinement, even for low degree and small problem size:
```bash
> fhhos4 -geo square4quadrants_tri_localref -no-cache -tc square -cs r -k 1 -n 32 -e-ogb 3
```
The parameter ``-e-ogb 3`` orthonormalizes the element bases. Replace it with ``-e-ogb 1`` (orthogonalization without normalization) to make the multigrid method converge properly.

Note that divergence occurs if the refinement ratio is >= 1e5. To see it, edit the file ``data/meshes/2D/square4quadrants_tri_localref.geo`` to play with the refinement ratio:
```bash
h_center = h*1e-5; # divergence
h_center = h*1e-4; # convergence
```

It is interesting to see that normalization does not seem to affect the conditioning of the matrix, and that using a direct solver seems to work. Indeed, execute the following command and check the L2-error obtained after the system is solved by a Cholesky factorization.
```bash
> fhhos4 -geo square4quadrants_tri_localref -no-cache -tc square -cs r -k 1 -n 256 -e-ogb 3 -s ch
```
Replace the last parameters with ``-e-ogb 1 -s mg`` to check that you obtain the same L2-error using the multigrid and without normalization.

### Figure 3 (advantages of high-order, smooth solution)
```bash
> fhhos4 -geo square -mesh cart -cs r -tol 1e-12 -s fcgmg -hp-config 2 -k 2 -n 1024
> fhhos4 -geo square -mesh cart -cs r -tol 1e-12 -s fcgmg -hp-config 2 -k 3 -n 512
> fhhos4 -geo square -mesh cart -cs r -tol 1e-12 -s fcgmg -hp-config 2 -k 4 -n 256
> fhhos4 -geo square -mesh cart -cs r -tol 1e-12 -s fcgmg -hp-config 2 -k 5 -n 128
```
### Figure 5 (equivalent test, non-smooth solution)
```bash
> fhhos4 -geo square4quadrants_tri_localref -tc kellogg -no-cache -cs r -s fcgmg -hp-config 1 -e-ogb 0 -f-ogb 0 -e-basis monomials -f-basis monomials -tol 1e-3 -k 1 -n 256
> fhhos4 -geo square4quadrants_tri_localref -tc kellogg -no-cache -cs r -s fcgmg -hp-config 1 -e-ogb 0 -f-ogb 0 -e-basis monomials -f-basis monomials -tol 1e-3 -k 2 -n 128
> fhhos4 -geo square4quadrants_tri_localref -tc kellogg -no-cache -cs r -s fcgmg -hp-config 1 -e-ogb 0 -f-ogb 0 -e-basis monomials -f-basis monomials -tol 1e-3 -k 3 -n 64 
> fhhos4 -geo square4quadrants_tri_localref -tc kellogg -no-cache -cs r -s fcgmg -hp-config 1 -e-ogb 0 -f-ogb 0 -e-basis monomials -f-basis monomials -tol 1e-3 -k 4 -n 64 
> fhhos4 -geo square4quadrants_tri_localref -tc kellogg -no-cache -cs r -s fcgmg -hp-config 1 -e-ogb 0 -f-ogb 0 -e-basis monomials -f-basis monomials -tol 1e-3 -k 5 -n 32 
```
### Figures 7 and 8 (square, Cart. mesh)
Set ``-s mg`` for Fig. 6 and ``-s fcgmg`` for Fig. 7.
```bash
> fhhos4 -geo square -mesh cart -cs r -k 5 -n 128 -tol 1e-12 -s {mg|fcgmg} -hp-config {1|2|3|4}
```
### Figure 9 and 10 (square, tri. mesh)
Set ``-s mg`` for Fig. 8 and ``-s fcgmg`` for Fig. 9.
```bash
> fhhos4 -geo square -cs r -k 5 -n 64 -tol 1e-12 -s {mg|fcgmg} -hp-config {1|2|3|4}
```
### Figure 11 and 12 (Kellogg)
Set ``-s mg`` for Fig. 10 and ``-s fcgmg`` for Fig. 11.
```bash
> fhhos4 -geo square4quadrants_tri_localref -tc kellogg -no-cache -cs r -k 5 -n 256 -tol 1e-12 -s {mg|fcgmg} -hp-config {1|2|3|4}
```
### Figure 13 and 14 (cube, Cart. mesh)
Set ``-s mg`` for Fig. 12 and ``-s fcgmg`` for Fig. 13.
```bash
> fhhos4 -geo cube -k 3 -n 64 -mesh cart -mesher inhouse -cs s -tol 1e-12 -s {mg|fcgmg} -hp-config {1|2|3|4}
```
### Figure 15 and 16 (cube, tet. mesh)
Set ``-s mg`` for Fig. 14 and ``-s fcgmg`` for Fig. 15.
```bash
> fhhos4 -geo cube -k 3 -n 32 -mesh stetra -mesher inhouse -cs s -tol 1e-12 -s {mg|fcgmg} -hp-config {1|2|3|4}
```
### Figure 17 and 18 (Geneva wheel)
Set ``-s mg`` for Fig. 16 and ``-s fcgmg`` for Fig. 17.
```bash
> fhhos4 -geo genevawheel.geo -tc default -k 2 -n 32 -cs m -tol 1e-12 -s {mg|fcgmg} -hp-config {1|2|3|4}
```

