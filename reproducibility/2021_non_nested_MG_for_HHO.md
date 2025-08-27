# Towards robust, fast solutions of elliptic equations on complex domains through hybrid high‐order discretizations and non‐nested multigrid methods

D. A. Di Pietro, F. Hülsemann, P. Matalon, P. Mycek, U. Rüde, D. Ruiz, *International Journal for Numerical Methods in Engineering*, 2021

The experiments can be reproduced with release 1.0.

## Figure 6 
`-prolong 7` enables the non-nested prolongation with exact L2-projection.\
`-prolong 9` enables the non-nested prolongation with approximate L2-projection computed by subtriangulation.

```bash
# (a)
-geo square -mesh tri -s mg -cs m -prolong 7 -k {0|1|2|3} -n {32|64|128|256|512|1024}
# (b)
-geo square -mesh tri -s mg -cs m -prolong 9 -k {0|1|2|3} -n {32|64|128|256|512|1024}
```

## Figure 8

```bash
-geo cube -mesh tetra -s mg -cs m -prolong 9 -k {0|1|2|3} -n {8|16|32|64|128}
```

## Figure 9

`-cs n` enables the agglomeration coarsening with face collapsing.\
`-prolong 9 -subtri-meth wco` enables the approximate L2-projection with the optimal subtriangulation ensuring that there is no overlap with the coarse elements.\
`-prolong 9 -subtri-meth bary` enables the approximate L2-projection with barycentric subtriangulation.\
Note that the parameter `-subtri-meth` requires Release 1.3.

```bash
# (a)
-geo square -mesh tri -s mg -cs n -prolong 9 [-subtri-meth wco] -k {0|1|2|3} -n {32|64|128|256|512|1024}
# (b)
-geo square -mesh tri -s mg -cs n -prolong 9 -subtri-meth bary -k {0|1|2|3} -n {32|64|128|256|512|1024}
```

# Figure 10 (b)
