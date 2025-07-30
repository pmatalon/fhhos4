# Iterative solution to the biharmonic equation in mixed form discretized by the Hybrid High-Order method
P. F. Antonietti, P. Matalon, M. Verani, *Computers & Mathematics with Applications*, 2023


The experiments can be reproduced with release 1.2.

### Convergence of the discrete normal derivative
Fig. 2

```bash
-pb diff -geo square -source sine -s ch -mesh {cart|tri} -normalder -k {0|1|2|3} -n {16|32|64|128|256|512}
```
<!-- 
### Well-posedness
```bash
# square cart
-pb bihar -geo square       -source sine -s ch -nbh-depth 8 -mesh cart -cs r -kc {0|1} -k {0|1|2|3|4} -n {8|16}
# polygonal(2)
-pb bihar -geo squarecircularhole -tc default -s ch -bihar-prec p -nbh-depth 8 -mesh poly -cs r -polymesh-bfc c -polymesh-n-pass 1 -n 16 -kc {0|1} -k {0|1|2|3} -ut
# polygonal(4) -> retry until you get a max of 4 boundary faces per element
-pb bihar -geo squarecircularhole -tc default -s ch -bihar-prec p -nbh-depth 8 -mesh poly -cs r -polymesh-bfc c -polymesh-n-pass 2 -n 32 -kc {0|1} -k {0|1|2|3} -ut
-pb bihar -geo magnetism2         -tc default -s ch -bihar-prec p -nbh-depth 8 -mesh poly -cs r -polymesh-bfc c -polymesh-n-pass 1 -kc 1 -k 5 -n 8 -export mesh -ut -fc-coplanar-tol 1e-5 -threads 0 -f-basis monomials -f-ogb 0
``` 
-->

### Convergence of the biharmonic scheme
Fig. 3, 4, 5
```bash
# Cartesian mesh
-pb bihar -geo square -source {exp|poly} -s ch -bihar-prec s -nbh-depth 8 -mesh cart -cs r -k {0|1|2|3} -n {16|32|64|128|256} -tol 1e-10
# Polygonal mesh
-pb bihar -geo square -source exp -s ch -bihar-prec s -nbh-depth 8 -bihar-prec-solver bicgstab -mesh poly -polymesh-init cart -polymesh-n-pass 1 -polymesh-fcs c -k {0|1|2|3} -n {16|32|64|128|256} -tol 1e-10
```

### Preconditioner convergence
Fig. 7
```bash
-pb bihar -geo square -source poly -s ch -mesh cart -cs r -k 1 -n 256 -tol 1e-14 -export iter -bihar-prec no #no preconditioner
-pb bihar -geo square -source poly -s ch -mesh cart -cs r -k 1 -n 256 -tol 1e-14 -export iter -nbh-depth {2|4|6|8|10} #with preconditioner
```

### Asymptotic behaviour
Tables 1, 2, 3
```bash
-pb bihar -geo square -source exp -mesh cart -cs r -s ch -nbh-depth 8 -tol 1e-8 -bihar-prec {s|no} -k {0,1,2,3} -n {32,64,128,256,512} 
-pb bihar -geo cube -source exp -mesh tetra -not-compute-errors -s fcguamg -hp-cs p_h -nbh-depth 2 -bihar-prec-solver bicgstab -tol 1e-8 -k 0 -bihar-prec {s|no} -n {8,16,32,64} 
```
<!-- -pb bihar -geo {disk|magnetism2} -tc default -mesh poly -polymesh-bfc c -polymesh-n-pass 2 -s ch -nbh-depth 8 -tol 1e-8 -k 0 -n {64,128,256,512,1024}  -->

<!-- 
### Heuristics
```bash
-pb bihar -geo square -source poly -s fcguamg -mesh tri -cs r -k 3 -n 32 -tol 1e-8 -bihar-prec p -nbh-depth 2 -export iter -opt2 {0|2} 
``` -->
