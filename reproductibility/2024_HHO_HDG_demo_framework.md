# Homogeneous multigrid for hybrid discretizations: application to HHO methods
D. A. Di Pietro, Z. Dong, G. Kanschat, P. Matalon, A. Rupp, *Numerical Methods for Partial Differential Equations*, 2025

The experiments can be reproduced with release 1.2.1.

### Square

#### Injection operator $I_\ell^1$
This injection operator is enabled by the parameters ``-prolong 2 -disable-hor``.

```bash
-pb diff -geo square -mesh stri -mesher inhouse -s mg -cs s -tol 1e-6 -smoothers ags -cycle {V,1,1|V,2,2} -prolong 2 -disable-hor -k {1|2|3} -n {32|64|128|256|512}
```

#### Injection operator $I_\ell^2$
This injection operator is enabled by the parameters ``-prolong 1 -disable-hor``.

```bash
-pb diff -geo square -mesh stri -mesher inhouse -s mg -cs s -tol 1e-6 -smoothers ags -cycle {V,1,1|V,2,2} -prolong 1 -disable-hor -k {1|2|3} -n {32|64|128|256|512}
```

#### Injection operator $I_\ell^3$
This injection operator is enabled by the parameter ``-prolong 1``.

```bash
-pb diff -geo square -mesh stri -mesher inhouse -s mg -cs s -tol 1e-6 -smoothers ags -cycle {V,1,1|V,2,2} -prolong 1 -k {1|2|3} -n {32|64|128|256|512}
```

#### As a preconditioner, with injection operator $I_\ell^1$

```bash
-pb diff -geo square -mesh stri -mesher inhouse -s cg -prec mg -cs s -tol 1e-6 -smoothers ags -cycle V,1,1 -prolong 2 -disable-hor -k {1|2|3} -n {32|64|128|256|512}
```

### L-shape

#### Injection operator $I_\ell^1$

```bash
-pb diff -geo L_shape -mesh tri -mesher gmsh -s mg -cs r -tol 1e-6 -smoothers ags -cycle {V,1,1|V,2,2} -prolong 2 -disable-hor -k {1|2|3} -n {16|32|64|128|256}
```

#### Injection operator $I_\ell^2$

```bash
-pb diff -geo L_shape -mesh tri -mesher gmsh -s mg -cs r -tol 1e-6 -smoothers ags -cycle {V,1,1|V,2,2} -prolong 1 -disable-hor -k {1|2|3} -n {16|32|64|128|256}
```

#### Injection operator $I_\ell^3$

```bash
-pb diff -geo L_shape -mesh tri -mesher gmsh -s mg -cs r -tol 1e-6 -smoothers ags -cycle {V,1,1|V,2,2} -prolong 1 -k {1|2|3} -n {16|32|64|128|256}
```

### Cube

```bash
-pb diff -geo cube -mesh stetra -mesher inhouse -s mg -cs s -tol 1e-6 -smoothers ags -cycle {V,1,1|V,2,2} -prolong 1 -k {1|2|3} -n {8|16|32}
```