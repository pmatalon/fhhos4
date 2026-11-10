# fhhos4 as a preconditioner of PETSc (-pc_type fhhos4, README_library.md "With PETSc"), from Python: the HHO
# discretization of -div(kappa grad u) = 1 on the unit square (n x n squares, kappa 1 or 10 on a checkerboard of 8 x 8
# squares, homogeneous Dirichlet conditions), with polynomials of degree k on the faces and l in the cells, assembled
# cell by cell into a MATNEST with petsc4py, and solved by PETSc's FCG preconditioned by fhhos4. Needs petsc4py and
# NumPy, and fhhos4 built with -DENABLE_PETSC_PC=ON (and -DENABLE_MPI=ON on several processes) against the same PETSc.
# The same as petsc_preconditioner.cpp.
#
#   mpiexec -n 4 python petsc_preconditioner.py -dll_append <prefix>/lib/libfhhos4_petsc.so [-n 64] [-k 1] [-l 1] [-ksp_view]
import sys
import numpy as np
from numpy.polynomial import legendre
import petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc

comm = PETSc.COMM_WORLD
rank, size = comm.getRank(), comm.getSize()
options = PETSc.Options()
n = options.getInt("n", 64)
k = options.getInt("k", 1)  # the degree on the faces
l = options.getInt("l", k)  # the degree in the cells


def L(m, x):
    """The Legendre polynomial of degree m"""
    return legendre.legval(x, np.eye(m + 1)[m])


def dL(m, x):
    return legendre.legval(x, legendre.legder(np.eye(m + 1)[m]))


def Basis(degree):
    """The functions L_a(x) L_b(y) of degree a + b <= degree, by degree: hierarchical, the first one 1"""
    return [(a, d - a) for d in range(degree + 1) for a in range(d, -1, -1)]


def LocalMatrixReferenceSquare():
    """The HHO matrix of the reference square [-1, 1]^2: the cell DoFs, then those of its faces, left, right, bottom,
    top (L_m(s), m <= k, s the coordinate along the face: the same functions from both cells of a face). The
    reconstruction R in P^{k+1}: (grad R, grad w) = (grad u_T, grad w) + sum_F (u_F - u_T, grad w.n)_F; the
    stabilization: sum_F h_F^-1 |pi_F(u_F - u_T - (R - pi_T R))|_F^2, h_F = 2. In 2D, the matrix of every square."""
    cellBasis, reconstruction = Basis(l), Basis(k + 1)[1:]  # R without its constant: its gradient
    nT, nF = len(cellBasis), k + 1
    x, w = legendre.leggauss(k + l + 3)
    X, Y = np.meshgrid(x, x, indexing="ij")
    W = np.outer(w, w)
    value = lambda a, b, X, Y: L(a, X) * L(b, Y)
    grad = lambda a, b, X, Y: (dL(a, X) * L(b, Y), L(a, X) * dL(b, Y))
    faces = [(0, -1), (0, 1), (1, -1), (1, 1)]  # (the axis of the normal, its sign)
    onFace = lambda axis, side: (np.full_like(x, side), x) if axis == 0 else (x, np.full_like(x, side))
    G = [grad(a, b, X, Y) for a, b in reconstruction]
    K = np.array([[np.sum(W * (gi[0] * gj[0] + gi[1] * gj[1])) for gj in G] for gi in G])
    B = np.zeros((len(reconstruction), nT + 4 * nF))
    for i, (a, b) in enumerate(reconstruction):
        for j, (c, d) in enumerate(cellBasis):
            gj = grad(c, d, X, Y)
            B[i, j] = np.sum(W * (G[i][0] * gj[0] + G[i][1] * gj[1]))
        for f, (axis, side) in enumerate(faces):
            gn = side * grad(a, b, *onFace(axis, side))[axis]
            for j, (c, d) in enumerate(cellBasis):
                B[i, j] -= np.sum(w * value(c, d, *onFace(axis, side)) * gn)
            for m in range(nF):
                B[i, nT + f * nF + m] = np.sum(w * L(m, x) * gn)
    R = np.linalg.solve(K, B)  # the coefficients of R on the reconstruction basis
    A = B.T @ R
    high = [i for i, (a, b) in enumerate(reconstruction) if a + b > l]  # R - pi_T R: the functions of R not in P^l
    for f, (axis, side) in enumerate(faces):
        mass = np.array([np.sum(w * L(m, x) ** 2) for m in range(nF)])
        project = lambda functions: np.array([[np.sum(w * L(m, x) * value(a, b, *onFace(axis, side))) for a, b in functions] for m in range(nF)]) / mass[:, None]
        D = -project(reconstruction)[:, high] @ R[high]
        D[:, :nT] -= project(cellBasis)
        D[:, nT + f * nF:nT + (f + 1) * nF] += np.eye(nF)
        A += D.T @ (mass[:, None] * D) / 2
    return A, nT, nF


A_local, nT, nF = LocalMatrixReferenceSquare()  # nT, nF: the DoFs of a cell, of a face

# The mesh, seen by the assembly as an unstructured mesh: its cells and faces by their global numbers, in rank order
# (those of the process of rank 0 first, etc.), each cell with its faces. Here, the cells of each process are rows of the
# mesh (a horizontal strip), and each face is owned by the process of its lower or left cell; numbered row by row: the
# cells of a row (the cell (i, j) is j n + i), and its faces: the n - 1 vertical ones between its cells, then the n
# horizontal ones above it (none above the last row).
facesBefore = lambda j: j * (2 * n - 1) - (n if j == n else 0)  # the faces of the rows below j
firstRow, lastRow = (rank * n + size - 1) // size, ((rank + 1) * n + size - 1) // size  # its rows: j size // n == rank
myCells = range(firstRow * n, lastRow * n)
nCells, nFaces = len(myCells), facesBefore(lastRow) - facesBefore(firstRow)


def Faces(T):
    """The faces of the cell T, in the order of the local matrix: left, right, bottom, top; None on the boundary
    (Dirichlet: no unknowns)"""
    i, j = T % n, T // n
    return [facesBefore(j) + i - 1 if i > 0 else None, facesBefore(j) + i if i < n - 1 else None,
            facesBefore(j - 1) + n - 1 + i if j > 0 else None, facesBefore(j) + n - 1 + i if j < n - 1 else None]


Kappa = lambda T: 10 if (8 * (T % n) // n + 8 * (T // n) // n) % 2 else 1  # the diffusion coefficient of the cell T
Measure = lambda T: (1 / n) ** 2  # the measure of the cell T (its area)

# The four blocks, preallocated (non-zeros per row in the columns of the process, and in the others)
AIJ = lambda rows, cols, nnz: PETSc.Mat().createAIJ(((rows, None), (cols, None)), nnz=nnz, comm=comm)
A_TT, A_TF = AIJ(nT * nCells, nT * nCells, (nT, 0)), AIJ(nT * nCells, nF * nFaces, (4 * nF, 4 * nF))
A_FT, A_FF = AIJ(nF * nFaces, nT * nCells, (2 * nT, 2 * nT)), AIJ(nF * nFaces, nF * nFaces, (7 * nF, 7 * nF))

# The local matrices of the cells of the process, added; the DoFs of the boundary faces (Dirichlet, homogeneous) dropped
for T in myCells:
    faces = Faces(T)
    local = [nT + f * nF + m for f, F in enumerate(faces) if F is not None for m in range(nF)]
    cellDoFs = nT * T + np.arange(nT, dtype=PETSc.IntType)
    faceDoFs = [nF * F + m for F in faces if F is not None for m in range(nF)]
    a = Kappa(T) * A_local
    A_TT.setValues(cellDoFs, cellDoFs, a[:nT, :nT], addv=True)
    A_TF.setValues(cellDoFs, faceDoFs, a[np.ix_(range(nT), local)], addv=True)
    A_FT.setValues(faceDoFs, cellDoFs, a[np.ix_(local, range(nT))], addv=True)
    A_FF.setValues(faceDoFs, faceDoFs, a[np.ix_(local, local)], addv=True)
for M in (A_TT, A_TF, A_FT, A_FF):
    M.assemble()
A = PETSc.Mat().createNest([[A_TT, A_TF], [A_FT, A_FF]], comm=comm)
A.setOption(PETSc.Mat.Option.SYMMETRIC, True)  # the symmetric setup of fhhos4

x, b = A.createVecs()
cellsIS, facesIS = A.getNestISs()[0]  # the cell and face DoFs in the vectors of the operator
cells, faces = cellsIS.getIndices(), facesIS.getIndices()
b.zeroEntries()
b.setValues(cells[::nT], [Measure(T) for T in myCells])  # the source 1: its integral on the first function of each cell
b.assemble()

# The interpolation of 1 on the DoFs, for fhhos4: 1 on the first function of each cell and face, 0 on the others
one = b.duplicate()
one.zeroEntries()
one.setValues(np.concatenate([cells[::nT], faces[::nF]]), np.ones(nCells + nFaces))
one.assemble()
one.normalize()  # PETSc's near null spaces are orthonormal
A.setNearNullSpace(PETSc.NullSpace().create(vectors=[one], comm=comm))

# PETSc's flexible CG, preconditioned by fhhos4 (one K-cycle per application), down to a relative residual of 1e-8
options["pc_fhhos4_dimension"] = 2
options["pc_fhhos4_face_degree"] = k
options["pc_fhhos4_cell_degree"] = l
ksp = PETSc.KSP().create(comm)
ksp.setOperators(A)
ksp.setType("fcg")
ksp.getPC().setType("fhhos4")
ksp.setNormType(PETSc.KSP.NormType.UNPRECONDITIONED)  # the residual b - A x (by default, FCG's is preconditioned)
ksp.setTolerances(rtol=1e-8)
ksp.setFromOptions()
ksp.solve(b, x)
PETSc.Sys.Print(f"{n} x {n} cells, degrees {k} (faces) and {l} (cells), {size} processes: {ksp.getType()} + {ksp.getPC().getType()}, "
                f"converged reason {ksp.getConvergedReason()}, {ksp.getIterationNumber()} iterations, relative residual {ksp.getResidualNorm() / b.norm():.1e}")
