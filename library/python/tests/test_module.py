# Tests of the Python module fhhos4 (library/python) on a hybrid system of degree 0 built here (Cartesian mesh of the
# unit square, random coefficients: HybridSystem of tests/UAMGLibraryTest.cpp), against SciPy's direct solver. The
# SciPy matrices of the system: CSR with int32 indices, the format of the library, read in place.
# Run by ctest (PythonModule), with the build's python directory in PYTHONPATH.
import gc
import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import fhhos4

# The module compiled with MPI: the solvers of these tests on one process (Comm = MPI.COMM_SELF), MPI initialized by
# mpi4py (skipped without it)
if fhhos4.WithMPI:
    try:
        from mpi4py import MPI
    except ImportError:
        print("The module is compiled with MPI, and mpi4py is not installed: skipped")
        sys.exit(77)


def NewSolver():
    solver = fhhos4.Solver()
    if fhhos4.WithMPI:
        solver.Comm = MPI.COMM_SELF
    return solver


class HybridSystem:
    """Hybrid system of degree 0 on an n x n Cartesian mesh, interior faces only. asymmetry: A_FT != A_TF^T."""

    def __init__(self, n, asymmetry=0.0):
        nCells = n * n
        nVertical = (n - 1) * n
        nFaces = 2 * nVertical
        vertical = lambda i, j: j * (n - 1) + i                # between (i, j) and (i+1, j)
        horizontal = lambda i, j: nVertical + j * n + i       # between (i, j) and (i, j+1)
        rng = np.random.default_rng(42)
        TT, TF, FT, FF = [], [], [], []
        for j in range(n):
            for i in range(n):
                T = j * n + i
                a = rng.uniform(1, 10)
                TT.append((T, T, 4 * a))
                faces = []  # face, s (1: outflow)
                if i > 0: faces.append((vertical(i - 1, j), -1))
                if i < n - 1: faces.append((vertical(i, j), 1))
                if j > 0: faces.append((horizontal(i, j - 1), -1))
                if j < n - 1: faces.append((horizontal(i, j), 1))
                for F, s in faces:
                    TF.append((T, F, -a * (1 + asymmetry * s)))
                    FT.append((F, T, -a * (1 - asymmetry * s)))
                    FF.append((F, F, a))
        coo = lambda t, shape: sp.csr_matrix(([v for _, _, v in t], ([r for r, _, _ in t], [c for _, c, _ in t])), shape=shape)
        self.A_TT = coo(TT, (nCells, nCells))
        self.A_TF = coo(TF, (nCells, nFaces))
        self.A_FT = coo(FT, (nFaces, nCells))
        self.A_FF = coo(FF, (nFaces, nFaces))
        invDiag = sp.diags(1 / self.A_TT.diagonal())
        self.A = (self.A_FF - self.A_FT @ invDiag @ self.A_TF).tocsr()
        self.A.sort_indices()  # SciPy's products don't sort the columns: the CSR format of the library requires it
        self.CellInterpOfOne = np.ones(nCells)
        self.FaceInterpOfOne = np.ones(nFaces)
        self.b_T = np.ones(nCells)  # unit source in the cells
        self.b_F = np.zeros(nFaces)
        self.b = -self.A_FT @ (invDiag @ self.b_T)  # right-hand side of the condensed system

    def Uncondensed(self):
        return sp.bmat([[self.A_TT, self.A_TF], [self.A_FT, self.A_FF]]).tocsc()


def ToySolver():
    solver = NewSolver()
    solver.Dimension = 2
    solver.FaceDegree = 0
    solver.CellDegree = 0
    solver.Tolerance = 1e-10
    return solver


def RelativeError(x, reference):
    return np.linalg.norm(x - reference) / np.linalg.norm(reference)


def ToSciPy(M):
    """A copy of a fhhos4.CSRMatrix as a SciPy matrix"""
    return sp.csr_matrix((M.Values.copy(), M.ColumnIndices.copy(), M.RowPointers.copy()), shape=(M.Rows, M.Cols))


def test_uncondensed_symmetric():
    system = HybridSystem(32)
    solver = ToySolver()
    solver.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne)
    assert solver.IsSetUp() and solver.NumberOfLevels() > 1
    x_T, x_F = np.zeros(system.b_T.size), np.zeros(system.b_F.size)
    result = solver.Solve(system.b_T, system.b_F, x_T, x_F)
    assert result.Converged and 0 < result.Iterations < 50, result
    x = spla.spsolve(system.Uncondensed(), np.concatenate([system.b_T, system.b_F]))
    assert RelativeError(np.concatenate([x_T, x_F]), x) < 1e-8
    # The condensed matrix computed by the library (a copy)
    assert abs(ToSciPy(solver.CondensedMatrix()) - system.A).max() < 1e-12 * abs(system.A).max()


def test_uncondensed_nonsymmetric():
    system = HybridSystem(32, asymmetry=0.3)
    solver = ToySolver()
    A_FT = system.A_FT.copy()
    solver.SetupNonsymmetric(system.A_TT, system.A_TF, A_FT, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne)
    # The solves read A_FT: the solver keeps its arrays, the caller's object may go
    del A_FT
    gc.collect()
    x_T, x_F = np.zeros(system.b_T.size), np.zeros(system.b_F.size)
    result = solver.Solve(system.b_T, system.b_F, x_T, x_F)
    assert result.Converged, result
    x = spla.spsolve(system.Uncondensed(), np.concatenate([system.b_T, system.b_F]))
    assert RelativeError(np.concatenate([x_T, x_F]), x) < 1e-8


def test_condensed():
    system = HybridSystem(32)
    solver = ToySolver()
    A = system.A.copy()
    solver.SetupSymmetricCondensed(A, system.A_TT, system.A_TF, system.CellInterpOfOne, system.FaceInterpOfOne)
    # Read in place: the solver keeps the arrays of A, the caller's object may go
    del A
    gc.collect()
    x = solver.Solve(system.b)
    assert RelativeError(x, spla.spsolve(system.A.tocsc(), system.b)) < 1e-8
    # From an initial guess, with A_FF
    solver.SetupSymmetricCondensed(system.A, system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne)
    y = np.zeros(system.b.size)
    result = solver.Solve(system.b, y)
    assert result.Converged and RelativeError(y, x) < 1e-8


def test_condensed_nonsymmetric():
    system = HybridSystem(32, asymmetry=0.3)
    solver = ToySolver()
    solver.SetupNonsymmetricCondensed(system.A, system.A_TT, system.A_TF, system.A_FT, system.CellInterpOfOne, system.FaceInterpOfOne)
    x = solver.Solve(system.b)
    assert RelativeError(x, spla.spsolve(system.A.tocsc(), system.b)) < 1e-8


def test_failed_setup_keeps_the_previous_one():
    system = HybridSystem(32)
    solver = ToySolver()
    A = system.A.copy()
    solver.SetupSymmetricCondensed(A, system.A_TT, system.A_TF, system.CellInterpOfOne, system.FaceInterpOfOne)
    del A
    x = solver.Solve(system.b)
    other = HybridSystem(8)
    try:
        solver.SetupSymmetricCondensed(other.A, system.A_TT, system.A_TF, system.CellInterpOfOne, system.FaceInterpOfOne)
        assert False, "A of the wrong size accepted"
    except fhhos4.Error:
        pass
    gc.collect()
    assert RelativeError(solver.Solve(system.b), x) < 1e-12


def ExpectTypeError(f, what):
    try:
        f()
        assert False, what + " accepted"
    except TypeError:
        pass


def test_other_formats():
    # Read in place: CSR matrices of int32 indices (the index type of the library). Refused: the other formats and index
    # types, which the helpers convert (copies) or view (TransposedView(): symmetric CSC matrices).
    system = HybridSystem(32)
    solver = ToySolver()
    x_T, x_F = np.zeros(system.b_T.size), np.zeros(system.b_F.size)
    csc, coo = system.A_FF.tocsc(), system.A_FF.tocoo()
    int64 = system.A_TF.copy()
    int64.indices, int64.indptr = int64.indices.astype(np.int64), int64.indptr.astype(np.int64)
    for A_TF, A_FF, what in [(system.A_TF, csc, "a CSC matrix"), (system.A_TF, coo, "a COO matrix"), (int64, system.A_FF, "int64 indices")]:
        ExpectTypeError(lambda: solver.SetupSymmetric(system.A_TT, A_TF, A_FF, system.CellInterpOfOne, system.FaceInterpOfOne), what)
    assert not solver.IsSetUp()
    # Converted, or viewed
    solver.SetupSymmetric(system.A_TT, fhhos4.ToCSR(int64), fhhos4.ToCSR(csc), system.CellInterpOfOne, system.FaceInterpOfOne)
    reference = solver.Solve(system.b_T, system.b_F, x_T, x_F)
    y_T, y_F = np.zeros(system.b_T.size), np.zeros(system.b_F.size)
    solver.SetupSymmetric(system.A_TT, system.A_TF, fhhos4.TransposedView(csc), system.CellInterpOfOne, system.FaceInterpOfOne)
    assert solver.Solve(system.b_T, system.b_F, y_T, y_F).Iterations == reference.Iterations
    assert RelativeError(y_F, x_F) < 1e-12
    # A CSRView of NumPy arrays
    A = system.A_TF
    view = fhhos4.CSRView(A.shape[0], A.shape[1], A.indptr, A.indices, A.data)
    solver.SetupSymmetric(system.A_TT, view, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne)
    y_T, y_F = np.zeros(system.b_T.size), np.zeros(system.b_F.size)
    assert solver.Solve(system.b_T, system.b_F, y_T, y_F).Iterations == reference.Iterations


def test_blocks_of_the_global_matrix():
    # The global matrix in its own numbering: the hybrid system with Dirichlet DoFs, scrambled, with NGSolve's index
    # types (uint64 row pointers, int32 columns). ExtractSymmetricBlocks() and ExtractNonsymmetricBlocks() give the
    # blocks of the lists, Gather() and Scatter() the vectors.
    for asymmetry in [0.0, 0.3]:
        system = HybridSystem(32, asymmetry)
        nT, nF, nD = system.A_TT.shape[0], system.A_FF.shape[0], 40
        n = nT + nF + nD
        rng = np.random.default_rng(7)
        perm = rng.permutation(n)  # perm[i]: the global number of the i-th DoF of (cells, faces, Dirichlet)
        D = sp.random(nD, nT + nF, density=0.01, random_state=1, format="csr")  # couplings with the Dirichlet DoFs
        M = sp.bmat([[system.A_TT, system.A_TF, D[:, :nT].T], [system.A_FT, system.A_FF, D[:, nT:].T], [D[:, :nT], D[:, nT:], sp.eye(nD)]]).tocsr()
        P = sp.csr_matrix((np.ones(n), (perm, np.arange(n))), shape=(n, n))
        M = (P @ M @ P.T).tocsr()
        M.sort_indices()
        view = fhhos4.CSRView(n, n, M.indptr.astype(np.uint64), M.indices.astype(np.int32), M.data)
        cells = fhhos4.DoFLists(nT, np.arange(nT + 1), perm[:nT])
        faces = fhhos4.DoFLists(nF, np.arange(nF + 1), perm[nT:nT + nF])
        solver = ToySolver()
        if asymmetry == 0:
            blocks = fhhos4.ExtractSymmetricBlocks(view, cells, faces)
            assert abs(ToSciPy(blocks.A_TT) - sp.tril(system.A_TT)).max() == 0
            assert abs(ToSciPy(blocks.A_FF) - sp.tril(system.A_FF)).max() == 0
            assert blocks.A_FT.Rows == 0
            solver.SetupSymmetric(blocks.A_TT, blocks.A_TF, blocks.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne)
        else:
            blocks = fhhos4.ExtractNonsymmetricBlocks(view, cells, faces)
            assert abs(ToSciPy(blocks.A_FT) - system.A_FT).max() == 0
            solver.SetupNonsymmetric(blocks.A_TT, blocks.A_TF, blocks.A_FT, blocks.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne)
        assert abs(ToSciPy(blocks.A_TF) - system.A_TF).max() == 0
        # The global vectors: b (homogeneous Dirichlet conditions), and the solution x (the Dirichlet values unchanged)
        b = np.zeros(n)
        fhhos4.Scatter(system.b_T, system.b_F, cells, faces, b)
        b_T, b_F = np.zeros(nT), np.zeros(nF)
        fhhos4.Gather(b, cells, faces, b_T, b_F)
        assert (b_T == system.b_T).all() and (b_F == system.b_F).all()
        x_T, x_F = np.zeros(nT), np.zeros(nF)
        assert solver.Solve(b_T, b_F, x_T, x_F).Converged
        x = np.full(n, 2.0)
        fhhos4.Scatter(x_T, x_F, cells, faces, x)
        assert (x[perm[nT + nF:]] == 2.0).all()
        reference = spla.spsolve(system.Uncondensed(), np.concatenate([system.b_T, system.b_F]))
        assert RelativeError(x[perm[:nT + nF]], reference) < 1e-8
        # The arrays of a block keep it alive
        values = blocks.A_TF.Values
        del blocks
        gc.collect()
        assert (values == system.A_TF.data).all()


def test_parameters():
    solver = fhhos4.Solver()
    assert solver.Krylov == fhhos4.Krylov.Default and solver.Cycle == fhhos4.Cycle.K
    assert solver.PCoarsening == fhhos4.PCoarsening.Minus2 and solver.Tolerance == 1e-8
    solver.PCoarsening = fhhos4.PCoarsening.ToZero
    solver.PreSmoother = fhhos4.Smoother.BlockGaussSeidel
    assert solver.PCoarsening == fhhos4.PCoarsening.ToZero and solver.PreSmoother == fhhos4.Smoother.BlockGaussSeidel
    assert fhhos4.Version().startswith("1.")


def test_errors():
    system = HybridSystem(8)
    # Dimension not set
    solver = NewSolver()
    try:
        solver.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne)
        assert False, "Dimension not set accepted"
    except fhhos4.Error as e:
        assert "Dimension" in str(e)
    assert isinstance(fhhos4.Error("message"), RuntimeError)
    # Vectors: float64 only, not converted
    solver = ToySolver()
    solver.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne)
    ExpectTypeError(lambda: solver.Solve(system.b_T, system.b_F, np.zeros(system.b_T.size, dtype=np.float32), np.zeros(system.b_F.size)), "a float32 output")
    ExpectTypeError(lambda: solver.Solve(system.b_T.astype(np.float32), system.b_F, np.zeros(system.b_T.size), np.zeros(system.b_F.size)), "a float32 input")
    # Not in the CSR format: unsorted columns
    A = system.A_TF.copy()
    A.indices[[0, 1]] = A.indices[[1, 0]]
    try:
        solver.SetupSymmetric(system.A_TT, A, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne)
        assert False, "unsorted columns accepted"
    except fhhos4.Error as e:
        assert "CSR" in str(e)


if __name__ == "__main__":
    tests = [(name, f) for name, f in list(globals().items()) if name.startswith("test_")]
    failures = 0
    for name, f in tests:
        try:
            f()
            print(f"[ OK ] {name}")
        except Exception as e:
            failures += 1
            print(f"[FAIL] {name}: {type(e).__name__}: {e}")
    print(f"{len(tests) - failures}/{len(tests)} passed")
    sys.exit(1 if failures else 0)
