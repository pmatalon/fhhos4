# Tests of the Python module fhhos4 (library/python) on a hybrid system of degree 0 built here (Cartesian mesh of the
# unit square, random coefficients: HybridSystem of tests/UAMGLibraryTest.cpp), against SciPy's direct solver.
# Run by ctest (PythonModule), with the build's python directory in PYTHONPATH.
import gc
import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import fhhos4


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
        self.CellInterpOfOne = np.ones(nCells)
        self.FaceInterpOfOne = np.ones(nFaces)
        self.b_T = np.ones(nCells)  # unit source in the cells
        self.b_F = np.zeros(nFaces)
        self.b = -self.A_FT @ (invDiag @ self.b_T)  # right-hand side of the condensed system

    def Uncondensed(self):
        return sp.bmat([[self.A_TT, self.A_TF], [self.A_FT, self.A_FF]]).tocsc()


def ToySolver():
    solver = fhhos4.Solver()
    solver.Dimension = 2
    solver.FaceDegree = 0
    solver.CellDegree = 0
    solver.Tolerance = 1e-10
    return solver


def RelativeError(x, reference):
    return np.linalg.norm(x - reference) / np.linalg.norm(reference)


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
    # The condensed matrix computed by the library
    assert abs(solver.CondensedMatrix() - system.A).max() < 1e-12 * abs(system.A).max()


def test_uncondensed_nonsymmetric():
    system = HybridSystem(32, asymmetry=0.3)
    solver = ToySolver()
    solver.SetupNonsymmetric(system.A_TT, system.A_TF, system.A_FT, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne)
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
    # The solver owns its copy of A: the caller's may go
    del A
    gc.collect()
    x = solver.Solve(system.b)
    assert RelativeError(x, spla.spsolve(system.A.tocsc(), system.b)) < 1e-8
    # In place, from an initial guess, with A_FF
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
    solver.SetupSymmetricCondensed(system.A, system.A_TT, system.A_TF, system.CellInterpOfOne, system.FaceInterpOfOne)
    x = solver.Solve(system.b)
    other = HybridSystem(8)
    try:
        solver.SetupSymmetricCondensed(other.A, system.A_TT, system.A_TF, system.CellInterpOfOne, system.FaceInterpOfOne)
        assert False, "A of the wrong size accepted"
    except fhhos4.Error:
        pass
    gc.collect()
    assert RelativeError(solver.Solve(system.b), x) < 1e-12


def test_other_formats():
    # CSC and COO matrices (converted to CSR), int64 indices
    system = HybridSystem(32)
    solver = ToySolver()
    A_TF = system.A_TF.tocsc()
    A_TF.indices, A_TF.indptr = A_TF.indices.astype(np.int64), A_TF.indptr.astype(np.int64)
    solver.SetupSymmetric(system.A_TT.tocoo(), A_TF, system.A_FF.tocsc(), system.CellInterpOfOne, system.FaceInterpOfOne)
    x_T, x_F = np.zeros(system.b_T.size), np.zeros(system.b_F.size)
    assert solver.Solve(system.b_T, system.b_F, x_T, x_F).Converged


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
    solver = fhhos4.Solver()
    try:
        solver.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne)
        assert False, "Dimension not set accepted"
    except fhhos4.Error as e:
        assert "Dimension" in str(e)
    assert isinstance(fhhos4.Error("message"), RuntimeError)
    # Output vectors: float64 only, not converted (they would not be written)
    solver = ToySolver()
    solver.SetupSymmetric(system.A_TT, system.A_TF, system.A_FF, system.CellInterpOfOne, system.FaceInterpOfOne)
    try:
        solver.Solve(system.b_T, system.b_F, np.zeros(system.b_T.size, dtype=np.float32), np.zeros(system.b_F.size))
        assert False, "float32 output accepted"
    except TypeError:
        pass


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
