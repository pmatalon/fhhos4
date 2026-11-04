# Test of the Python module fhhos4 compiled with MPI, on several processes (ctest PythonModuleMPI: mpiexec): the hybrid
# system of test_module.py distributed by rows, through mpi4py, against SciPy's direct solver on the whole system. The
# cells of each process: horizontal strips; each face owned by the process of the smaller rank of its cells; the cells
# and faces renumbered in rank order; the ghost faces (and cells of A_FT) of each process numbered in reverse order (any
# order). Also with the process of rank 0 empty, non-symmetric, and the blocks numbered locally by LocalColumnsSym() and
# LocalColumnsNonSym() (the ghosts in increasing order). The condensed system (CondenseRHS(), Solve(b, x),
# SolveCellUnknowns(), ApplyCondensedMatrix()) too.
# Skipped (exit code 77) if the module is compiled without MPI, or without mpi4py.
import os
import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import fhhos4

if not fhhos4.WithMPI:
    print("The module is compiled without MPI: skipped")
    sys.exit(77)
try:
    from mpi4py import MPI
except ImportError:
    print("mpi4py is not installed: skipped")
    sys.exit(77)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from test_module import HybridSystem, RelativeError, ToSciPy


def Renumbering(ranks, nRanks):
    """The new numbers: by rank, then in the old order"""
    order = np.argsort(ranks, kind="stable")
    renumbering = np.empty(len(ranks), dtype=np.int64)
    renumbering[order] = np.arange(len(ranks))
    return renumbering


def DistributedSystem(emptyRank0, asymmetry=0.0, helpers=False):
    comm = MPI.COMM_WORLD
    rank, P = comm.Get_rank(), comm.Get_size()
    n = 48
    system = HybridSystem(n, asymmetry)
    nCells, nFaces = system.A_TT.shape[0], system.A_FF.shape[0]
    if emptyRank0:
        cellRank = np.array([1 + (T // n) * (P - 1) // n for T in range(nCells)])
    else:
        cellRank = np.array([(T // n) * P // n for T in range(nCells)])
    faceRank = np.full(nFaces, P)
    A_TF = system.A_TF.tocsr()
    for T in range(nCells):
        for F in A_TF.indices[A_TF.indptr[T]:A_TF.indptr[T + 1]]:
            faceRank[F] = min(faceRank[F], cellRank[T])
    newCell, newFace = Renumbering(cellRank, P), Renumbering(faceRank, P)
    # The system in the new numbering, and the rows of this process
    perm = lambda new, size: sp.csr_matrix((np.ones(size), (new, np.arange(size))), shape=(size, size))
    PT, PF = perm(newCell, nCells), perm(newFace, nFaces)
    A_TT, A_TF, A_FF = (PT @ system.A_TT @ PT.T).tocsr(), (PT @ system.A_TF @ PF.T).tocsr(), (PF @ system.A_FF @ PF.T).tocsr()
    A_FT = (PF @ system.A_FT @ PT.T).tocsr()
    b_T, b_F = PT @ system.b_T, PF @ system.b_F
    cellBegin, faceBegin = np.sum(cellRank < rank), np.sum(faceRank < rank)
    nOwnCells, nOwnFaces = np.sum(cellRank == rank), np.sum(faceRank == rank)
    ownCells, ownFaces = slice(cellBegin, cellBegin + nOwnCells), slice(faceBegin, faceBegin + nOwnFaces)
    rowsTF, rowsFF, rowsFT = A_TF[ownCells], A_FF[ownFaces], A_FT[ownFaces]

    def ghostsOf(referenced, begin, nOwn):
        """The entities of referenced not owned, in reverse order (any order)"""
        return np.array(sorted([e for e in np.unique(referenced) if not begin <= e < begin + nOwn], reverse=True), dtype=np.int64)

    def localColumns(rows, ghosts, begin, nOwn):
        """The columns of rows numbered locally: the owned entities, then the ghosts"""
        local = {e: e - begin for e in range(begin, begin + nOwn)}
        local.update({e: nOwn + g for g, e in enumerate(ghosts)})
        coo = rows.tocoo()
        M = sp.csr_matrix((coo.data, (coo.row, [local[c] for c in coo.col])), shape=(rows.shape[0], nOwn + len(ghosts)))
        M.sort_indices()
        return M

    ghosts = ghostsOf(np.concatenate([rowsTF.indices, rowsFF.indices]), faceBegin, nOwnFaces)
    cellGhostGlobalNumbers = ghostsOf(rowsFT.indices, cellBegin, nOwnCells)
    if helpers:
        ghosts, cellGhostGlobalNumbers = np.sort(ghosts), np.sort(cellGhostGlobalNumbers)
    faceColumns = lambda rows: localColumns(rows, ghosts, faceBegin, nOwnFaces)
    localA_TF, localA_FF, localA_FT = faceColumns(rowsTF), faceColumns(rowsFF), localColumns(rowsFT, cellGhostGlobalNumbers, cellBegin, nOwnCells)
    if helpers:
        tf, ft, ff = (M.sorted_indices() for M in (rowsTF, rowsFT, rowsFF))
        if asymmetry == 0:
            blocks = fhhos4.LocalColumnsSym(tf, ff, faceBegin, 1)
        else:
            blocks = fhhos4.LocalColumnsNonSym(tf, ft, ff, faceBegin, cellBegin, 1, 1)
        assert np.array_equal(blocks.FaceGhostGlobalNumbers, ghosts) and np.array_equal(blocks.CellGhostGlobalNumbers, cellGhostGlobalNumbers if asymmetry != 0 else [])
        assert (ToSciPy(blocks.A_TF) != localA_TF).nnz == 0 and (ToSciPy(blocks.A_FF) != localA_FF).nnz == 0
        if asymmetry != 0:
            assert (ToSciPy(blocks.A_FT) != localA_FT).nnz == 0
            localA_FT, cellGhostGlobalNumbers = blocks.A_FT, blocks.CellGhostGlobalNumbers
        localA_TF, localA_FF, ghosts = blocks.A_TF, blocks.A_FF, blocks.FaceGhostGlobalNumbers

    solver = fhhos4.Solver()
    solver.Comm = comm
    solver.Dimension = 2
    solver.FaceDegree = 0
    solver.CellDegree = 0
    solver.Tolerance = 1e-10
    localA_TT = A_TT[ownCells][:, ownCells].tocsr()
    localA_TT.sort_indices()
    if asymmetry == 0:
        solver.SetupSym(localA_TT, localA_TF, localA_FF, np.ones(nOwnCells), np.ones(nOwnFaces), ghosts)
    else:
        solver.SetupNonSym(localA_TT, localA_TF, localA_FT, localA_FF, np.ones(nOwnCells), np.ones(nOwnFaces), ghosts, cellGhostGlobalNumbers)
    assert solver.IsSetUp() and solver.NumberOfLevels() > 1
    x_T, x_F = np.zeros(nOwnCells), np.zeros(nOwnFaces)
    result = solver.Solve(b_T[ownCells].copy(), b_F[ownFaces].copy(), x_T, x_F)
    assert result.Converged and 0 < result.Iterations < 50, result
    allX = np.concatenate(comm.allgather(x_T) + comm.allgather(x_F))
    reference = spla.spsolve(sp.bmat([[A_TT, A_TF], [A_FT, A_FF]]).tocsc(), np.concatenate([b_T, b_F]))
    assert RelativeError(allX, reference) < 1e-8
    # The condensed system: the same solution; ApplyCondensedMatrix(), the rows of the process of the condensed matrix of
    # the whole system
    b, y_F, y_T = np.empty(nOwnFaces), np.zeros(nOwnFaces), np.empty(nOwnCells)
    solver.CondenseRHS(b_T[ownCells].copy(), b_F[ownFaces].copy(), b)
    assert solver.Solve(b, y_F).Converged
    solver.SolveCellUnknowns(b_T[ownCells].copy(), y_F, y_T)
    assert RelativeError(np.concatenate(comm.allgather(y_T) + comm.allgather(y_F)), reference) < 1e-8
    x, Ax = np.linspace(-1, 2, nFaces), np.empty(nOwnFaces)
    solver.ApplyCondensedMatrix(x[ownFaces].copy(), Ax)
    assert RelativeError(np.concatenate(comm.allgather(Ax)), (PF @ system.A @ PF.T) @ x) < 1e-13
    if asymmetry == 0 and P > 1:
        # A_FF lower triangular on the last process (here diagonal: one coupling added), read in full on several
        # processes: refused on all of them
        lowerA_FF = faceColumns(rowsFF).tolil()
        if rank == P - 1:
            lowerA_FF[1, 0] = -0.1
        lowerA_FF = lowerA_FF.tocsr()
        lowerA_FF.sort_indices()
        try:
            solver.SetupSym(localA_TT, faceColumns(rowsTF), lowerA_FF, np.ones(nOwnCells), np.ones(nOwnFaces), ghosts)
            assert False, "a lower triangular A_FF accepted on several processes"
        except fhhos4.Error as e:
            assert "lower triangular" in str(e), e
    # An error on one process only is raised on all of them
    wrong = np.ones(nOwnFaces + (1 if rank == P - 1 else 0))
    try:
        solver.SetupSym(localA_TT, faceColumns(rowsTF), faceColumns(rowsFF), np.ones(nOwnCells), wrong, ghosts)
        assert False, "a wrong interpolation of 1 accepted"
    except fhhos4.Error:
        pass


def test_distributed_system():
    DistributedSystem(emptyRank0=False)


def test_empty_process():
    """The process of rank 0 without cells nor faces (e.g. NGSolve's distributed meshes)"""
    DistributedSystem(emptyRank0=True)


def test_nonsymmetric_system():
    """SetupNonSym(): A_FT with the ghost cells of its rows (cellGhostGlobalNumbers)"""
    DistributedSystem(emptyRank0=False, asymmetry=0.3)
    DistributedSystem(emptyRank0=True, asymmetry=0.3)


def test_local_columns_helpers():
    """The blocks with global columns numbered locally by LocalColumnsSym() and LocalColumnsNonSym()"""
    for asymmetry in (0.0, 0.3):
        DistributedSystem(emptyRank0=False, asymmetry=asymmetry, helpers=True)
        DistributedSystem(emptyRank0=True, asymmetry=asymmetry, helpers=True)


if __name__ == "__main__":
    failures = 0
    for name, f in [(name, f) for name, f in list(globals().items()) if name.startswith("test_")]:
        try:
            f()
            ok = True
        except Exception as e:
            ok = False
            print(f"[FAIL] {name} (process {MPI.COMM_WORLD.Get_rank()}): {type(e).__name__}: {e}")
        ok = MPI.COMM_WORLD.allreduce(ok, op=MPI.LAND)
        failures += not ok
        if MPI.COMM_WORLD.Get_rank() == 0:
            print(f"[ {'OK' if ok else 'FAIL'} ] {name} on {MPI.COMM_WORLD.Get_size()} processes")
    sys.exit(1 if failures else 0)
