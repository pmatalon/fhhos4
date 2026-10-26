# Test of the Python module fhhos4 compiled with MPI, on several processes (ctest PythonModuleMPI: mpiexec): the hybrid
# system of test_module.py distributed by rows, through mpi4py, against SciPy's direct solver on the whole system. The
# cells of each process: horizontal strips; each face owned by the process of the smaller rank of its cells; the cells
# and faces renumbered in rank order; the ghost faces (and cells of A_FT) of each process numbered in reverse order (any
# order). Also with the process of rank 0 empty, and non-symmetric.
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
from test_module import HybridSystem, RelativeError


def Renumbering(ranks, nRanks):
    """The new numbers: by rank, then in the old order"""
    order = np.argsort(ranks, kind="stable")
    renumbering = np.empty(len(ranks), dtype=np.int64)
    renumbering[order] = np.arange(len(ranks))
    return renumbering


def DistributedSystem(emptyRank0, asymmetry=0.0):
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
    cellGhosts = ghostsOf(rowsFT.indices, cellBegin, nOwnCells)
    faceColumns = lambda rows: localColumns(rows, ghosts, faceBegin, nOwnFaces)

    solver = fhhos4.Solver()
    solver.Comm = comm
    solver.Dimension = 2
    solver.FaceDegree = 0
    solver.CellDegree = 0
    solver.Tolerance = 1e-10
    localA_TT = A_TT[ownCells][:, ownCells].tocsr()
    localA_TT.sort_indices()
    if asymmetry == 0:
        solver.SetupSymmetric(localA_TT, faceColumns(rowsTF), faceColumns(rowsFF), np.ones(nOwnCells), np.ones(nOwnFaces), ghosts)
    else:
        localA_FT = localColumns(rowsFT, cellGhosts, cellBegin, nOwnCells)
        solver.SetupNonsymmetric(localA_TT, faceColumns(rowsTF), localA_FT, faceColumns(rowsFF), np.ones(nOwnCells), np.ones(nOwnFaces), ghosts, cellGhosts)
    assert solver.IsSetUp() and solver.NumberOfLevels() > 1
    x_T, x_F = np.zeros(nOwnCells), np.zeros(nOwnFaces)
    result = solver.Solve(b_T[ownCells].copy(), b_F[ownFaces].copy(), x_T, x_F)
    assert result.Converged and 0 < result.Iterations < 50, result
    allX = np.concatenate(comm.allgather(x_T) + comm.allgather(x_F))
    reference = spla.spsolve(sp.bmat([[A_TT, A_TF], [A_FT, A_FF]]).tocsc(), np.concatenate([b_T, b_F]))
    assert RelativeError(allX, reference) < 1e-8
    # An error on one process only is raised on all of them
    wrong = np.ones(nOwnFaces + (1 if rank == P - 1 else 0))
    try:
        solver.SetupSymmetric(localA_TT, faceColumns(rowsTF), faceColumns(rowsFF), np.ones(nOwnCells), wrong, ghosts)
        assert False, "a wrong interpolation of 1 accepted"
    except fhhos4.Error:
        pass


def test_distributed_system():
    DistributedSystem(emptyRank0=False)


def test_empty_process():
    """The process of rank 0 without cells nor faces (e.g. NGSolve's distributed meshes)"""
    DistributedSystem(emptyRank0=True)


def test_nonsymmetric_system():
    """SetupNonsymmetric(): A_FT with the ghost cells of its rows (cellGhosts)"""
    DistributedSystem(emptyRank0=False, asymmetry=0.3)
    DistributedSystem(emptyRank0=True, asymmetry=0.3)


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
