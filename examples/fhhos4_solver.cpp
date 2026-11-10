// fhhos4 as the solver (README_library.md "How to use"), without PETSc: the HHO discretization of HHO_diffusion.h
// assembled into the blocks of the uncondensed system, Eigen matrices in the format of the library, and solved by
// fhhos4 (its FCG preconditioned by the AMG). With MPI (FHHOS4_MPI, README_library.md "With MPI"), each process gives
// the rows of its own cells and faces, the columns of the faces numbered locally: its own faces, then the others met in
// its rows, its ghosts, whose global numbers it gives too. Nothing sends the contributions of a cell to the rows of the
// faces of other processes: each process also computes the local matrices of the other cells of its faces (of other
// processes), for the rows of its faces. The same system as petsc_preconditioner.cpp.
//
//   mpiexec -n 4 ./fhhos4_solver [-n 64] [-k 1] [-l 1]    (fhhos4 compiled without MPI: ./fhhos4_solver ...)
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <fhhos4>
#include "HHO_diffusion.h"
#if FHHOS4_MPI
#include <algorithm>
#include <map>
#endif

int main(int argc, char** argv)
{
	int rank = 0;
	int size = 1;
#if FHHOS4_MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
	int n = 64;
	int k = 1;
	int l = -1;
	for (int a = 1; a + 1 < argc; a += 2)
	{
		int value = std::atoi(argv[a + 1]);
		if (!std::strcmp(argv[a], "-n"))
			n = value;
		else if (!std::strcmp(argv[a], "-k"))
			k = value; // the degree on the faces
		else if (!std::strcmp(argv[a], "-l"))
			l = value; // the degree in the cells
	}
	l = l < 0 ? k : l;
	int nT, nF;
	HHO::Matrix A_local = HHO::LocalMatrixReferenceSquare(k, l, nT, nF);
	HHO::Mesh mesh{ n, size };
	auto [firstCell, lastCell] = mesh.CellOwnershipRange(rank);
	auto [firstFace, lastFace] = mesh.FaceOwnershipRange(rank);
	std::int64_t nCells = lastCell - firstCell;
	std::int64_t nFaces = lastFace - firstFace;

	// The cells whose local matrices it computes: its own; with MPI, then the other cells of its faces
	std::vector<std::int64_t> cells;
	for (std::int64_t T = firstCell; T < lastCell; T++)
		cells.push_back(T);
#if FHHOS4_MPI
	for (std::int64_t F = firstFace; F < lastFace; F++)
		for (std::int64_t T : mesh.Cells(F))
			if (T < firstCell || T >= lastCell)
				cells.push_back(T);
	std::sort(cells.begin() + nCells, cells.end());
	cells.erase(std::unique(cells.begin() + nCells, cells.end()), cells.end());

	// The local number of a face: its own faces first (from firstFace), then its ghosts, in the order met
	std::vector<std::int64_t> ghostGlobalNumbers; // the global number of the local face nFaces + i
	std::map<std::int64_t, std::int64_t> ghostGlobalToLocal;
	auto local = [&](std::int64_t F) -> std::int64_t
	{
		if (F >= firstFace && F < lastFace)
			return F - firstFace;
		auto [it, isNew] = ghostGlobalToLocal.insert({ F, nFaces + (std::int64_t)ghostGlobalNumbers.size() });
		if (isNew)
			ghostGlobalNumbers.push_back(F);
		return it->second;
	};
#else
	auto local = [](std::int64_t F) { return F; }; // one process: all the faces are its own
#endif

	// The local matrices of the cells, added to the rows of its cells (A_TT, A_TF) and of its faces (A_FF)
	std::vector<Eigen::Triplet<double>> TT, TF, FF;
	for (std::int64_t T : cells)
	{
		std::array<std::int64_t, 4> faces = mesh.Faces(T);
		HHO::Matrix a = mesh.Kappa(T) * A_local;
		bool ownCell = T >= firstCell && T < lastCell;
		std::int64_t t0 = (T - firstCell) * nT; // its first row of A_TT and A_TF
		if (ownCell)
		{
			for (int t = 0; t < nT; t++)
				for (int s = 0; s < nT; s++)
					TT.emplace_back(t0 + t, t0 + s, a(t, s));
		}
		for (int f = 0; f < 4; f++)
		{
			for (int m = 0; faces[f] >= 0 && m < nF; m++)
			{
				int dof = nT + f * nF + m; // in the local matrix
				if (ownCell)
				{
					for (int t = 0; t < nT; t++)
						TF.emplace_back(t0 + t, local(faces[f]) * nF + m, a(t, dof));
				}
				if (faces[f] >= firstFace && faces[f] < lastFace) // its face: the rows of A_FF
				{
					for (int g = 0; g < 4; g++)
						for (int p = 0; faces[g] >= 0 && p < nF; p++)
							FF.emplace_back((faces[f] - firstFace) * nF + m, local(faces[g]) * nF + p, a(dof, nT + g * nF + p));
				}
			}
		}
	}
#if FHHOS4_MPI
	std::int64_t nCols = (nFaces + (std::int64_t)ghostGlobalNumbers.size()) * nF; // its faces, then its ghosts
#else
	std::int64_t nCols = nFaces * nF;
#endif
	fhhos4::SparseMatrix A_TT(nCells * nT, nCells * nT);
	fhhos4::SparseMatrix A_TF(nCells * nT, nCols);
	fhhos4::SparseMatrix A_FF(nFaces * nF, nCols);
	A_TT.setFromTriplets(TT.begin(), TT.end());
	A_TF.setFromTriplets(TF.begin(), TF.end());
	A_FF.setFromTriplets(FF.begin(), FF.end());

	// The interpolations of 1, one coefficient per cell and face: that of its first function, 1. The source 1: its
	// integral on the first function of each cell (1): its measure.
	Eigen::VectorXd cellInterpOfOne = Eigen::VectorXd::Ones(nCells);
	Eigen::VectorXd faceInterpOfOne = Eigen::VectorXd::Ones(nFaces);
	Eigen::VectorXd b_T = Eigen::VectorXd::Zero(nCells * nT);
	Eigen::VectorXd b_F = Eigen::VectorXd::Zero(nFaces * nF);
	for (std::int64_t T = firstCell; T < lastCell; T++)
		b_T[(T - firstCell) * nT] = mesh.Measure(T);

	bool converged = false;
	try
	{
		fhhos4::Solver solver;
		solver.Dimension = 2;
		solver.FaceDegree = k;
		solver.CellDegree = l;
		solver.Tolerance = 1e-8;
#if FHHOS4_MPI
		solver.Comm = MPI_COMM_WORLD;
		solver.SetupSym(A_TT, A_TF, A_FF, cellInterpOfOne, faceInterpOfOne, ghostGlobalNumbers);
#else
		solver.SetupSym(A_TT, A_TF, A_FF, cellInterpOfOne, faceInterpOfOne);
#endif
		Eigen::VectorXd x_T(nCells * nT);
		Eigen::VectorXd x_F = Eigen::VectorXd::Zero(nFaces * nF);
		fhhos4::Result result = solver.Solve(b_T, b_F, x_T, x_F);
		converged = result.Converged;
		if (rank == 0)
			std::printf("%d x %d cells, degrees %d (faces) and %d (cells), %d processes: fhhos4, %s, %d iterations, relative residual %.1e\n",
			            n, n, k, l, size, result.Converged ? "converged" : "not converged", result.Iterations, result.RelativeResidual);
	}
	catch (const fhhos4::Error& e)
	{
		std::fprintf(stderr, "%s\n", e.what());
	}
#if FHHOS4_MPI
	MPI_Finalize();
#endif
	return converged ? EXIT_SUCCESS : EXIT_FAILURE;
}
