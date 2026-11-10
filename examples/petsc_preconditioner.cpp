// fhhos4 as a preconditioner of PETSc (-pc_type fhhos4, README_library.md "With PETSc"), in C++: the HHO
// discretization of -div(kappa grad u) = 1 on the unit square (n x n squares, kappa 1 or 10 on a checkerboard of 8 x 8
// squares, homogeneous Dirichlet conditions), with polynomials of degree k on the faces and l in the cells, assembled
// cell by cell into a MATNEST, and solved by PETSc's FCG preconditioned by fhhos4. The program links fhhos4::petsc and
// registers the preconditioner itself (PetscDLLibraryRegister_fhhos4_petsc()); the discretization and the mesh:
// HHO_diffusion.h. The same as petsc_preconditioner.py, its local matrix up to rounding: the same iterations, though
// not the same residuals (on this regular mesh, many couplings are equal, and the aggregation breaks their ties by the
// rounding errors).
//
//   mpiexec -n 4 ./petsc_preconditioner [-n 64] [-k 1] [-l 1] [-ksp_view]
#include <string>
#include <vector>
#include <petscksp.h>
#include <fhhos4_petsc.h>
#include "HHO_diffusion.h"

using RowMajorMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

int main(int argc, char** argv)
{
	PetscInt n = 64;
	PetscInt k = 1;
	PetscInt l, iterations;
	PetscMPIInt rank, size;
	int nT, nF;
	Mat A_TT, A_TF, A_FT, A_FF, A;
	Vec x, b, one;
	IS rows[2];
	MatNullSpace nearNullSpace;
	KSP ksp;
	PC pc;
	KSPConvergedReason reason;
	PetscReal residual, bNorm;

	PetscCall(PetscInitialize(&argc, &argv, nullptr, nullptr));
	PetscCall(PetscDLLibraryRegister_fhhos4_petsc()); // -pc_type fhhos4
	PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
	PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
	PetscCall(PetscOptionsGetInt(nullptr, nullptr, "-n", &n, nullptr));
	PetscCall(PetscOptionsGetInt(nullptr, nullptr, "-k", &k, nullptr)); // the degree on the faces
	l = k;
	PetscCall(PetscOptionsGetInt(nullptr, nullptr, "-l", &l, nullptr)); // the degree in the cells
	HHO::Matrix A_local = HHO::LocalMatrixReferenceSquare((int)k, (int)l, nT, nF);
	HHO::Mesh mesh{ n, size };
	auto [firstCell, lastCell] = mesh.CellOwnershipRange(rank);
	auto [firstFace, lastFace] = mesh.FaceOwnershipRange(rank);
	PetscInt nCells = lastCell - firstCell;
	PetscInt nFaces = lastFace - firstFace;

	// The four blocks, preallocated (non-zeros per row in the columns of the process, and in the others)
	PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, nT * nCells, nT * nCells, PETSC_DETERMINE, PETSC_DETERMINE, nT, nullptr, 0, nullptr, &A_TT));
	PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, nT * nCells, nF * nFaces, PETSC_DETERMINE, PETSC_DETERMINE, 4 * nF, nullptr, 4 * nF, nullptr, &A_TF));
	PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, nF * nFaces, nT * nCells, PETSC_DETERMINE, PETSC_DETERMINE, 2 * nT, nullptr, 2 * nT, nullptr, &A_FT));
	PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, nF * nFaces, nF * nFaces, PETSC_DETERMINE, PETSC_DETERMINE, 7 * nF, nullptr, 7 * nF, nullptr, &A_FF));

	// The local matrices of the cells of the process, added; the DoFs of the boundary faces (Dirichlet, homogeneous)
	// dropped
	for (PetscInt T = firstCell; T < lastCell; T++)
	{
		std::array<std::int64_t, 4> faces = mesh.Faces(T);
		std::vector<int> cellLocal, faceLocal;
		std::vector<PetscInt> cellDoFs, faceDoFs;
		for (int t = 0; t < nT; t++)
		{
			cellLocal.push_back(t);
			cellDoFs.push_back(nT * T + t);
		}
		for (int f = 0; f < 4; f++)
		{
			for (int m = 0; faces[f] >= 0 && m < nF; m++)
			{
				faceLocal.push_back(nT + f * nF + m);
				faceDoFs.push_back(nF * faces[f] + m);
			}
		}
		HHO::Matrix a = mesh.Kappa(T) * A_local;
		RowMajorMatrix TT = a(cellLocal, cellLocal);
		RowMajorMatrix TF = a(cellLocal, faceLocal);
		RowMajorMatrix FT = a(faceLocal, cellLocal);
		RowMajorMatrix FF = a(faceLocal, faceLocal);
		PetscInt cT = (PetscInt)cellDoFs.size();
		PetscInt cF = (PetscInt)faceDoFs.size();
		PetscCall(MatSetValues(A_TT, cT, cellDoFs.data(), cT, cellDoFs.data(), TT.data(), ADD_VALUES));
		PetscCall(MatSetValues(A_TF, cT, cellDoFs.data(), cF, faceDoFs.data(), TF.data(), ADD_VALUES));
		PetscCall(MatSetValues(A_FT, cF, faceDoFs.data(), cT, cellDoFs.data(), FT.data(), ADD_VALUES));
		PetscCall(MatSetValues(A_FF, cF, faceDoFs.data(), cF, faceDoFs.data(), FF.data(), ADD_VALUES));
	}
	for (Mat M : { A_TT, A_TF, A_FT, A_FF })
	{
		PetscCall(MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY));
		PetscCall(MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY));
	}
	Mat blocks[4] = { A_TT, A_TF, A_FT, A_FF };
	PetscCall(MatCreateNest(PETSC_COMM_WORLD, 2, nullptr, 2, nullptr, blocks, &A));
	PetscCall(MatSetOption(A, MAT_SYMMETRIC, PETSC_TRUE)); // the symmetric setup of fhhos4

	// The source 1: its integral on the first function of each cell. The interpolation of 1 on the DoFs, for fhhos4: 1 on
	// the first function of each cell and face, 0 on the others (the near null space; orthonormal in PETSc)
	PetscCall(MatCreateVecs(A, &x, &b));
	PetscCall(VecDuplicate(b, &one));
	PetscCall(VecSet(b, 0));
	PetscCall(VecSet(one, 0));
	PetscCall(MatNestGetISs(A, rows, nullptr)); // the cell and face DoFs in the vectors of the operator
	for (int part = 0; part < 2; part++)
	{
		const PetscInt* indices;
		PetscInt nDoFs;
		PetscInt blockSize = part == 0 ? nT : nF;
		PetscCall(ISGetLocalSize(rows[part], &nDoFs));
		PetscCall(ISGetIndices(rows[part], &indices));
		for (PetscInt p = 0; p < nDoFs; p += blockSize)
		{
			if (part == 0)
				PetscCall(VecSetValue(b, indices[p], mesh.Measure(firstCell + p / nT), INSERT_VALUES));
			PetscCall(VecSetValue(one, indices[p], 1, INSERT_VALUES));
		}
		PetscCall(ISRestoreIndices(rows[part], &indices));
	}
	for (Vec v : { b, one })
	{
		PetscCall(VecAssemblyBegin(v));
		PetscCall(VecAssemblyEnd(v));
	}
	PetscCall(VecNormalize(one, nullptr));
	PetscCall(MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_FALSE, 1, &one, &nearNullSpace));
	PetscCall(MatSetNearNullSpace(A, nearNullSpace));

	// PETSc's flexible CG, preconditioned by fhhos4 (one K-cycle per application), down to a relative residual of 1e-8
	PetscCall(PetscOptionsSetValue(nullptr, "-pc_fhhos4_dimension", "2"));
	PetscCall(PetscOptionsSetValue(nullptr, "-pc_fhhos4_face_degree", std::to_string(k).c_str()));
	PetscCall(PetscOptionsSetValue(nullptr, "-pc_fhhos4_cell_degree", std::to_string(l).c_str()));
	PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
	PetscCall(KSPSetOperators(ksp, A, A));
	PetscCall(KSPSetType(ksp, KSPFCG));
	PetscCall(KSPGetPC(ksp, &pc));
	PetscCall(PCSetType(pc, PCFHHOS4));
	PetscCall(KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED)); // the residual b - A x (by default, FCG's is preconditioned)
	PetscCall(KSPSetTolerances(ksp, 1e-8, PETSC_CURRENT, PETSC_CURRENT, PETSC_CURRENT));
	PetscCall(KSPSetFromOptions(ksp));
	PetscCall(KSPSolve(ksp, b, x));
	PetscCall(KSPGetIterationNumber(ksp, &iterations));
	PetscCall(KSPGetConvergedReason(ksp, &reason));
	PetscCall(KSPGetResidualNorm(ksp, &residual));
	PetscCall(VecNorm(b, NORM_2, &bNorm));
	PetscCall(PetscPrintf(PETSC_COMM_WORLD, "%" PetscInt_FMT " x %" PetscInt_FMT " cells, degrees %" PetscInt_FMT " (faces) and %" PetscInt_FMT " (cells), %d processes: fcg + fhhos4, converged reason %d, %" PetscInt_FMT " iterations, relative residual %.1e\n",
	                      n, n, k, l, size, (int)reason, iterations, (double)(residual / bNorm)));

	PetscCall(KSPDestroy(&ksp));
	PetscCall(MatNullSpaceDestroy(&nearNullSpace));
	for (Vec* v : { &x, &b, &one })
		PetscCall(VecDestroy(v));
	for (Mat* M : { &A, &A_TT, &A_TF, &A_FT, &A_FF })
		PetscCall(MatDestroy(M));
	PetscCall(PetscFinalize());
	return 0;
}
