// fhhos4 as a preconditioner of PETSc (-pc_type fhhos4, README_library.md "With PETSc"), in C++: the HHO
// discretization of -div(kappa grad u) = 1 on the unit square (n x n squares, kappa 1 or 10 on a checkerboard of 8 x 8
// squares, homogeneous Dirichlet conditions), with polynomials of degree k on the faces and l in the cells, assembled
// cell by cell into a MATNEST, and solved by PETSc's FCG preconditioned by fhhos4. The program links fhhos4::petsc and
// registers the preconditioner itself (PetscDLLibraryRegister_fhhos4_petsc()); Eigen computes the local matrix. The
// same as petsc_preconditioner.py, its local matrix up to rounding: the same iterations, though not the same residuals
// (on this regular mesh, many couplings are equal, and the aggregation breaks their ties by the rounding errors).
//
//   mpiexec -n 4 ./petsc_preconditioner [-n 64] [-k 1] [-l 1] [-ksp_view]
#include <cmath>
#include <string>
#include <utility>
#include <vector>
#include <Eigen/Dense>
#include <petscksp.h>
#include <fhhos4_petsc.h>

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using RowMajorMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using Function = std::pair<int, int>; // L_a(x) L_b(y)

// The Legendre polynomial of degree m at x, and its derivative
std::pair<double, double> Legendre(int m, double x)
{
	double L = 1, previous = 0, dL = 0, dPrevious = 0;
	for (int i = 0; i < m; i++)
	{
		double next = ((2 * i + 1) * x * L - i * previous) / (i + 1), dNext = dPrevious + (2 * i + 1) * L;
		previous = L;
		L = next;
		dPrevious = dL;
		dL = dNext;
	}
	return { L, dL };
}

// The Gauss-Legendre quadrature of q points on [-1, 1] (Golub-Welsch)
void GaussLegendre(int q, Vector& x, Vector& w)
{
	Matrix J = Matrix::Zero(q, q);
	for (int i = 1; i < q; i++)
		J(i, i - 1) = J(i - 1, i) = i / std::sqrt(4.0 * i * i - 1);
	Eigen::SelfAdjointEigenSolver<Matrix> eigen(J);
	x = eigen.eigenvalues();
	w = 2 * eigen.eigenvectors().row(0).transpose().array().square();
}

// The functions L_a(x) L_b(y) of degree a + b <= degree, by degree: hierarchical, the first one 1
std::vector<Function> Basis(int degree)
{
	std::vector<Function> basis;
	for (int d = 0; d <= degree; d++)
		for (int a = d; a >= 0; a--)
			basis.push_back({ a, d - a });
	return basis;
}

// The HHO matrix of the reference square [-1, 1]^2: the cell DoFs, then those of its faces, left, right, bottom, top
// (L_m(s), m <= k, s the coordinate along the face: the same functions from both cells of a face). The reconstruction R
// in P^{k+1}: (grad R, grad w) = (grad u_T, grad w) + sum_F (u_F - u_T, grad w.n)_F; the stabilization:
// sum_F h_F^-1 |pi_F(u_F - u_T - (R - pi_T R))|_F^2, h_F = 2. In 2D, the matrix of every square. nT, nF: the DoFs of a
// cell, of a face.
Matrix LocalMatrix(int k, int l, int& nT, int& nF)
{
	std::vector<Function> cell = Basis(l), reconstruction = Basis(k + 1);
	reconstruction.erase(reconstruction.begin()); // R without its constant: its gradient
	nT = (int)cell.size();
	nF = k + 1;
	int nR = (int)reconstruction.size(), n = nT + 4 * nF;
	Vector x, w;
	GaussLegendre(k + l + 3, x, w);
	auto value = [](Function f, double X, double Y) { return Legendre(f.first, X).first * Legendre(f.second, Y).first; };
	auto grad = [](Function f, double X, double Y)
	{
		auto [Lx, dLx] = Legendre(f.first, X);
		auto [Ly, dLy] = Legendre(f.second, Y);
		return Eigen::Vector2d(dLx * Ly, Lx * dLy);
	};
	const int faces[4][2] = { { 0, -1 }, { 0, 1 }, { 1, -1 }, { 1, 1 } }; // (the axis of the normal, its sign)
	auto onFace = [&](int f, int p) { return faces[f][0] == 0 ? Eigen::Vector2d(faces[f][1], x[p]) : Eigen::Vector2d(x[p], faces[f][1]); };

	Matrix K = Matrix::Zero(nR, nR), B = Matrix::Zero(nR, n);
	for (int p = 0; p < x.size(); p++)
		for (int q = 0; q < x.size(); q++)
			for (int i = 0; i < nR; i++)
			{
				Eigen::Vector2d gi = grad(reconstruction[i], x[p], x[q]);
				for (int j = 0; j < nR; j++)
					K(i, j) += w[p] * w[q] * gi.dot(grad(reconstruction[j], x[p], x[q]));
				for (int j = 0; j < nT; j++)
					B(i, j) += w[p] * w[q] * gi.dot(grad(cell[j], x[p], x[q]));
			}
	for (int f = 0; f < 4; f++)
		for (int p = 0; p < x.size(); p++)
		{
			Eigen::Vector2d point = onFace(f, p);
			for (int i = 0; i < nR; i++)
			{
				double gn = faces[f][1] * grad(reconstruction[i], point[0], point[1])[faces[f][0]];
				for (int j = 0; j < nT; j++)
					B(i, j) -= w[p] * value(cell[j], point[0], point[1]) * gn;
				for (int m = 0; m < nF; m++)
					B(i, nT + f * nF + m) += w[p] * Legendre(m, x[p]).first * gn;
			}
		}
	Matrix R = K.ldlt().solve(B); // the coefficients of the reconstruction
	Matrix A = B.transpose() * R;

	std::vector<int> high; // R - pi_T R: the functions of R not in P^l
	for (int i = 0; i < nR; i++)
		if (reconstruction[i].first + reconstruction[i].second > l)
			high.push_back(i);
	Vector mass(nF);
	for (int m = 0; m < nF; m++)
		mass[m] = 2.0 / (2 * m + 1);
	for (int f = 0; f < 4; f++)
	{
		auto project = [&](const std::vector<Function>& functions)
		{
			Matrix P = Matrix::Zero(nF, functions.size());
			for (int m = 0; m < nF; m++)
				for (size_t j = 0; j < functions.size(); j++)
					for (int p = 0; p < x.size(); p++)
						P(m, j) += w[p] * Legendre(m, x[p]).first * value(functions[j], onFace(f, p)[0], onFace(f, p)[1]) / mass[m];
			return P;
		};
		Matrix projectedR = project(reconstruction), D = Matrix::Zero(nF, n);
		for (int i : high)
			D -= projectedR.col(i) * R.row(i);
		D.leftCols(nT) -= project(cell);
		D.block(0, nT + f * nF, nF, nF) += Matrix::Identity(nF, nF);
		A += D.transpose() * mass.asDiagonal() * D / 2;
	}
	return A;
}

int main(int argc, char** argv)
{
	PetscInt n = 64, k = 1, l, iterations;
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
	Matrix A_local = LocalMatrix((int)k, (int)l, nT, nF);

	// The partition: the cells of each process, rows of the mesh (a horizontal strip); each face owned by the process of
	// its lower or left cell. Numbered in rank order: row by row, the cells of a row, and its faces: the n - 1 vertical
	// ones between its cells, then the n horizontal ones above it (none above the last row).
	auto facesBefore = [&](PetscInt j) { return j * (2 * n - 1); };                // the faces of the rows below j
	auto vertical = [&](PetscInt i, PetscInt j) { return facesBefore(j) + i; };         // between the cells (i, j) and (i + 1, j)
	auto horizontal = [&](PetscInt i, PetscInt j) { return facesBefore(j) + n - 1 + i; }; // between the cells (i, j) and (i, j + 1)
	std::vector<PetscInt> myRows;
	for (PetscInt j = 0; j < n; j++)
		if (j * size / n == rank)
			myRows.push_back(j);
	PetscInt nCells = n * (PetscInt)myRows.size(), nFaces = 0;
	for (PetscInt j : myRows)
		nFaces += j < n - 1 ? 2 * n - 1 : n - 1;

	// The four blocks, preallocated (non-zeros per row in the columns of the process, and in the others)
	PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, nT * nCells, nT * nCells, PETSC_DETERMINE, PETSC_DETERMINE, nT, nullptr, 0, nullptr, &A_TT));
	PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, nT * nCells, nF * nFaces, PETSC_DETERMINE, PETSC_DETERMINE, 4 * nF, nullptr, 4 * nF, nullptr, &A_TF));
	PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, nF * nFaces, nT * nCells, PETSC_DETERMINE, PETSC_DETERMINE, 2 * nT, nullptr, 2 * nT, nullptr, &A_FT));
	PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, nF * nFaces, nF * nFaces, PETSC_DETERMINE, PETSC_DETERMINE, 7 * nF, nullptr, 7 * nF, nullptr, &A_FF));

	// The local matrices of the cells of the process, added; the DoFs of the boundary faces (Dirichlet, homogeneous)
	// dropped
	for (PetscInt j : myRows)
	{
		for (PetscInt i = 0; i < n; i++)
		{
			PetscInt T = j * n + i;
			PetscInt faces[4] = { i > 0 ? vertical(i - 1, j) : -1, i < n - 1 ? vertical(i, j) : -1, j > 0 ? horizontal(i, j - 1) : -1, j < n - 1 ? horizontal(i, j) : -1 };
			double kappa = (8 * i / n + 8 * j / n) % 2 ? 10 : 1; // the diffusion coefficient of the cell
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
			Matrix a = kappa * A_local;
			RowMajorMatrix TT = a(cellLocal, cellLocal), TF = a(cellLocal, faceLocal), FT = a(faceLocal, cellLocal), FF = a(faceLocal, faceLocal);
			PetscInt cT = (PetscInt)cellDoFs.size(), cF = (PetscInt)faceDoFs.size();
			PetscCall(MatSetValues(A_TT, cT, cellDoFs.data(), cT, cellDoFs.data(), TT.data(), ADD_VALUES));
			PetscCall(MatSetValues(A_TF, cT, cellDoFs.data(), cF, faceDoFs.data(), TF.data(), ADD_VALUES));
			PetscCall(MatSetValues(A_FT, cF, faceDoFs.data(), cT, cellDoFs.data(), FT.data(), ADD_VALUES));
			PetscCall(MatSetValues(A_FF, cF, faceDoFs.data(), cF, faceDoFs.data(), FF.data(), ADD_VALUES));
		}
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
		PetscInt nDoFs, blockSize = part == 0 ? nT : nF;
		PetscCall(ISGetLocalSize(rows[part], &nDoFs));
		PetscCall(ISGetIndices(rows[part], &indices));
		for (PetscInt p = 0; p < nDoFs; p += blockSize)
		{
			if (part == 0)
				PetscCall(VecSetValue(b, indices[p], 1.0 / (n * n), INSERT_VALUES));
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
