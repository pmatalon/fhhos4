// The HHO discretization of the C++ examples (petsc_preconditioner.cpp, fhhos4_solver.cpp): -div(kappa grad u) = 1 on
// the unit square, n x n squares, kappa 1 or 10 on a checkerboard of 8 x 8 squares, homogeneous Dirichlet conditions,
// with polynomials of degree k on the faces and l in the cells (Legendre bases, ordered by degree); the mesh distributed
// by rows among the processes.
#pragma once
#include <array>
#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>
#include <Eigen/Dense>

namespace HHO
{
	using Matrix = Eigen::MatrixXd;
	using Vector = Eigen::VectorXd;
	using Function = std::pair<int, int>; // L_a(x) L_b(y)

	// The Legendre polynomial of degree m at x, and its derivative
	inline std::pair<double, double> Legendre(int m, double x)
	{
		double L = 1;
		double previous = 0;
		double dL = 0;
		double dPrevious = 0;
		for (int i = 0; i < m; i++)
		{
			double next = ((2 * i + 1) * x * L - i * previous) / (i + 1);
			double dNext = dPrevious + (2 * i + 1) * L;
			previous = L;
			L = next;
			dPrevious = dL;
			dL = dNext;
		}
		return { L, dL };
	}

	// The Gauss-Legendre quadrature of q points on [-1, 1] (Golub-Welsch)
	inline void GaussLegendre(int q, Vector& x, Vector& w)
	{
		Matrix J = Matrix::Zero(q, q);
		for (int i = 1; i < q; i++)
			J(i, i - 1) = J(i - 1, i) = i / std::sqrt(4.0 * i * i - 1);
		Eigen::SelfAdjointEigenSolver<Matrix> eigen(J);
		x = eigen.eigenvalues();
		w = 2 * eigen.eigenvectors().row(0).transpose().array().square();
	}

	// The functions L_a(x) L_b(y) of degree a + b <= degree, by degree: hierarchical, the first one 1
	inline std::vector<Function> Basis(int degree)
	{
		std::vector<Function> basis;
		for (int d = 0; d <= degree; d++)
			for (int a = d; a >= 0; a--)
				basis.push_back({ a, d - a });
		return basis;
	}

	// The HHO matrix of the reference square [-1, 1]^2: the cell DoFs, then those of its faces, left, right, bottom, top
	// (L_m(s), m <= k, s the coordinate along the face: the same functions from both cells of a face). The reconstruction
	// R in P^{k+1}: (grad R, grad w) = (grad u_T, grad w) + sum_F (u_F - u_T, grad w.n)_F; the stabilization:
	// sum_F h_F^-1 |pi_F(u_F - u_T - (R - pi_T R))|_F^2, h_F = 2. In 2D, the matrix of every square. nT, nF: the DoFs of a
	// cell, of a face. The interpolation of 1: 1 on the first function of each cell and face, 0 on the others.
	inline Matrix LocalMatrixReferenceSquare(int k, int l, int& nT, int& nF)
	{
		std::vector<Function> cell = Basis(l);
		std::vector<Function> reconstruction = Basis(k + 1);
		reconstruction.erase(reconstruction.begin()); // R without its constant: its gradient
		nT = (int)cell.size();
		nF = k + 1;
		int nR = (int)reconstruction.size();
		int n = nT + 4 * nF;
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

		Matrix K = Matrix::Zero(nR, nR);
		Matrix B = Matrix::Zero(nR, n);
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
			Matrix projectedR = project(reconstruction);
			Matrix D = Matrix::Zero(nF, n);
			for (int i : high)
				D -= projectedR.col(i) * R.row(i);
			D.leftCols(nT) -= project(cell);
			D.block(0, nT + f * nF, nF, nF) += Matrix::Identity(nF, nF);
			A += D.transpose() * mass.asDiagonal() * D / 2;
		}
		return A;
	}

	// The mesh, n x n squares, distributed among the processes, seen by the assemblies as an unstructured mesh: its cells
	// and faces by their global numbers, in rank order (those of the process of rank 0 first, etc.), each cell with its
	// faces, each face with its cells. Here, the cells of each process are rows of the mesh (a horizontal strip), and each
	// face is owned by the process of its lower or left cell; numbered row by row: the cells of a row (the cell (i, j) is
	// j n + i), and its faces: the n - 1 vertical ones between its cells, then the n horizontal ones above it (none above
	// the last row).
	struct Mesh
	{
		std::int64_t n;
		int Processes;

		// The faces of the cell T, in the order of the local matrix: left, right, bottom, top; -1 on the boundary
		// (Dirichlet: no unknowns)
		std::array<std::int64_t, 4> Faces(std::int64_t T) const
		{
			std::int64_t i = T % n;
			std::int64_t j = T / n;
			return { i > 0 ? FacesBefore(j) + i - 1 : -1, i < n - 1 ? FacesBefore(j) + i : -1,
			         j > 0 ? FacesBefore(j - 1) + n - 1 + i : -1, j < n - 1 ? FacesBefore(j) + n - 1 + i : -1 };
		}

		// The two cells of the face F (not on the boundary)
		std::array<std::int64_t, 2> Cells(std::int64_t F) const
		{
			std::int64_t j = F / (2 * n - 1);
			std::int64_t i = F % (2 * n - 1);
			return i < n - 1 ? std::array<std::int64_t, 2>{ j * n + i, j * n + i + 1 }      // vertical
			                 : std::array<std::int64_t, 2>{ j * n + i - (n - 1), (j + 1) * n + i - (n - 1) }; // horizontal
		}

		// The diffusion coefficient of the cell T
		double Kappa(std::int64_t T) const
		{
			return (8 * (T % n) / n + 8 * (T / n) / n) % 2 ? 10 : 1;
		}

		// The measure of the cell T (its area)
		double Measure(std::int64_t) const
		{
			return 1.0 / ((double)n * n);
		}

		// The cells of the process of rank r, [first, last)
		std::pair<std::int64_t, std::int64_t> CellOwnershipRange(int r) const
		{
			auto [first, last] = Rows(r);
			return { first * n, last * n };
		}

		// The faces of the process of rank r, [first, last)
		std::pair<std::int64_t, std::int64_t> FaceOwnershipRange(int r) const
		{
			auto [first, last] = Rows(r);
			return { FacesBefore(first), FacesBefore(last) };
		}

		// The rows of the process of rank r, [first, last): those j with j Processes / n = r
		std::pair<std::int64_t, std::int64_t> Rows(int r) const
		{
			return { (r * n + Processes - 1) / Processes, ((r + 1) * n + Processes - 1) / Processes };
		}

		// The faces of the rows below j (2 n - 1 each, n - 1 for the last row)
		std::int64_t FacesBefore(std::int64_t j) const
		{
			return j * (2 * n - 1) - (j == n ? n : 0);
		}
	};
}
