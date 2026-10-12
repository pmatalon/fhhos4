#pragma once
#include <vector>
#include "Utils/Types.h"

namespace fhhos4_tests
{
	// Hybrid system of degree 0 on the n x n Cartesian mesh of the unit square (interior faces), as in
	// UAMGLibraryTest.cpp: cell coefficients a_T, cell-face couplings -a_T (1 + asymmetry s) in A_TF and
	// -a_T (1 - asymmetry s) in A_FT (s = 1 on the right and top faces of the cell)
	struct CartesianHybridSystem
	{
		SparseMatrix A_TT, A_TF, A_FT, A_FF, S; // S: the condensed operator A_FF - A_FT A_TT^-1 A_TF

		CartesianHybridSystem(int n, double asymmetry)
		{
			int nVertical = (n - 1) * n;
			auto vertical = [n](int i, int j) { return j * (n - 1) + i; };
			auto horizontal = [n, nVertical](int i, int j) { return nVertical + j * n + i; };
			std::vector<Eigen::Triplet<double, SparseMatrixIndex>> TT, TF, FT, FF;
			for (int j = 0; j < n; j++)
			{
				for (int i = 0; i < n; i++)
				{
					int T = j * n + i;
					double a = 1 + (T * 7 % 5);
					TT.push_back({ T, T, 4 * a });
					std::vector<std::pair<int, int>> faces;
					if (i > 0)     faces.push_back({ vertical(i - 1, j), -1 });
					if (i < n - 1) faces.push_back({ vertical(i, j), 1 });
					if (j > 0)     faces.push_back({ horizontal(i, j - 1), -1 });
					if (j < n - 1) faces.push_back({ horizontal(i, j), 1 });
					for (auto [F, s] : faces)
					{
						TF.push_back({ T, F, -a * (1 + asymmetry * s) });
						FT.push_back({ F, T, -a * (1 - asymmetry * s) });
						FF.push_back({ F, F, a });
					}
				}
			}
			int nCells = n * n, nFaces = 2 * nVertical;
			A_TT.resize(nCells, nCells);
			A_TT.setFromTriplets(TT.begin(), TT.end());
			A_TF.resize(nCells, nFaces);
			A_TF.setFromTriplets(TF.begin(), TF.end());
			A_FT.resize(nFaces, nCells);
			A_FT.setFromTriplets(FT.begin(), FT.end());
			A_FF.resize(nFaces, nFaces);
			A_FF.setFromTriplets(FF.begin(), FF.end());
			Vector invDiag = A_TT.diagonal().cwiseInverse();
			S = A_FF - SparseMatrix(A_FT * invDiag.asDiagonal() * A_TF);
			S.makeCompressed();
		}
	};
}
