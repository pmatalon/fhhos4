// Unit tests of the local Galerkin products of the coarsening passes of U-AMG (src/Solver/Multigrid/UncondensedAMG/
// LocalMatrices.h): the sum of the local matrices of the aggregates, computed from those of the cells, must be the global
// product P^T S P, for a symmetric operator (lower triangular parts) and a non-symmetric one (full).
#include <gtest/gtest.h>
#include <vector>
#include "Solver/Multigrid/UncondensedAMG/LocalMatrices.h"

namespace
{
	// Hybrid system of degree 0 on the n x n Cartesian mesh of the unit square (interior faces), as in
	// UAMGLibraryTest.cpp: cell coefficients a_T, cell-face couplings -a_T (1 + asymmetry s) in A_TF and
	// -a_T (1 - asymmetry s) in A_FT (s = 1 on the right and top faces of the cell)
	struct System
	{
		SparseMatrix A_TT, A_TF, A_FT, A_FF, S; // S: the condensed operator A_FF - A_FT A_TT^-1 A_TF

		System(int n, double asymmetry)
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

	// Face prolongation with interface collapsing (as UncondensedLevel::BuildQ_F()): a kept face takes the value of its
	// coarse face, a removed one the average of those of its aggregate
	SparseMatrix FaceProlongation(const HybridAlgebraicMesh& mesh)
	{
		std::vector<Eigen::Triplet<double, SparseMatrixIndex>> P;
		for (const HybridAlgebraicFace& f : mesh.Faces)
		{
			if (f.IsRemovedOnCoarseMesh)
			{
				const HybridElementAggregate* K = f.Elements[0]->CoarseElement;
				for (const HybridFaceAggregate* cf : K->CoarseFaces)
					P.push_back({ (SparseMatrixIndex)f.Number, (SparseMatrixIndex)cf->Number, 1.0 / K->CoarseFaces.size() });
			}
			else
				P.push_back({ (SparseMatrixIndex)f.Number, (SparseMatrixIndex)f.CoarseFace->Number, 1.0 });
		}
		SparseMatrix Q(mesh.Faces.size(), mesh.CoarseFaces.size());
		Q.setFromTriplets(P.begin(), P.end());
		return Q;
	}

	void CheckGalerkinProduct(double asymmetry)
	{
		System system(16, asymmetry);
		bool symmetric = asymmetry == 0;
		HybridAlgebraicMesh mesh(&system.A_TT, &system.A_TF, nullptr, 1, 1, 0.25);
		if (!symmetric)
		{
			mesh.A_F_T = &system.A_FT;
			mesh.Symmetric = false;
		}
		mesh.Build();
		bool coarsest = false;
		mesh.Coarsen(H_CoarsStgy::MultiplePairwiseAggregation, FaceCoarseningStrategy::InterfaceCollapsing, coarsest);
		ASSERT_FALSE(coarsest);
		SparseMatrix P = FaceProlongation(mesh);

		LocalMatrices cells, aggregates;
		cells.Decompose(system.S, mesh, 1, symmetric);
		EXPECT_LT((cells.Assemble(mesh.Faces.size()) - system.S).norm(), 1e-12 * system.S.norm()) << "asymmetry " << asymmetry;
		aggregates.GalerkinProduct(cells, mesh, P);
		SparseMatrix local = aggregates.Assemble(mesh.CoarseFaces.size());
		SparseMatrix global = SparseMatrixOps::Multiply(SparseMatrixOps::Multiply(SparseMatrixOps::Transpose(P), system.S), P);
		EXPECT_LT((local - global).norm(), 1e-12 * global.norm()) << "asymmetry " << asymmetry;
		if (!symmetric)
			EXPECT_GT((global - SparseMatrix(global.transpose())).norm(), 0.1 * global.norm()); // the test is meaningful
	}
}

TEST(LocalMatricesTest, GalerkinProductOfASymmetricOperator)
{
	CheckGalerkinProduct(0);
}

TEST(LocalMatricesTest, GalerkinProductOfANonsymmetricOperator)
{
	CheckGalerkinProduct(0.5);
}
