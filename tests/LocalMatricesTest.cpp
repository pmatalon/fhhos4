// Unit tests of the local Galerkin products of the coarsening passes of U-AMG (src/Solver/Multigrid/UncondensedAMG/
// LocalMatrices.h): the sum of the local matrices of the aggregates, computed from those of the cells, must be the global
// product P^T S P, for a symmetric operator (lower triangular parts) and a non-symmetric one (full).
#include <gtest/gtest.h>
#include <vector>
#include "Solver/Multigrid/UncondensedAMG/LocalMatrices.h"
#include "support/CartesianHybridSystem.h"

namespace
{
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
		fhhos4_tests::CartesianHybridSystem system(16, asymmetry);
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
