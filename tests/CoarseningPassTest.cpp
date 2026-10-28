// Tests of the coarsening passes of U-AMG on flat arrays (src/Solver/Multigrid/UncondensedAMG/CoarseningPass.h), which
// run with the default options: they must build exactly the hierarchy of the reference implementation
// (HybridAlgebraicMesh, with UncondensedLevel::UseFlatPasses = false).
// - CoarseningPassProgram: the program on HHO systems (simplicial, Cartesian and polygonal meshes, heterogeneous and
//   anisotropic problems, mixed Dirichlet-Neumann conditions (faces of a single cell in the system), k = 0 to 2, the
//   h-coarsening of the degree-k blocks, orthonormal bases, the non-symmetric path on a symmetric matrix, the
//   coarsening prolongations 3 to 6, the face prolongations 1 and 2, the multigrid
//   prolongations other than the chained ones, the aggregation with all the neighbours, the inverses of the diagonal
//   blocks, p-levels after h-levels): same iteration count and L2 error, bit for bit. On 1 thread: on several, the L2 error
//   (an OpenMP reduction) and the polygonal meshes (their agglomeration) vary from one run to the next.
// - CoarseningPass: the operators and prolongations of all the levels, on a symmetric and a non-symmetric system built
//   here, on all the threads; also with the coarse operators Q_F^T A Q_F (not Galerkin).
// The flat passes then keep the product P of their prolongations, as the reference (UncondensedLevel::
// KeepPassProlongations = false). PassProlongations: the transfers by the prolongations of the passes (the default)
// against P: the same operators, the same transfers up to rounding, the same iterations.
#include <gtest/gtest.h>
#include <cstring>
#include <string>
#include "Solver/Multigrid/UncondensedAMG/UncondensedAMG.h"
#include "support/CartesianHybridSystem.h"
#include "support/RunHelper.h"

using namespace fhhos4_tests;

namespace
{
	// Runs the coarsening passes of HybridAlgebraicMesh (the reference) during its lifetime
	struct ReferencePasses
	{
		ReferencePasses() { UncondensedLevel::UseFlatPasses = false; }
		~ReferencePasses() { UncondensedLevel::UseFlatPasses = true; }
	};

	// The levels keep the product P of the prolongations of their passes during its lifetime, as the reference
	struct ChainedProlongations
	{
		ChainedProlongations() { UncondensedLevel::KeepPassProlongations = false; }
		~ChainedProlongations() { UncondensedLevel::KeepPassProlongations = true; }
	};

	// The levels keep the prolongations of their passes during its lifetime, also at degree 0 (by default, from blocks of
	// 2 DoFs only)
	struct PassProlongationsAtDegree0
	{
		int MinBlockSize = UncondensedLevel::PassProlongationsMinFaceBlockSize;
		PassProlongationsAtDegree0() { UncondensedLevel::PassProlongationsMinFaceBlockSize = 1; }
		~PassProlongationsAtDegree0() { UncondensedLevel::PassProlongationsMinFaceBlockSize = MinBlockSize; }
	};
}

//------------------------------------------------------------------------------------------------------------------//
//                                                Through the program                                                 //
//------------------------------------------------------------------------------------------------------------------//

struct FlatPassCase
{
	std::string Name;
	std::string Geo;
	std::string Mesh;
	std::string Mesher;                // "" for the default
	int N;
	int K;
	bool HCoarseningOfBlocks = false;  // -hp-cs h: the coarsening passes run on the degree-k blocks
	bool Orthonormal = false;          // -e-ogb 3 -f-ogb 3: the constant does not have coordinate 1
	bool Nonsymmetric = false;         // -s gcruamg -uamg-nonsym
	double Heterogeneity = 1;          // -heterog
	double Anisotropy = 1;             // -aniso
	bool PolyMeshSharedFaces = false;  // -mesh poly -polymesh-fcs n: neighbouring cells share several faces
	int CoarseningProlong = 0;         // -coarsening-prolong (0: default)
	int FaceProlong = 0;               // -face-prolong (0: default)
	bool DiagInverse = false;          // -diag-block-solve inv: the inverses of the diagonal blocks (block Jacobi of P_F)
	int MeshesThenP = 0;               // -hp-cs h_p -num-meshes N: N meshes, then p-levels (the coarse A_F_F is computed)
	int Prolong = 0;                   // -prolong: the multigrid prolongation (0: default, the chained coarsening prolongations)
	bool AllNeighbours = false;        // -cs mn: each cell aggregated with all its neighbours
	std::string BC;                    // -bc ("": the default, Dirichlet; m: mixed Dirichlet-Neumann)
	bool CollapseBoundaryFaces = false; // -uamg-bfc: the Neumann faces collapsed into one coarse face per aggregate
};

inline std::ostream& operator<<(std::ostream& os, const FlatPassCase& c)
{
	return os << c.Name;
}

class CoarseningPassProgramTest : public ::testing::TestWithParam<FlatPassCase>
{
};

// E.g. Cube_cart_n8_k1_h: ./bin/fhhos4 -geo cube -mesh cart -mesher inhouse -n 8 -k 1 -hp-cs h -s fcguamg -smoothers bgs -no-cache -threads 1
TEST_P(CoarseningPassProgramTest, SameAsReference)
{
	FlatPassCase c = GetParam();
	SequentialExecution sequential;
	auto run = [&]()
	{
		ProgramArguments args;
		args.Problem.GeoCode = c.Geo;
		args.Problem.HeterogeneityRatio = c.Heterogeneity;
		args.Problem.AnisotropyRatio = c.Anisotropy;
		if (!c.BC.empty())
			args.Problem.BCCode = c.BC;
		args.Solver.MG.UAMGCollapseBoundaryFaces = c.CollapseBoundaryFaces;
		args.Discretization.MeshCode = c.Mesh;
		if (!c.Mesher.empty())
			args.Discretization.Mesher = c.Mesher;
		if (c.PolyMeshSharedFaces)
			args.Discretization.PolyMeshFaceCoarseningStgy = FaceCoarseningStrategy::None;
		args.Discretization.N = c.N;
		args.Discretization.PolyDegree = c.K + 1;
		if (c.Orthonormal)
		{
			args.Discretization.OrthogonalizeElemBasesCode = 3;
			args.Discretization.OrthogonalizeFaceBasesCode = 3;
		}
		args.Solver.SolverCode = c.Nonsymmetric ? "gcruamg" : "fcguamg";
		args.Solver.MG.PreSmoother = SmootherType::BlockGaussSeidel; // -smoothers bgs
		args.Solver.MG.UAMGNonsymmetric = c.Nonsymmetric;
		args.Solver.MG.CoarseningProlongationCode = c.CoarseningProlong;
		args.Solver.MG.FaceProlongationCode = c.FaceProlong;
		args.Solver.MG.ProlongationCode = c.Prolong;
		if (c.AllNeighbours)
			args.Solver.MG.H_CS = H_CoarsStgy::MultipleAgglomerationCoarseningByFaceNeighbours;
		if (c.HCoarseningOfBlocks)
			args.Solver.MG.HP_CS = HP_CoarsStgy::H_only;
		if (c.MeshesThenP > 0)
		{
			args.Solver.MG.HP_CS = HP_CoarsStgy::H_then_P;
			args.Solver.MG.NumberOfMeshes = c.MeshesThenP;
		}
		args.Actions.UseCache = false;
		BlockDiagonalSolver::DefaultMethod = c.DiagInverse ? DiagBlockSolveMethod::Inverse : DiagBlockSolveMethod::LU;
		ProgramResults results = RunDiffusionHHO(args, true, !c.HCoarseningOfBlocks && c.MeshesThenP == 0);
		BlockDiagonalSolver::DefaultMethod = DiagBlockSolveMethod::LU;
		return results;
	};

	ChainedProlongations chained;
	ProgramResults reference;
	{
		ReferencePasses referencePasses;
		reference = run();
	}
	ProgramResults flat = run();
	EXPECT_GT(reference.IterationCount, 0);
	EXPECT_LT(reference.IterationCount, 200); // converged (default maximum: 200)
	EXPECT_EQ(flat.IterationCount, reference.IterationCount);
	EXPECT_EQ(flat.L2Error, reference.L2Error); // the same solution
}

std::vector<FlatPassCase> FlatPassCases()
{
	std::vector<FlatPassCase> cases;
#ifdef ENABLE_2D
	cases.push_back({ "Square_tri_n16_k0", "square", "tri", "", 16, 0 });
	cases.push_back({ "Square_tri_n16_k1", "square", "tri", "", 16, 1 });
	cases.push_back({ "Square_tri_n8_k2_h", "square", "tri", "", 8, 2, true });
	for (int coarseningProlong : { 3, 4, 6 })
	{
		FlatPassCase c = { "Square_tri_n16_k0_fp2_cp" + std::to_string(coarseningProlong), "square", "tri", "", 16, 0 };
		c.FaceProlong = 2;
		c.CoarseningProlong = coarseningProlong;
		cases.push_back(c);
	}
	{
		FlatPassCase c = { "Square_poly_n32_k0_sharedfaces", "square", "poly", "", 32, 0 };
		c.PolyMeshSharedFaces = true;
		cases.push_back(c);
	}
	{
		FlatPassCase c = { "Square4quadrants_cart_n32_k0_heterog1e4", "square4quadrants", "cart", "inhouse", 32, 0 };
		c.Heterogeneity = 1e4;
		cases.push_back(c);
	}
	{
		FlatPassCase c = { "Square_cart_n32_k0_aniso100", "square", "cart", "inhouse", 32, 0 };
		c.Anisotropy = 100;
		cases.push_back(c);
	}
	// Mixed Dirichlet-Neumann conditions: the Neumann faces, of a single cell, removed by the first pass, or collapsed
	// into one coarse face per aggregate (-uamg-bfc)
	for (bool collapse : { false, true })
	{
		for (int k : { 0, 1 })
		{
			FlatPassCase c = { "Square_tri_n16_k" + std::to_string(k) + "_mixed" + (collapse ? "_bfc" : ""), "square", "tri", "", 16, k };
			c.BC = "m";
			c.CollapseBoundaryFaces = collapse;
			cases.push_back(c);
		}
		FlatPassCase c = { "Square_cart_n32_k0_mixed" + std::string(collapse ? "_bfc" : ""), "square", "cart", "inhouse", 32, 0 };
		c.BC = "m";
		c.CollapseBoundaryFaces = collapse;
		cases.push_back(c);
	}
	{
		FlatPassCase c = { "Square_tri_n16_k0_mixed_cp6", "square", "tri", "", 16, 0 };
		c.BC = "m";
		c.CoarseningProlong = 6;
		cases.push_back(c);
	}
#endif
#ifdef ENABLE_3D
	cases.push_back({ "Cube_tetra_n8_k0", "cube", "tetra", "", 8, 0 });
	{
		FlatPassCase c = { "Cube_tetra_n8_k1_orthonormal", "cube", "tetra", "", 8, 1 };
		c.Orthonormal = true;
		cases.push_back(c);
	}
	cases.push_back({ "Cube_cart_n8_k1_h", "cube", "cart", "inhouse", 8, 1, true });
	{
		FlatPassCase c = { "Cube_tetra_n8_k0_nonsym", "cube", "tetra", "", 8, 0 };
		c.Nonsymmetric = true;
		cases.push_back(c);
	}
	for (int coarseningProlong : { 3, 4, 5 })
	{
		FlatPassCase c = { "Cube_tetra_n8_k0_cp" + std::to_string(coarseningProlong), "cube", "tetra", "", 8, 0 };
		c.CoarseningProlong = coarseningProlong;
		cases.push_back(c);
	}
	{
		FlatPassCase c = { "Cube_cart_n8_k1_h_cp4_diaginverse", "cube", "cart", "inhouse", 8, 1, true };
		c.CoarseningProlong = 4;
		c.DiagInverse = true;
		cases.push_back(c);
	}
	{
		FlatPassCase c = { "Cube_cart_n8_k1_2meshes_then_p", "cube", "cart", "inhouse", 8, 1 };
		c.MeshesThenP = 2;
		cases.push_back(c);
	}
	{
		FlatPassCase c = { "Cube_tetra_n8_k0_mn", "cube", "tetra", "", 8, 0 };
		c.AllNeighbours = true;
		cases.push_back(c);
	}
	for (int prolong : { 3, 5, 6, 8 })
	{
		FlatPassCase c = { "Cube_tetra_n8_k0_prolong" + std::to_string(prolong), "cube", "tetra", "", 8, 0 };
		c.Prolong = prolong;
		cases.push_back(c);
	}
	{
		FlatPassCase c = { "Cube_cart_n8_k1_h_prolong5", "cube", "cart", "inhouse", 8, 1, true };
		c.Prolong = 5;
		cases.push_back(c);
	}
#endif
	return cases;
}

INSTANTIATE_TEST_SUITE_P(CoarseningPassProgram, CoarseningPassProgramTest, ::testing::ValuesIn(FlatPassCases()),
	[](const ::testing::TestParamInfo<FlatPassCase>& info) { return info.param.Name; });

//------------------------------------------------------------------------------------------------------------------//
//                                         The levels, on a non-symmetric system                                      //
//------------------------------------------------------------------------------------------------------------------//

namespace
{
	class HierarchyAccess : public UncondensedAMG
	{
	public:
		using UncondensedAMG::UncondensedAMG;
		Level* Fine() { return this->_fineLevel; }
	};

	class LevelAccess : public Level
	{
	public:
		static const SparseMatrix& GetP(Level* l) { return static_cast<LevelAccess*>(l)->P; }
	};

	bool BitIdentical(const SparseMatrix& X, const SparseMatrix& Y)
	{
		return X.rows() == Y.rows() && X.cols() == Y.cols() && X.nonZeros() == Y.nonZeros()
			&& std::equal(X.outerIndexPtr(), X.outerIndexPtr() + X.rows() + 1, Y.outerIndexPtr())
			&& std::equal(X.innerIndexPtr(), X.innerIndexPtr() + X.nonZeros(), Y.innerIndexPtr())
			&& std::memcmp(X.valuePtr(), Y.valuePtr(), X.nonZeros() * sizeof(double)) == 0;
	}

	// U-AMG as the library builds it, for the system (2D, degree 0), with the interpolations of 1 given (c = 1 or not)
	void Setup(HierarchyAccess& amg, const CartesianHybridSystem& system, bool symmetric, double constantCoordinate, bool galerkin)
	{
		amg.UseGalerkinOperator = galerkin;
		amg.MatrixMaxSizeForCoarsestLevel = 30;
		amg.CoarseningFactor = 3.8;
		amg.Cycle = CycleType::K;
		amg.KCycleKrylov = symmetric ? KrylovMethod::FCG : KrylovMethod::GCR;
		amg.PreSmoother = SmootherType::BlockGaussSeidel;
		amg.PostSmoother = SmootherType::BackwardBlockGaussSeidel;
		Vector cellInterpOfOne = Vector::Constant(system.A_TF.rows(), constantCoordinate);
		Vector faceInterpOfOne = Vector::Constant(system.A_TF.cols(), constantCoordinate);
		for (BigNumber i = 0; i < faceInterpOfOne.rows(); i += 3)
			faceInterpOfOne[i] *= 2;
		amg.Setup(system.S, system.A_TT, system.A_TF, symmetric ? SparseMatrix() : system.A_FT, SparseMatrix(), cellInterpOfOne, faceInterpOfOne);
	}

	void CheckLevels(double asymmetry, double constantCoordinate, bool galerkin = true)
	{
		CartesianHybridSystem system(32, asymmetry);
		bool symmetric = asymmetry == 0;
		auto create = []()
		{
			return HierarchyAccess(2, 0, false, 1, 1, 0.25, UAMGFaceProlongation::BoundaryAggregatesInteriorAverage,
				UAMGProlongation::ReconstructSmoothedTraceOrInject, UAMGProlongation::ChainedCoarseningProlongations);
		};
		ChainedProlongations chained;
		HierarchyAccess reference = create();
		{
			ReferencePasses referencePasses;
			Setup(reference, system, symmetric, constantCoordinate, galerkin);
		}
		HierarchyAccess flat = create();
		Setup(flat, system, symmetric, constantCoordinate, galerkin);

		ASSERT_EQ(flat.NumberOfLevels(), reference.NumberOfLevels());
		ASSERT_GE(reference.NumberOfLevels(), 3); // several coarsening passes and levels
		Level* r = reference.Fine();
		Level* f = flat.Fine();
		for (int n = 0; r; n++, r = r->CoarserLevel, f = f->CoarserLevel)
		{
			EXPECT_TRUE(BitIdentical(*f->OperatorMatrix, *r->OperatorMatrix)) << "operator of level " << n;
			if (r->CoarserLevel)
				EXPECT_TRUE(BitIdentical(LevelAccess::GetP(f), LevelAccess::GetP(r))) << "prolongation of level " << n;
		}
	}
}

TEST(CoarseningPass, SymmetricSystem)
{
	CheckLevels(0, 1);
}

TEST(CoarseningPass, SymmetricSystemRescaledBases)
{
	CheckLevels(0, 0.5);
}

TEST(CoarseningPass, NonsymmetricSystem)
{
	CheckLevels(0.5, 1);
}

// The coarse operators Q_F^T A Q_F: the face prolongation chained over the coarsening passes
TEST(CoarseningPass, NotGalerkin)
{
	CheckLevels(0, 1, false);
}

//------------------------------------------------------------------------------------------------------------------//
//                     The transfers by the prolongations of the passes (KeepPassProlongations)                      //
//------------------------------------------------------------------------------------------------------------------//

namespace
{
	// The prolongations of the passes against their product P: the same levels and operators (the Galerkin products use
	// the passes' prolongations either way), the same transfers up to rounding. At degree 0: the passes' prolongations
	// forced.
	void CheckPassProlongations(double asymmetry, double constantCoordinate)
	{
		PassProlongationsAtDegree0 atDegree0;
		CartesianHybridSystem system(32, asymmetry);
		bool symmetric = asymmetry == 0;
		auto create = []()
		{
			return HierarchyAccess(2, 0, false, 1, 1, 0.25, UAMGFaceProlongation::BoundaryAggregatesInteriorAverage,
				UAMGProlongation::HarmonicExtensionOrInject, UAMGProlongation::ChainedCoarseningProlongations);
		};
		HierarchyAccess chained = create();
		{
			ChainedProlongations chainedProlongations;
			Setup(chained, system, symmetric, constantCoordinate, true);
		}
		HierarchyAccess passes = create();
		Setup(passes, system, symmetric, constantCoordinate, true);

		ASSERT_EQ(passes.NumberOfLevels(), chained.NumberOfLevels());
		ASSERT_GE(chained.NumberOfLevels(), 3);
		Level* c = chained.Fine();
		Level* p = passes.Fine();
		for (int n = 0; c; n++, c = c->CoarserLevel, p = p->CoarserLevel)
		{
			EXPECT_TRUE(BitIdentical(*p->OperatorMatrix, *c->OperatorMatrix)) << "operator of level " << n;
			if (!c->CoarserLevel)
				continue;
			EXPECT_EQ(LevelAccess::GetP(p).nonZeros(), 0) << "level " << n << ": P not computed";
			EXPECT_GT(p->ProlongCost(), 0);
			Vector e = Vector::LinSpaced(c->CoarserLevel->OperatorMatrix->rows(), -1, 2).array().sin();
			Vector r = Vector::LinSpaced(c->OperatorMatrix->rows(), 0, 3).array().cos();
			Vector Pe = c->Prolong(e), passesPe = p->Prolong(e);
			Vector Rr = c->Restrict(r), passesRr = p->Restrict(r);
			EXPECT_LE((passesPe - Pe).norm(), 1e-13 * Pe.norm()) << "prolongation of level " << n;
			EXPECT_LE((passesRr - Rr).norm(), 1e-13 * Rr.norm()) << "restriction of level " << n;
		}
	}
}

TEST(CoarseningPass, PassProlongations)
{
	CheckPassProlongations(0, 1);
}

TEST(CoarseningPass, PassProlongationsRescaledBases)
{
	CheckPassProlongations(0, 0.5);
}

TEST(CoarseningPass, PassProlongationsNonsymmetric)
{
	CheckPassProlongations(0.5, 1);
}

// The program: the same iterations with the transfers by the prolongations of the passes as with P (also at degree 0,
// forced). On 1 thread.
//   ./bin/fhhos4 -geo square -mesh tri -n 16 -k 1 -s fcguamg -smoothers bgs -no-cache -threads 1 (etc.)
class PassProlongationsProgramTest : public ::testing::TestWithParam<FlatPassCase>
{
};

TEST_P(PassProlongationsProgramTest, SameAsChained)
{
	FlatPassCase c = GetParam();
	SequentialExecution sequential;
	PassProlongationsAtDegree0 atDegree0;
	auto run = [&]()
	{
		ProgramArguments args;
		args.Problem.GeoCode = c.Geo;
		args.Discretization.MeshCode = c.Mesh;
		if (!c.Mesher.empty())
			args.Discretization.Mesher = c.Mesher;
		args.Discretization.N = c.N;
		args.Discretization.PolyDegree = c.K + 1;
		if (c.Orthonormal)
		{
			args.Discretization.OrthogonalizeElemBasesCode = 3;
			args.Discretization.OrthogonalizeFaceBasesCode = 3;
		}
		args.Solver.SolverCode = c.Nonsymmetric ? "gcruamg" : "fcguamg";
		args.Solver.MG.PreSmoother = SmootherType::BlockGaussSeidel; // -smoothers bgs
		args.Solver.MG.UAMGNonsymmetric = c.Nonsymmetric;
		args.Solver.MG.CoarseningProlongationCode = c.CoarseningProlong;
		if (c.HCoarseningOfBlocks)
			args.Solver.MG.HP_CS = HP_CoarsStgy::H_only;
		args.Actions.UseCache = false;
		return RunDiffusionHHO(args, true, !c.HCoarseningOfBlocks);
	};
	ProgramResults chained;
	{
		ChainedProlongations chainedProlongations;
		chained = run();
	}
	ProgramResults passes = run();
	EXPECT_GT(chained.IterationCount, 0);
	EXPECT_LT(chained.IterationCount, 200);
	EXPECT_EQ(passes.IterationCount, chained.IterationCount);
	EXPECT_NEAR(passes.L2Error, chained.L2Error, 1e-8 * chained.L2Error);
}

std::vector<FlatPassCase> PassProlongationsCases()
{
	std::vector<FlatPassCase> cases;
#ifdef ENABLE_2D
	cases.push_back({ "Square_tri_n16_k1", "square", "tri", "", 16, 1 });
#endif
#ifdef ENABLE_3D
	{
		FlatPassCase c = { "Cube_tetra_n8_k0_orthonormal", "cube", "tetra", "", 8, 0 };
		c.Orthonormal = true;
		cases.push_back(c);
	}
	cases.push_back({ "Cube_cart_n8_k1_h", "cube", "cart", "inhouse", 8, 1, true });
	{
		FlatPassCase c = { "Cube_tetra_n8_k0_nonsym", "cube", "tetra", "", 8, 0 };
		c.Nonsymmetric = true;
		cases.push_back(c);
	}
	{
		FlatPassCase c = { "Cube_tetra_n8_k0_cp10", "cube", "tetra", "", 8, 0 };
		c.CoarseningProlong = 10;
		cases.push_back(c);
	}
#endif
	return cases;
}

INSTANTIATE_TEST_SUITE_P(PassProlongationsProgram, PassProlongationsProgramTest, ::testing::ValuesIn(PassProlongationsCases()),
	[](const ::testing::TestParamInfo<FlatPassCase>& info) { return info.param.Name; });
