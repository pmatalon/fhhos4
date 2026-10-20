#include <gtest/gtest.h>
#include <Eigen/SparseCholesky>
#include "ProgramArgumentsDefaults.h"
#include "Solver/AlgebraicSolverFactory.h"
#include "support/CartesianHybridSystem.h"
#include "support/PartitionedSystem.h"

using namespace fhhos4_tests;

namespace
{
	const Communicator& World()
	{
		return Communicator::Current();
	}

	// U-AMG preconditioning FCG, with the library's defaults (library/Solver.cpp): K-cycle, block Chebyshev of the fourth
	// kind of degree 2, Cholesky on the coarsest level
	struct UAMGWithFCG
	{
		unique_ptr<UncondensedAMG> Multigrid;
		unique_ptr<FlexibleConjugateGradient> FCG;

		UAMGWithFCG(int dim, int faceDegree, int cellBlockSize, int faceBlockSize, double tolerance)
		{
			ProgramArguments args;
			args.Solver.Tolerance = tolerance;
			args.Solver.MaxIterations = 200;
			args.Solver.PrintIterationResults = false;
			MultigridArguments& mg = args.Solver.MG;
			mg.Cycle = CycleType::K;
			mg.KCycleKrylov = KrylovMethod::FCG;
			mg.PreSmoother = SmootherType::BlockChebyshevFourthKind;
			mg.PostSmoother = SmootherType::Default;
			mg.PreSmoothingIterations = 2;
			mg.PostSmoothingIterations = 2;
			mg.CoarseningFactor = 3.8;
			mg.CoarseSolverCode = "ch";
			ApplyUncondensedAMGDefaults(mg, faceDegree, true, false, true);
			ExportModule out;
			Multigrid.reset(AlgebraicSolverFactory::CreateUncondensedAMG(args, dim, faceDegree, cellBlockSize, faceBlockSize, faceBlockSize, out, AlgebraicSolverFactory::AlgebraicCoarseSolvers(dim, out)));
			FCG = make_unique<FlexibleConjugateGradient>(1);
			FCG->ResidualRecomputation = 20;
			FCG->Precond = SolverPreconditioner(Multigrid.get());
			FCG->Tolerance = tolerance;
			FCG->MaxIterations = 200;
			FCG->PrintIterationResults = false;
		}
	};

	Vector DirectSolution(const SparseMatrix& S, const Vector& b)
	{
		ColMajorSparseMatrix colMajor = S;
		return Eigen::SimplicialLDLT<ColMajorSparseMatrix>(colMajor).solve(b);
	}
}

// The distributed U-AMG (each process aggregates its cells) on the Cartesian system partitioned in strips or blocks:
// the solution of the direct solver, in about as many iterations as on one process
TEST(DistributedUAMG, CartesianSystem)
{
	ostream silent(nullptr);
	ostream* savedLog = Utils::LogStream;
	Utils::LogStream = &silent;
	for (PartitionShape shape : { PartitionShape::Strips, PartitionShape::Blocks })
	{
		int n = 100;
		CartesianHybridSystem global(n, 0);
		PartitionedSystem system(global, n, shape, World());
		Vector b = Vector::Ones(system.Global.S.rows()); // the permuted system

		// One process: the sequential code on the whole system
		int sequentialIterations;
		Vector sequentialSolution;
		{
			CurrentCommunicator self(system.Self);
			UAMGWithFCG solver(2, 0, 1, 1, 1e-10);
			Vector ones = Vector::Ones(system.Global.A_TT.rows()), facesOnes = Vector::Ones(system.Global.S.rows());
			solver.FCG->Setup(system.Global.S, system.Global.A_TT, system.Global.A_TF, SparseMatrix(), SparseMatrix(), ones, facesOnes);
			sequentialSolution = static_cast<IterativeSolver&>(*solver.FCG).Solve(b);
			sequentialIterations = solver.FCG->IterationCount;
		}

		UAMGWithFCG solver(2, 0, 1, 1, 1e-10);
		Vector cellOnes = Vector::Ones(system.Cells->NOwned), faceOnes = Vector::Ones(system.Faces->NOwned);
		SparseMatrix A_TF = system.A_TF.LocalColumns();
		solver.FCG->Setup(system.S, system.A_TT, A_TF, system.A_TF.ColMap, cellOnes, faceOnes);
		Vector localB = b.segment(system.Faces->Begin, system.Faces->NOwned);
		Vector x = static_cast<IterativeSolver&>(*solver.FCG).Solve(localB);
		int iterations = solver.FCG->IterationCount;

		vector<double> local(x.data(), x.data() + x.rows());
		vector<double> all = World().AllGatherV(local);
		Vector solution = Eigen::Map<Vector>(all.data(), all.size());
		Vector exact = DirectSolution(system.Global.S, b);
		EXPECT_LE((solution - exact).norm(), 1e-8 * exact.norm()) << "shape " << (int)shape;
		if (World().Size() == 1) // the distributed code on one process: the sequential one, up to the product by the lower triangle of S
		{
			EXPECT_EQ(iterations, sequentialIterations);
			EXPECT_LE((solution - sequentialSolution).norm(), 1e-12 * exact.norm());
		}
		EXPECT_LE(iterations, sequentialIterations + 4) << "shape " << (int)shape << ", " << World().Size() << " processes";
		EXPECT_GE(solver.Multigrid->NumberOfLevels(), 3);
		if (World().IsRoot())
			cout << "\t" << World().Size() << " processes, " << (shape == PartitionShape::Strips ? "strips" : "blocks") << ": " << iterations << " iterations (" << sequentialIterations << " on one), " << solver.Multigrid->NumberOfLevels() << " levels" << endl;
	}
	Utils::LogStream = savedLog;
}
