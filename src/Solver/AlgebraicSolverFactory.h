#pragma once
#include <functional>
#include "Direct/EigenLU.h"
#include "Direct/EigenCholesky.h"
#include "Krylov/ConjugateGradient.h"
#include "Krylov/FlexibleConjugateGradient.h"
#include "Multigrid/UncondensedAMG/UncondensedAMG.h"
#include "Multigrid/AggregAMG/AggregAMG.h"
#include "FixedPoint/BlockJacobi.h"
#include "Krylov/EigenCG.h"
#include "Krylov/EigenBiCGSTAB.h"
#include "Multigrid/AggregAMG/HighOrderAggregAMG.h"
#include "Multigrid/AGMG.h"
using namespace std;

// The solvers that only need the matrix of the linear system (and, for U-AMG, the sizes of the polynomial bases): no
// mesh, no discretization. Used by SolverFactory, and by the library fhhos4_AMG (library/), which must not depend on
// the meshes and discretizations.
class AlgebraicSolverFactory
{
public:
	// Creates the coarse solver of a multigrid
	using CoarseSolverFactory = function<Solver*(const ProgramArguments&)>;

	// The coarse solvers among the solvers of this factory
	static CoarseSolverFactory AlgebraicCoarseSolvers(int dim, const ExportModule& out)
	{
		return [dim, out](const ProgramArguments& args) { return CreateSolver(args, dim, 1, out, AlgebraicCoarseSolvers(dim, out)); };
	}

	static Solver* CreateSolver(const ProgramArguments& args, int dim, int blockSize, const ExportModule& out, const CoarseSolverFactory& createCoarseSolver)
	{
		Solver* solver = nullptr;
		string GaussSeidelRelaxation1 = "The relaxation parameter of Gauss-Seidel is 1. Delete -relax arg to remove this warning, or use -s sor.";

		if (args.Solver.SolverCode.compare("lu") == 0)
			solver = new EigenSparseLU();
		else if (args.Solver.SolverCode.compare("ch") == 0)
			solver = new EigenSparseCholesky();
		else if (args.Solver.SolverCode.compare("j") == 0)
			solver = new BlockJacobi(1, args.Solver.RelaxationParameter);
		else if (args.Solver.SolverCode.compare("gs") == 0)
		{
			if (args.Solver.RelaxationParameter != 1)
				Utils::Warning(GaussSeidelRelaxation1);
			solver = new GaussSeidel(Direction::Forward);
		}
		else if (args.Solver.SolverCode.compare("rgs") == 0)
		{
			if (args.Solver.RelaxationParameter != 1)
				Utils::Warning(GaussSeidelRelaxation1);
			solver = new GaussSeidel(Direction::Backward);
		}
		else if (args.Solver.SolverCode.compare("sgs") == 0)
		{
			if (args.Solver.RelaxationParameter != 1)
				Utils::Warning(GaussSeidelRelaxation1);
			solver = new GaussSeidel(Direction::Symmetric);
		}
		else if (args.Solver.SolverCode.compare("sor") == 0)
			solver = new BlockSOR(1, args.Solver.RelaxationParameter, Direction::Forward);
		else if (args.Solver.SolverCode.compare("rsor") == 0)
			solver = new BlockSOR(1, args.Solver.RelaxationParameter, Direction::Backward);
		else if (args.Solver.SolverCode.compare("ssor") == 0)
			solver = new BlockSOR(1, args.Solver.RelaxationParameter, Direction::Symmetric);
		else if (args.Solver.SolverCode.compare("bj") == 0)
			solver = new BlockJacobi(blockSize, args.Solver.RelaxationParameter);
		else if (args.Solver.SolverCode.compare("bj23") == 0)
			solver = new BlockJacobi(blockSize, 2.0 / 3.0);
		else if (args.Solver.SolverCode.compare("bgs") == 0)
		{
			if (args.Solver.RelaxationParameter != 1)
				Utils::Warning(GaussSeidelRelaxation1);
			solver = new BlockSOR(blockSize, 1, Direction::Forward);
		}
		else if (args.Solver.SolverCode.compare("rbgs") == 0)
		{
			if (args.Solver.RelaxationParameter != 1)
				Utils::Warning(GaussSeidelRelaxation1);
			solver = new BlockSOR(blockSize, 1, Direction::Backward);
		}
		else if (args.Solver.SolverCode.compare("sbgs") == 0)
		{
			if (args.Solver.RelaxationParameter != 1)
				Utils::Warning(GaussSeidelRelaxation1);
			solver = new BlockSOR(blockSize, 1, Direction::Symmetric);
		}
		else if (args.Solver.SolverCode.compare("bsor") == 0)
			solver = new BlockSOR(blockSize, args.Solver.RelaxationParameter, Direction::Forward);
		else if (args.Solver.SolverCode.compare("rbsor") == 0)
			solver = new BlockSOR(blockSize, args.Solver.RelaxationParameter, Direction::Backward);
		else if (args.Solver.SolverCode.compare("sbsor") == 0)
			solver = new BlockSOR(blockSize, args.Solver.RelaxationParameter, Direction::Symmetric); 
		else if (args.Solver.SolverCode.compare("cg") == 0)
			solver = new ConjugateGradient();
		else if (args.Solver.SolverCode.compare("eigencg") == 0)
			solver = new EigenCG();
		else if (args.Solver.SolverCode.compare("bicgstab") == 0)
			solver = new EigenBiCGSTAB();
		else if (args.Solver.SolverCode.compare("agmg") == 0)
			solver = new AGMG();
		else if (args.Solver.SolverCode.compare("aggregamg") == 0)
		{
			AggregAMG* mg = new AggregAMG(blockSize, 0.25, args.Solver.MG.Levels);
			SetMultigridParameters(mg, args, dim, blockSize, out, createCoarseSolver);
			mg->UseGalerkinOperator = 1;
			solver = mg;
		}
		else if (args.Solver.SolverCode.compare("hoaggregamg") == 0)
		{
			HighOrderAggregAMG* mg = new HighOrderAggregAMG(blockSize, 0.25, 2.0 / 3.0, 2);
			SetMultigridParameters(mg, args, dim, blockSize, out, createCoarseSolver);
			mg->UseGalerkinOperator = 1;
			solver = mg;
		}
		else
			Utils::FatalError("Unknown solver '" + args.Solver.SolverCode + "' or not applicable.");

		IterativeSolver* iterativeSolver = dynamic_cast<IterativeSolver*>(solver);
		if (iterativeSolver)
		{
			iterativeSolver->StoppingCrit = args.Solver.StoppingCrit;
			iterativeSolver->Tolerance = args.Solver.Tolerance;
			iterativeSolver->StagnationConvRate = args.Solver.StagnationConvRate;
			iterativeSolver->MaxIterations = args.Solver.MaxIterations;
			iterativeSolver->PrintIterationResults = args.Solver.PrintIterationResults;
		}

		return solver;
	}

	// U-AMG for a hybrid discretization of face degree faceDegree in dimension dim, with cellBlockSize (resp.
	// faceBlockSize) DoFs per cell (resp. face)
	static UncondensedAMG* CreateUncondensedAMG(const ProgramArguments& args, int dim, int faceDegree, int cellBlockSize, int faceBlockSize, int blockSize, const ExportModule& out, const CoarseSolverFactory& createCoarseSolver)
	{
		UncondensedAMG* mg = new UncondensedAMG(dim, faceDegree, cellBlockSize, faceBlockSize, 0.25, args.Solver.MG.UAMGFaceProlong, args.Solver.MG.UAMGCoarseningProlong, args.Solver.MG.UAMGMultigridProlong, args.Solver.MG.Levels);
		SetMultigridParameters(mg, args, dim, blockSize, out, createCoarseSolver);
		mg->CoarsePolyDegree = 0;
		mg->ManageAnisotropy = args.Solver.MG.ManageAnisotropy;
		return mg;
	}

	static void SetMultigridParameters(Multigrid* mg, const ProgramArguments& args, int dim, int blockSize, const ExportModule& out, const CoarseSolverFactory& createCoarseSolver)
	{
		mg->Out = ExportModule(out);
		mg->MatrixMaxSizeForCoarsestLevel = args.Solver.MG.MatrixMaxSizeForCoarsestLevel;
		mg->Cycle = args.Solver.MG.CycleLetter;
		mg->WLoops = args.Solver.MG.WLoops;
		mg->UseGalerkinOperator = args.Solver.MG.UseGalerkinOperator;
		mg->PreSmootherCode = args.Solver.MG.PreSmootherCode;
		mg->PostSmootherCode = args.Solver.MG.PostSmootherCode;
		mg->PreSmoothingIterations = args.Solver.MG.PreSmoothingIterations;
		mg->PostSmoothingIterations = args.Solver.MG.PostSmoothingIterations;
		mg->RelaxationParameter = args.Solver.RelaxationParameter;
		mg->BlockSizeForBlockSmoothers = blockSize;
		mg->CoarseLevelChangeSmoothingCoeff = args.Solver.MG.CoarseLevelChangeSmoothingCoeff;
		mg->CoarseLevelChangeSmoothingOperator = args.Solver.MG.CoarseLevelChangeSmoothingOperator;
		mg->HP_CS = args.Solver.MG.HP_CS;
		mg->H_CS = args.Solver.MG.H_CS;
		mg->P_CS = args.Solver.MG.P_CS;
		mg->FaceCoarseningStgy = args.Solver.MG.FaceCoarseningStgy;
		mg->BdryFaceCollapsing = args.Solver.MG.BoundaryFaceCollapsing;
		mg->NumberOfMeshes = args.Solver.MG.NumberOfMeshes;
		mg->CoarseningFactor = args.Solver.MG.CoarseningFactor;
		mg->ExportComponents = args.Actions.Export.MultigridComponents;
		mg->ExportIterationVectors = args.Actions.Export.MultigridIterationVectors;

		// Symmetric post-smoother
		if (mg->PostSmootherCode == "<symmetric smoother>")
		{
			if (mg->PreSmootherCode == "bj23" || mg->PreSmootherCode == "bj" || mg->PreSmootherCode == "sgs" || mg->PreSmootherCode == "sbgs") // symmetric smoothers
				mg->PostSmootherCode = mg->PreSmootherCode;
			else if (mg->PreSmootherCode == "gs")
				mg->PostSmootherCode = "rgs";
			else if (mg->PreSmootherCode == "bgs")
				mg->PostSmootherCode = "rbgs";
			else if (mg->PreSmootherCode == "hbgs")
				mg->PostSmootherCode = "hrbgs";
			else if (mg->PreSmootherCode == "ags" || mg->PreSmootherCode == "abgs" || mg->PreSmootherCode == "asor") // alternating smoothers
			{
				if (mg->CoarseLevelChangeSmoothingCoeff != 0)
					Utils::FatalError("The automatic determination of the symmetric version of the smoother '" + mg->PreSmootherCode + "' is not implemented when the cycle is variable.");

				if (mg->PreSmoothingIterations % 2 == 0) // even number of iterations
				{
					mg->PostSmootherCode = mg->PreSmootherCode;
				}
				else // odd number of iterations
				{
					mg->PostSmootherCode = "r" + mg->PreSmootherCode;
				}
			}
			else
				Utils::FatalError("The automatic determination of the symmetric version of the smoother '" + mg->PreSmootherCode + "' is not implemented. Please complete the post-smoother in the argument '-smoothers " + mg->PreSmootherCode + ",<post>'.");
		}

		// Coarse solver
		ProgramArguments argsCoarseSolver;
		argsCoarseSolver.Solver.SolverCode = args.Solver.MG.CoarseSolverCode;
		argsCoarseSolver.Solver.Tolerance = args.Solver.Tolerance;
		argsCoarseSolver.Solver.PrintIterationResults = false;
		argsCoarseSolver.Actions.Export.MultigridIterationVectors = args.Actions.Export.MultigridIterationVectors;
		if (Utils::EndsWith(args.Solver.MG.CoarseSolverCode, "aggregamg"))
		{
			argsCoarseSolver.Solver.MG.H_CS = H_CoarsStgy::AgglomerationCoarseningByFaceNeighbours;
			argsCoarseSolver.Solver.MG.CycleLetter = 'K';
		}
		else if (args.Solver.MG.CoarseSolverCode.compare("mg") == 0 || args.Solver.MG.CoarseSolverCode.compare("fcgmg") == 0)
		{
			argsCoarseSolver.Solver.MaxIterations = 1;
			argsCoarseSolver.Solver.MG.GMG_H_Prolong = args.Solver.MG.GMG_H_Prolong;
			argsCoarseSolver.Solver.MG.FaceProlongationCode = args.Solver.MG.FaceProlongationCode;
			argsCoarseSolver.Solver.MG.H_CS = args.Solver.MG.H_CS;
			argsCoarseSolver.Solver.MG.FaceCoarseningStgy = args.Solver.MG.FaceCoarseningStgy;
			argsCoarseSolver.Solver.MG.PreSmoothingIterations = 0;
			argsCoarseSolver.Solver.MG.PostSmoothingIterations = dim == 2 ? 3 : 6;
		}
		mg->CoarseSolver = createCoarseSolver(argsCoarseSolver);
	}
};
