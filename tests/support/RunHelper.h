#pragma once
#include "Program.h"
#include "Program/Program_Diffusion_HHO.h"
#include "Program/Program_BiHarmonic_HHO.h"
#include "ProgramArgumentsDefaults.h"
#include <cmath>
#include <ostream>
#include <stdexcept>
#include <vector>

namespace fhhos4_tests
{
	// Expected iteration count for the mesh size N: the current code's, and the paper's for
	// comparison (they can differ a little, e.g. when the meshes are built by GMSH).
	struct ExpectedIterations
	{
		int N;
		int Iterations;
		int PaperIterations;
	};

	// Readable parameter values in the test names
	inline std::ostream& operator<<(std::ostream& os, const ExpectedIterations& e)
	{
		return os << "N" << e.N << ":" << e.Iterations;
	}

	// Sets the global state that ProgramDim<Dim>::Start() (src/Program.cpp) sets before dispatching
	// to a Program_*::Execute(): a lot of mesh-construction code (the GMSH mesher, polyhedral
	// coarsening) reads the global Utils::ProgramArgs directly instead of the args reference
	// passed around, so skipping this would silently run with stale/leftover settings from
	// whatever test happened to run previously in this process.
	template <int Dim>
	inline void SyncGlobalProgramState(const ProgramArguments& args)
	{
		ProgramDim<Dim>::InitGlobalState(args);
	}

	// Applies the same argument-defaulting cascade as the CLI (main.cpp), then runs the
	// diffusion HHO problem in-process and returns the iteration count / L2 error.
	// A documented divergence (Utils::FatalError) terminates the process (exit(EXIT_FAILURE));
	// callers that expect that must wrap the call in GTest's ASSERT_EXIT/EXPECT_EXIT.
	// Pass defaultCycle = false when the test sets the MG cycle itself (the CLI equivalent of
	// -cycle); otherwise the defaults overwrite the pre/post-smoothing iterations.
	inline ProgramResults RunDiffusionHHO(ProgramArguments args, bool defaultCycle = true)
	{
		ApplyProgramArgumentDefaults(args, true, true, defaultCycle);

		ProgramResults results;
		switch (args.Problem.Dimension)
		{
#ifdef ENABLE_2D
		case 2:
			SyncGlobalProgramState<2>(args);
			Program_Diffusion_HHO<2>::Execute(args, &results);
			break;
#endif // ENABLE_2D
#ifdef ENABLE_3D
		case 3:
			SyncGlobalProgramState<3>(args);
			Program_Diffusion_HHO<3>::Execute(args, &results);
			break;
#endif // ENABLE_3D
		default:
			throw std::runtime_error("RunDiffusionHHO: unsupported (or disabled at compile time) dimension " + std::to_string(args.Problem.Dimension));
		}
		return results;
	}

	// Same as RunDiffusionHHO, for the biharmonic problem (-pb bihar): returns the iteration
	// count of the biharmonic solver and the L2 error of the solution.
	inline ProgramResults RunBiHarmonicHHO(ProgramArguments args)
	{
		args.Problem.Equation = EquationType::BiHarmonic;
		ApplyProgramArgumentDefaults(args);

		ProgramResults results;
		switch (args.Problem.Dimension)
		{
#ifdef ENABLE_2D
		case 2:
			SyncGlobalProgramState<2>(args);
			Program_BiHarmonic_HHO<2>::Execute(args, &results);
			break;
#endif // ENABLE_2D
#ifdef ENABLE_3D
		case 3:
			SyncGlobalProgramState<3>(args);
			Program_BiHarmonic_HHO<3>::Execute(args, &results);
			break;
#endif // ENABLE_3D
		default:
			throw std::runtime_error("RunBiHarmonicHHO: unsupported (or disabled at compile time) dimension " + std::to_string(args.Problem.Dimension));
		}
		return results;
	}

	// Sequential execution (the CLI's -threads 1) while in scope. Needed where the result depends
	// on the thread scheduling, e.g. the agglomeration coarsening (-cs n).
	class SequentialExecution
	{
	public:
		SequentialExecution() { Parallelism::SetNThreads(1); }
		~SequentialExecution() { Parallelism::SetNThreads(0); } // 0: back to the automatic default
	};

	// Least-squares slope of log(errors) vs. log(h): the empirical convergence order.
	inline double EstimateConvergenceOrder(const std::vector<double>& h, const std::vector<double>& errors)
	{
		if (h.size() != errors.size() || h.size() < 2)
			throw std::runtime_error("EstimateConvergenceOrder: need at least 2 matching (h, error) points");

		size_t n = h.size();
		double sumX = 0, sumY = 0, sumXY = 0, sumXX = 0;
		for (size_t i = 0; i < n; i++)
		{
			double x = std::log(h[i]);
			double y = std::log(errors[i]);
			sumX += x;
			sumY += y;
			sumXY += x * y;
			sumXX += x * x;
		}
		double denom = static_cast<double>(n) * sumXX - sumX * sumX;
		return (static_cast<double>(n) * sumXY - sumX * sumY) / denom;
	}
}
