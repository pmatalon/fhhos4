#pragma once
#include "Program/Program_Diffusion_HHO.h"
#include "ProgramArgumentsDefaults.h"
#include <cmath>
#include <stdexcept>
#include <vector>

namespace fhhos4_tests
{
	// Mirrors the side effects of ProgramDim<Dim>::Start() (src/Program.h) that the CLI relies
	// on before dispatching to a Program_*::Execute(): a lot of mesh-construction code (the GMSH
	// mesher, polyhedral coarsening) reads the global Utils::ProgramArgs directly instead of the
	// args reference passed around, so skipping this would silently run with stale/leftover
	// settings from whatever test happened to run previously in this process.
	template <int Dim>
	inline void SyncGlobalProgramState(const ProgramArguments& args)
	{
		Utils::ProgramArgs = args;
		Mesh<Dim>::SetDirectories();
		GMSHMesh<Dim>::GMSHLogEnabled = args.Actions.GMSHLogEnabled;
		GMSHMesh<Dim>::UseCache = args.Actions.UseCache;
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
		case 2:
			SyncGlobalProgramState<2>(args);
			Program_Diffusion_HHO<2>::Execute(args, &results);
			break;
		case 3:
			SyncGlobalProgramState<3>(args);
			Program_Diffusion_HHO<3>::Execute(args, &results);
			break;
		default:
			throw std::runtime_error("RunDiffusionHHO: unsupported dimension " + std::to_string(args.Problem.Dimension));
		}
		return results;
	}

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
