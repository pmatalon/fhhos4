#pragma once
#include "ProgramArguments.h"
#include "Utils/Utils.h"
using namespace std;

// Prints the message as a fatal argument error and terminates the process.
// Shared by main()'s CLI parsing and by ApplyProgramArgumentDefaults() below.
inline void argument_error(string msg)
{
	cout << Utils::BeginRed << "Argument error: " << msg << Utils::EndColor << endl;
	cout << "------------------------- FAILURE -------------------------" << endl;
	exit(EXIT_FAILURE);
}

// The defaults of U-AMG (-s uamg, fcguamg, libuamg...) for the face degree k: part of ApplyProgramArgumentDefaults(),
// also used by the library fhhos4_AMG (library/Solver.cpp). The booleans: see ApplyProgramArgumentDefaults().
inline void ApplyUncondensedAMGDefaults(MultigridArguments& mg, int k, bool defaultCoarseOperator = true, bool defaultCycle = true, bool defaultHPCoarseningStgy = true)
{
	if (mg.ProlongationCode == 0)
		mg.UAMGMultigridProlong = UAMGProlongation::ChainedCoarseningProlongations;
	else
		mg.UAMGMultigridProlong = static_cast<UAMGProlongation>(mg.ProlongationCode);
	mg.UAMGFaceProlong = mg.FaceProlongationCode == 0 ? UAMGFaceProlongation::BoundaryAggregatesInteriorAverage : static_cast<UAMGFaceProlongation>(mg.FaceProlongationCode);
	mg.UAMGCoarseningProlong = mg.CoarseningProlongationCode == 0 ? UAMGProlongation::ReconstructSmoothedTraceOrInject : static_cast<UAMGProlongation>(mg.CoarseningProlongationCode);

	if (mg.H_CS == H_CoarsStgy::None)
		mg.H_CS = H_CoarsStgy::MultiplePairwiseAggregation;

	if ((mg.H_CS == H_CoarsStgy::MultiplePairwiseAggregation ||
		 mg.H_CS == H_CoarsStgy::MultipleAgglomerationCoarseningByFaceNeighbours)
		&& mg.CoarseningFactor == 0)
		mg.CoarseningFactor = 3.8;

	if (defaultCoarseOperator)
		mg.UseGalerkinOperator = true;

	if (defaultCycle)
		mg.CycleLetter = 'K';

	// k >= 1: p-levels down to k=0, then the h-coarsening of the paper, designed for k=0. The h-coarsening of the
	// degree-k blocks (-hp-cs h) only transfers their higher modes by plain aggregation: never faster in the tests.
	if (defaultHPCoarseningStgy && k >= 1)
		mg.HP_CS = HP_CoarsStgy::P_then_H;
}

// Resolves every sentinel/"default" value left in a ProgramArguments (Dimension == -1,
// SolverCode/Mesher/MeshCode == "default", etc.) into a concrete, consistent configuration,
// exactly as main()'s CLI parsing does after getopt. Callers that build a ProgramArguments
// directly (e.g. tests) can call this instead of going through the CLI to get the same
// defaulting behavior.
//
// The boolean flags mirror whether the corresponding CLI option was explicitly passed by the
// user (main() tracks them while parsing getopt); left at their default of `true` (= "not
// explicitly set"), the same default-inference cascade as the CLI's own defaults applies.
void ApplyProgramArgumentDefaults(ProgramArguments& args,
	bool defaultRelativeCellPolyDegree = true,
	bool defaultTol2 = true,
	bool defaultCycle = true,
	bool defaultCoarseOperator = true,
	bool defaultCoarseSolver = true,
	bool defaultHPCoarseningStgy = true);
