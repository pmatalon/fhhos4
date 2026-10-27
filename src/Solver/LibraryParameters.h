#pragma once
#include <optional>
#include "fhhos4.h"
#include "../Utils/Utils.h"
#include "Multigrid/SmootherFactory.h"
using namespace std;

// The parameters of the library fhhos4 (fhhos4::Solver) and their values in the program, in both directions: from the
// library's to the program's in the library (library/Solver.cpp), the other way around in the program's -s libuamg
// (Solver/LibraryUAMG.cpp). Only these two files include it: it includes the public header of the library.
class LibraryParameters
{
private:
	template <typename LibraryValue, typename ProgramValue>
	using Table = vector<pair<LibraryValue, ProgramValue>>;

	// The smoothers of the library, all in the program; Default: the adjoint of the pre-smoother as post-smoother (the
	// library chooses its default pre-smoother itself)
	static const Table<fhhos4::Smoother, SmootherType>& Smoothers()
	{
		static const Table<fhhos4::Smoother, SmootherType> smoothers = {
			{ fhhos4::Smoother::Default,                        SmootherType::Default },
			{ fhhos4::Smoother::BlockGaussSeidel,               SmootherType::BlockGaussSeidel },
			{ fhhos4::Smoother::BackwardBlockGaussSeidel,       SmootherType::BackwardBlockGaussSeidel },
			{ fhhos4::Smoother::SymmetricBlockGaussSeidel,      SmootherType::SymmetricBlockGaussSeidel },
			{ fhhos4::Smoother::HybridBlockGaussSeidel,         SmootherType::HybridBlockGaussSeidel },
			{ fhhos4::Smoother::HybridBackwardBlockGaussSeidel, SmootherType::HybridBackwardBlockGaussSeidel },
			{ fhhos4::Smoother::DampedBlockJacobi,              SmootherType::DampedBlockJacobi },
			{ fhhos4::Smoother::BlockChebyshevFirstKind,        SmootherType::BlockChebyshevFirstKind },
			{ fhhos4::Smoother::BlockChebyshevFourthKind,       SmootherType::BlockChebyshevFourthKind }
		};
		return smoothers;
	}

	static const Table<fhhos4::Cycle, CycleType>& Cycles()
	{
		static const Table<fhhos4::Cycle, CycleType> cycles = {
			{ fhhos4::Cycle::V, CycleType::V },
			{ fhhos4::Cycle::W, CycleType::W },
			{ fhhos4::Cycle::K, CycleType::K }
		};
		return cycles;
	}

	static const Table<fhhos4::PCoarsening, P_CoarsStgy>& PCoarsenings()
	{
		static const Table<fhhos4::PCoarsening, P_CoarsStgy> pCoarsenings = {
			{ fhhos4::PCoarsening::Minus2,    P_CoarsStgy::Minus2 },
			{ fhhos4::PCoarsening::Minus1,    P_CoarsStgy::Minus1 },
			{ fhhos4::PCoarsening::DivideBy2, P_CoarsStgy::DivideBy2 },
			{ fhhos4::PCoarsening::ToZero,    P_CoarsStgy::DirectToLow }
		};
		return pCoarsenings;
	}

	// The prolongations of the h-coarsening of the library, the coarsening prolongations of the program
	// (-coarsening-prolong; the multigrid prolongation chains them)
	static const Table<fhhos4::Prolongation, UAMGProlongation>& Prolongations()
	{
		static const Table<fhhos4::Prolongation, UAMGProlongation> prolongations = {
			{ fhhos4::Prolongation::Reconstruction,    UAMGProlongation::ReconstructSmoothedTraceOrInject },
			{ fhhos4::Prolongation::HarmonicExtension, UAMGProlongation::HarmonicExtensionOrInject }
		};
		return prolongations;
	}

	// The polynomial spaces of the library, the program's -poly-space (UsePolynomialSpaceQ)
	static const Table<fhhos4::PolynomialSpace, bool>& PolynomialSpaces()
	{
		static const Table<fhhos4::PolynomialSpace, bool> polynomialSpaces = {
			{ fhhos4::PolynomialSpace::P, false },
			{ fhhos4::PolynomialSpace::Q, true }
		};
		return polynomialSpaces;
	}

	// The coarse solvers of the library (Default: chosen by the library) and their codes in the program
	// (-coarse-solver)
	static const Table<fhhos4::CoarseSolver, string>& CoarseSolvers()
	{
		static const Table<fhhos4::CoarseSolver, string> coarseSolvers = {
			{ fhhos4::CoarseSolver::Cholesky, "ch" },
			{ fhhos4::CoarseSolver::LU,       "lu" },
			{ fhhos4::CoarseSolver::MUMPS,    "mumpsch" },
			{ fhhos4::CoarseSolver::MUMPS,    "mumpslu" }
		};
		return coarseSolvers;
	}

	// The program's value of a library's value; error if the table doesn't have it (a value of the enum out of its
	// range, by a cast)
	template <typename LibraryValue, typename ProgramValue>
	static ProgramValue ToProgram(const Table<LibraryValue, ProgramValue>& table, LibraryValue value, const string& parameter)
	{
		for (const auto& [l, p] : table)
		{
			if (l == value)
				return p;
		}
		Utils::FatalError("fhhos4::Solver: invalid " + parameter + " (" + to_string((int)value) + ").");
		return ProgramValue();
	}

	// The library's value of a program's value; none if the library doesn't have it
	template <typename LibraryValue, typename ProgramValue>
	static optional<LibraryValue> ToLibrary(const Table<LibraryValue, ProgramValue>& table, const ProgramValue& value)
	{
		for (const auto& [l, p] : table)
		{
			if (p == value)
				return l;
		}
		return nullopt;
	}

public:
	//--------------- Library -> program ---------------//

	static SmootherType ProgramSmoother(fhhos4::Smoother smoother)
	{
		return ToProgram(Smoothers(), smoother, "smoother");
	}

	static CycleType ProgramCycle(fhhos4::Cycle cycle)
	{
		return ToProgram(Cycles(), cycle, "Cycle");
	}

	static P_CoarsStgy ProgramPCoarsening(fhhos4::PCoarsening pCoarsening)
	{
		return ToProgram(PCoarsenings(), pCoarsening, "PCoarsening");
	}

	static UAMGProlongation ProgramProlongation(fhhos4::Prolongation prolongation)
	{
		return ToProgram(Prolongations(), prolongation, "Prolongation");
	}

	// UsePolynomialSpaceQ
	static bool ProgramPolynomialSpace(fhhos4::PolynomialSpace polynomialSpace)
	{
		return ToProgram(PolynomialSpaces(), polynomialSpace, "PolynomialSpace");
	}

	// The code of the coarse solver (not Default); MUMPS: its factorization for the symmetry of the matrix
	static string ProgramCoarseSolver(fhhos4::CoarseSolver coarseSolver, bool symmetric)
	{
		if (coarseSolver == fhhos4::CoarseSolver::MUMPS)
			return symmetric ? "mumpsch" : "mumpslu";
		return ToProgram(CoarseSolvers(), coarseSolver, "CoarseSolver");
	}

	//--------------- Program -> library ---------------//

	static fhhos4::Smoother LibrarySmoother(SmootherType smoother)
	{
		optional<fhhos4::Smoother> librarySmoother = ToLibrary(Smoothers(), smoother);
		if (!librarySmoother)
			Utils::FatalError("The library fhhos4 has no smoother '" + SmootherFactory::Code(smoother) + "' (-smoothers).");
		return *librarySmoother;
	}

	static fhhos4::Cycle LibraryCycle(CycleType cycle)
	{
		return *ToLibrary(Cycles(), cycle);
	}

	static fhhos4::PCoarsening LibraryPCoarsening(P_CoarsStgy pCoarsening)
	{
		return *ToLibrary(PCoarsenings(), pCoarsening);
	}

	static fhhos4::Prolongation LibraryProlongation(UAMGProlongation coarseningProlongation)
	{
		optional<fhhos4::Prolongation> prolongation = ToLibrary(Prolongations(), coarseningProlongation);
		if (!prolongation)
			Utils::FatalError("The library fhhos4 has no coarsening prolongation " + to_string((unsigned)coarseningProlongation) + " (-coarsening-prolong): " + to_string((unsigned)UAMGProlongation::ReconstructSmoothedTraceOrInject) + " or " + to_string((unsigned)UAMGProlongation::HarmonicExtensionOrInject) + ".");
		return *prolongation;
	}

	static fhhos4::PolynomialSpace LibraryPolynomialSpace(bool usePolynomialSpaceQ)
	{
		return *ToLibrary(PolynomialSpaces(), usePolynomialSpaceQ);
	}

	static fhhos4::CoarseSolver LibraryCoarseSolver(const string& coarseSolverCode)
	{
		optional<fhhos4::CoarseSolver> coarseSolver = ToLibrary(CoarseSolvers(), coarseSolverCode);
		if (!coarseSolver)
			Utils::FatalError("The library fhhos4 has no coarse solver '" + coarseSolverCode + "' (-coarse-solver): ch, lu, mumpsch or mumpslu.");
		return *coarseSolver;
	}
};
