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

	// The coarse solvers of the library (Default: chosen by the library) and their codes in the program
	// (-coarse-solver)
	static const Table<fhhos4::CoarseSolver, string>& CoarseSolvers()
	{
		static const Table<fhhos4::CoarseSolver, string> coarseSolvers = {
			{ fhhos4::CoarseSolver::Cholesky, "ch" },
			{ fhhos4::CoarseSolver::LU,       "lu" }
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

	// The code of the coarse solver (not Default)
	static string ProgramCoarseSolver(fhhos4::CoarseSolver coarseSolver)
	{
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

	static fhhos4::CoarseSolver LibraryCoarseSolver(const string& coarseSolverCode)
	{
		optional<fhhos4::CoarseSolver> coarseSolver = ToLibrary(CoarseSolvers(), coarseSolverCode);
		if (!coarseSolver)
			Utils::FatalError("The library fhhos4 has no coarse solver '" + coarseSolverCode + "' (-coarse-solver): ch or lu.");
		return *coarseSolver;
	}
};
