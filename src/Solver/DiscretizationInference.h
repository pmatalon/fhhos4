#pragma once
#include <algorithm>
#include <cstdint>
#include <limits>
#include <numeric>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include "fhhos4.h"
#include "../Utils/PolynomialSpaceDimension.h"

// The discretization inferred from the blocks of the uncondensed system (ADR-0023): the parameters of fhhos4::Solver
// that are not set (Dimension, FaceDegree, CellDegree -1, PolynomialSpace::Inferred) determined from the DoFs per cell
// and per face, which the stored entries of A_TT and A_TF give. Header-only, without MPI nor PETSc: called by the setups
// of the library (library/Solver.cpp) and by the preconditioner for PETSc (library/petsc, into which it is compiled),
// each with its own gather of the analyses of the processes and its own names of the parameters in the messages. The
// messages start with "fhhos4: ".
namespace DiscretizationInference
{
	// The parameters of a discretization, as the members of fhhos4::Solver: -1 and PolynomialSpace::Inferred are not set
	struct Discretization
	{
		int Dimension = -1;
		fhhos4::PolynomialSpace PolynomialSpace = fhhos4::PolynomialSpace::Inferred;
		int FaceDegree = -1;
		int CellDegree = -1;
	};

	// The parameters of a Discretization, in the order of the proposals of InferDiscretization()
	enum class Parameter { Dimension, PolynomialSpace, FaceDegree, CellDegree };
	inline constexpr Parameter Parameters[] = { Parameter::Dimension, Parameter::PolynomialSpace, Parameter::FaceDegree, Parameter::CellDegree };

	// The names of the parameters in the messages: the members of fhhos4::Solver, or the options of the PETSc preconditioner
	struct ParameterNames
	{
		std::string Dimension, PolynomialSpace, FaceDegree, CellDegree, VectorComponents;
		std::string P, Q; // the values of PolynomialSpace
	};

	inline const ParameterNames MemberNames = { "Dimension", "PolynomialSpace", "FaceDegree", "CellDegree", "VectorComponents", "P", "Q" };

	//------------------------------------------------------------------------------------------------------------------//
	//                                   The DoFs per cell and face, from the blocks                                    //
	//------------------------------------------------------------------------------------------------------------------//

	// What the blocks of the cells of a process show (AnalyzeBlocks()); those of all the processes give the DoFs per cell
	// and per face (CombineAnalyses()). Integers only: gathered as bytes.
	struct BlockAnalysis
	{
		std::int64_t CellBlockSize = 0;    // the rows of its first cell in A_TT; 0: no cell
		std::int64_t FaceBoundaryGCD = 0;  // the GCD of the boundaries of the runs of consecutive columns of the first row of each cell in A_TF; 0: none
		std::int64_t InteriorFaceSeen = 0; // 1 if a column of A_TF is in the rows of two of its cells (a face between them)
	};

	// The analysis of the blocks of the cells of a process from their stored entries, not their values (the zeros inside
	// the blocks of a cell must be stored: exact zeros there are common, e.g. Legendre bases on a square): A_TT (full or
	// lower-only) and A_TF (its columns: the faces of the process, then its ghosts on several processes), row-major,
	// compressed. Reads the rows of its first and last cells in A_TT, and the first row of each cell in A_TF.
	// - The DoFs per cell: the first closed diagonal block of A_TT, the rows from 0 until no stored entry links them to the
	//   next one (none of them beyond the block, none of the next row within it: full or lower-only storage; with
	//   uncoupled components interleaved, the row of phi_0 e_x still reaches the last e_x of the cell). The last block,
	//   found the same way from the last row, must have its size (a zero stored between two cells merges them: none of
	//   the clients stores one).
	// - The DoFs per face, n_F: the first row of each cell in A_TF holds the DoFs of all its faces, in runs of faces
	//   numbered consecutively, whose boundaries are multiples of n_F. Their GCD is n_F as soon as a face is between two
	//   cells (reading g faces as one would take g faces shared by the same two cells), which a column in the rows of two
	//   cells shows, whatever n_F.
	template <typename TT, typename TF>
	BlockAnalysis AnalyzeBlocks(const Eigen::SparseCompressedBase<TT>& A_TT, const Eigen::SparseCompressedBase<TF>& A_TF)
	{
		static_assert(TT::IsRowMajor && TF::IsRowMajor, "DiscretizationInference::AnalyzeBlocks(): row-major blocks");
		std::int64_t n = A_TT.rows();
		auto size = [](std::int64_t rows, std::int64_t cols) { return std::to_string(rows) + " x " + std::to_string(cols); };
		if (A_TT.cols() != n || A_TF.rows() != n)
			throw fhhos4::Error("fhhos4: A_TT (" + size(n, A_TT.cols()) + ") must be square, with the rows of A_TF (" + size(A_TF.rows(), A_TF.cols()) + ").");
		if (!A_TT.isCompressed() || !A_TF.isCompressed())
			throw fhhos4::Error("fhhos4: A_TT and A_TF must be compressed.");
		BlockAnalysis analysis;
		if (n == 0)
			return analysis;

		const auto* rows = A_TT.outerIndexPtr();
		const auto* columns = A_TT.innerIndexPtr();
		auto empty = [&](std::int64_t i) { return rows[i] == rows[i + 1]; };
		auto first = [&](std::int64_t i) { return (std::int64_t)columns[rows[i]]; };
		auto last = [&](std::int64_t i) { return (std::int64_t)columns[rows[i + 1] - 1]; };
		std::int64_t end = 1; // the first cell: the rows [0, end)
		for (std::int64_t i = 0; i < n && (i < end || (!empty(i) && first(i) < end)); i++)
			end = std::max({ end, i + 1, empty(i) ? 0 : last(i) + 1 });
		std::int64_t begin = n - 1; // the last cell: [begin, n)
		for (std::int64_t i = n - 1; i >= 0 && (i >= begin || (!empty(i) && last(i) >= begin)); i--)
			begin = std::min({ begin, i, empty(i) ? n : first(i) });
		if (n - begin != end || n % end != 0)
			throw fhhos4::Error("fhhos4: the first cell of A_TT has " + std::to_string(end) + " DoFs and its last " + std::to_string(n - begin) + " (the closed diagonal blocks of its stored entries), for " + std::to_string(n) + " rows: the degrees must be the same in every cell, and the zeros inside the blocks of each cell stored.");
		analysis.CellBlockSize = end;

		// A face between two cells: a column in the first rows of two cells (usually among the first cells)
		const auto* faceRows = A_TF.outerIndexPtr();
		const auto* faceColumns = A_TF.innerIndexPtr();
		std::vector<bool> seen(A_TF.cols());
		for (std::int64_t r = 0; r < n && !analysis.InteriorFaceSeen; r += end)
		{
			for (std::int64_t p = faceRows[r]; p < (std::int64_t)faceRows[r + 1] && !analysis.InteriorFaceSeen; p++)
			{
				analysis.InteriorFaceSeen = seen[faceColumns[p]];
				seen[faceColumns[p]] = true;
			}
		}
		// The GCD of the boundaries of the runs of the first row of each cell, the cells in parallel where the header is
		// compiled with OpenMP (the library): the threads overlap the latency of these strided reads. A boundary that is a
		// multiple of the GCD so far, nearly all after the first cells, leaves it: a remainder, no gcd().
		#pragma omp parallel
		{
			std::int64_t gcd = 0;
			auto add = [&](std::int64_t boundary)
			{
				if (gcd == 0 || boundary % gcd != 0)
					gcd = std::gcd(gcd, boundary);
			};
			#pragma omp for schedule(static) nowait
			for (std::int64_t cell = 0; cell < n / end; cell++)
			{
				std::int64_t rowBegin = faceRows[cell * end], rowEnd = faceRows[cell * end + 1];
				for (std::int64_t p = rowBegin; p < rowEnd; p++)
				{
					std::int64_t previousEnd = p == rowBegin ? 0 : (std::int64_t)faceColumns[p - 1] + 1;
					if (p == rowBegin || faceColumns[p] != previousEnd) // a run starts: its first column, and the end of the previous one
					{
						add(faceColumns[p]);
						add(previousEnd);
					}
				}
				if (rowEnd > rowBegin) // the end of its last run
					add((std::int64_t)faceColumns[rowEnd - 1] + 1);
			}
			#pragma omp critical
			analysis.FaceBoundaryGCD = std::gcd(analysis.FaceBoundaryGCD, gcd);
		}
		return analysis;
	}

	// The DoFs per cell and per face of the blocks of all the processes, from their analyses: the processes with cells
	// must agree on the DoFs per cell (a process without cells has no opinion), and one of them must have a face between
	// two of its cells, otherwise the DoFs per face are undetermined (e.g. one cell on each process). The DoFs of the
	// cells and faces of each process must then be multiples of them: the checks of the setup.
	inline std::pair<int, int> CombineAnalyses(const std::vector<BlockAnalysis>& processes, const ParameterNames& names)
	{
		std::int64_t largest = 0, smallest = std::numeric_limits<std::int64_t>::max(), gcd = 0;
		bool interior = false;
		for (const BlockAnalysis& process : processes)
		{
			if (process.CellBlockSize > 0)
			{
				largest = std::max(largest, process.CellBlockSize);
				smallest = std::min(smallest, process.CellBlockSize);
			}
			gcd = std::gcd(gcd, process.FaceBoundaryGCD);
			interior = interior || process.InteriorFaceSeen != 0;
		}
		if (largest > 0 && smallest != largest)
			throw fhhos4::Error("fhhos4: the cells of the processes have different numbers of DoFs (" + std::to_string(smallest) + " and " + std::to_string(largest) + "): the degrees must be the same in every cell, and the zeros inside the blocks of each cell stored.");
		if (!interior)
			throw fhhos4::Error(std::string("fhhos4: the blocks don't give the DoFs per face: ") + (processes.size() > 1 ? "no process has a face between two of its cells" : "no face is between two cells") + " (a column of A_TF in the rows of two cells). Set " + names.Dimension + ", " + names.PolynomialSpace + ", " + names.FaceDegree + " and " + names.CellDegree + ".");
		if (largest > std::numeric_limits<int>::max() || gcd > std::numeric_limits<int>::max())
			throw fhhos4::Error("fhhos4: the blocks have " + std::to_string(largest) + " DoFs per cell and " + std::to_string(gcd) + " per face: too many.");
		return { (int)largest, (int)gcd };
	}

	//------------------------------------------------------------------------------------------------------------------//
	//                               The discretizations of the DoFs per cell and face                                  //
	//------------------------------------------------------------------------------------------------------------------//

	inline bool Q(const Discretization& d)
	{
		return d.PolynomialSpace == fhhos4::PolynomialSpace::Q;
	}

	inline int Value(const Discretization& d, Parameter p)
	{
		switch (p)
		{
			case Parameter::Dimension:       return d.Dimension;
			case Parameter::PolynomialSpace: return d.PolynomialSpace == fhhos4::PolynomialSpace::Inferred ? -1 : (int)d.PolynomialSpace;
			case Parameter::FaceDegree:      return d.FaceDegree;
			default:                         return d.CellDegree;
		}
	}

	// The degree k whose P^k (or Q^k) has n functions in m variables (their number increases with k), -1 if none
	inline int Degree(int m, bool q, std::int64_t n)
	{
		if (m == 1)
			return n >= 1 && n - 1 <= std::numeric_limits<int>::max() ? (int)(n - 1) : -1; // k + 1 functions
		for (int k = 0;; k++)
		{
			std::int64_t functions = ::PolynomialSpaceDimension(m, k, q);
			if (functions >= n)
				return functions == n ? k : -1;
		}
	}

	// The discretizations that the library takes (2D or 3D, P or Q, CellDegree k - 1, k or k + 1) whose cells and faces
	// have these DoFs (with the components), with the values of the parameters set: by face degree, dimension, space (P
	// first) and cell degree. Finite and exact: n_F / components DoFs per face give the face degree in each dimension
	// and space.
	inline std::vector<Discretization> Readings(int cellBlockSize, int faceBlockSize, int components, const Discretization& set)
	{
		std::vector<Discretization> readings;
		if (components < 1 || cellBlockSize < 1 || faceBlockSize < 1 || cellBlockSize % components != 0 || faceBlockSize % components != 0)
			return readings;
		for (int dimension : { 2, 3 })
		{
			for (fhhos4::PolynomialSpace space : { fhhos4::PolynomialSpace::P, fhhos4::PolynomialSpace::Q })
			{
				int k = Degree(dimension - 1, space == fhhos4::PolynomialSpace::Q, faceBlockSize / components);
				for (int offset : { -1, 0, 1 })
				{
					Discretization d{ dimension, space, k, k + offset };
					bool hasTheSetValues = true;
					for (Parameter p : Parameters)
						hasTheSetValues = hasTheSetValues && (Value(set, p) == -1 || Value(set, p) == Value(d, p));
					if (k >= 0 && d.CellDegree >= 0 && hasTheSetValues && ::PolynomialSpaceDimension(dimension, d.CellDegree, Q(d)) == cellBlockSize / components)
						readings.push_back(d);
				}
			}
		}
		auto order = [](const Discretization& d) { return std::tuple(d.FaceDegree, d.Dimension, d.PolynomialSpace, d.CellDegree); };
		std::sort(readings.begin(), readings.end(), [&](const Discretization& a, const Discretization& b) { return order(a) < order(b); });
		return readings;
	}

	// The degree of the p-level below one of degree k (0: the h-levels)
	inline int CoarserDegree(int k, fhhos4::PCoarsening pCoarsening)
	{
		switch (pCoarsening)
		{
			case fhhos4::PCoarsening::Minus2: return std::max(k - 2, 0);
			case fhhos4::PCoarsening::Minus1: return std::max(k - 1, 0);
			default:                          return 0;
		}
	}

	// What gives the solver of a discretization among those of the same DoFs per cell and face: whether it has p-levels
	// (k >= 1), and the DoFs per cell and face of its p-levels of degree >= 1 below the finest one (those of the finest
	// level and of the levels of degree 0 are the same). The rest of U-AMG doesn't depend on the dimension, the space nor
	// the degrees.
	inline std::vector<std::int64_t> SolverKey(const Discretization& d, fhhos4::PCoarsening pCoarsening)
	{
		std::vector<std::int64_t> key = { d.FaceDegree >= 1 };
		for (int degree = CoarserDegree(d.FaceDegree, pCoarsening); degree >= 1; degree = CoarserDegree(degree, pCoarsening))
		{
			key.push_back(::PolynomialSpaceDimension(d.Dimension, degree, Q(d)));
			key.push_back(::PolynomialSpaceDimension(d.Dimension - 1, degree, Q(d)));
		}
		return key;
	}

	// The number of different solvers of the discretizations
	inline int Solvers(const std::vector<Discretization>& discretizations, fhhos4::PCoarsening pCoarsening)
	{
		std::vector<std::vector<std::int64_t>> keys;
		for (const Discretization& d : discretizations)
			keys.push_back(SolverKey(d, pCoarsening));
		std::sort(keys.begin(), keys.end());
		return (int)(std::unique(keys.begin(), keys.end()) - keys.begin());
	}

	// The values of the parameter p among the discretizations, in increasing order
	inline std::vector<int> Values(const std::vector<Discretization>& discretizations, Parameter p)
	{
		std::vector<int> values;
		for (const Discretization& d : discretizations)
			values.push_back(Value(d, p));
		std::sort(values.begin(), values.end());
		values.erase(std::unique(values.begin(), values.end()), values.end());
		return values;
	}

	// The discretizations of the DoFs per cell and face read from the blocks (with the components), with the parameters
	// set, grouped by solver (InferDiscretization())
	struct Inference
	{
		int CellBlockSize = 0;
		int FaceBlockSize = 0;
		int Components = 0;
		Discretization Set;                        // the parameters set
		std::vector<Discretization> Readings;      // the discretizations of these DoFs with the parameters set (Readings())
		int Solvers = 0;                           // their different solvers: 1, resolved; 0, no reading; more, ambiguous
		Discretization Result;                     // resolved: the first reading
		std::vector<Parameter> NotNeeded;          // resolved: the parameters not set whose values differ among the readings (any gives the solver)
		Parameter Proposed = Parameter::Dimension; // ambiguous: the parameter to set, to one of its values among the readings
	};

	// The discretizations of cellBlockSize DoFs per cell and faceBlockSize per face (with the components), with the
	// parameters set (the others -1, PolynomialSpace::Inferred), grouped by solver: resolved if they have one solver (the
	// parameters it doesn't determine, from any reading: they don't change it); ambiguous if they have several (any
	// ambiguity is an error), the parameter proposed the first of Dimension, PolynomialSpace, FaceDegree and CellDegree
	// that, set to any of its values among the readings, leaves one solver (otherwise the first of several values).
	inline Inference InferDiscretization(int cellBlockSize, int faceBlockSize, int components, const Discretization& set, fhhos4::PCoarsening pCoarsening)
	{
		Inference inference;
		inference.CellBlockSize = cellBlockSize;
		inference.FaceBlockSize = faceBlockSize;
		inference.Components = components;
		inference.Set = set;
		inference.Readings = Readings(cellBlockSize, faceBlockSize, components, set);
		inference.Solvers = Solvers(inference.Readings, pCoarsening);
		if (inference.Solvers == 1)
		{
			inference.Result = inference.Readings.front();
			for (Parameter p : Parameters)
				if (Values(inference.Readings, p).size() > 1)
					inference.NotNeeded.push_back(p);
		}
		else if (inference.Solvers > 1)
		{
			auto resolves = [&](Parameter p)
			{
				for (int value : Values(inference.Readings, p))
				{
					std::vector<Discretization> withValue;
					std::copy_if(inference.Readings.begin(), inference.Readings.end(), std::back_inserter(withValue), [&](const Discretization& d) { return Value(d, p) == value; });
					if (Solvers(withValue, pCoarsening) != 1)
						return false;
				}
				return true;
			};
			auto several = [&](Parameter p) { return Values(inference.Readings, p).size() > 1; };
			auto proposed = std::find_if(std::begin(Parameters), std::end(Parameters), [&](Parameter p) { return several(p) && resolves(p); });
			if (proposed == std::end(Parameters))
				proposed = std::find_if(std::begin(Parameters), std::end(Parameters), several);
			inference.Proposed = *proposed;
		}
		return inference;
	}

	//------------------------------------------------------------------------------------------------------------------//
	//                                                     Messages                                                     //
	//------------------------------------------------------------------------------------------------------------------//

	inline const std::string& Name(Parameter p, const ParameterNames& names)
	{
		switch (p)
		{
			case Parameter::Dimension:       return names.Dimension;
			case Parameter::PolynomialSpace: return names.PolynomialSpace;
			case Parameter::FaceDegree:      return names.FaceDegree;
			default:                         return names.CellDegree;
		}
	}

	inline std::string ValueText(Parameter p, int value, const ParameterNames& names)
	{
		if (p == Parameter::PolynomialSpace)
			return value == (int)fhhos4::PolynomialSpace::Q ? names.Q : names.P;
		return std::to_string(value);
	}

	// "x", "x or y", "x, y or z"
	inline std::string Alternatives(const std::vector<std::string>& values)
	{
		std::string text;
		for (size_t i = 0; i < values.size(); i++)
			text += (i == 0 ? "" : i + 1 == values.size() ? " or " : ", ") + values[i];
		return text;
	}

	// "64 DoFs per cell and 9 per face (1 component)", and ", with Dimension 3," if it is set
	inline std::string Sizes(const Inference& inference, const ParameterNames& names, bool withSet)
	{
		std::string text = std::to_string(inference.CellBlockSize) + (inference.CellBlockSize == 1 ? " DoF" : " DoFs") + " per cell and " + std::to_string(inference.FaceBlockSize)
		                   + " per face (" + std::to_string(inference.Components) + (inference.Components == 1 ? " component)" : " components)");
		std::vector<std::string> set;
		for (Parameter p : Parameters)
			if (Value(inference.Set, p) != -1)
				set.push_back(Name(p, names) + " " + ValueText(p, Value(inference.Set, p), names));
		if (withSet && !set.empty())
		{
			text += ", with";
			for (size_t i = 0; i < set.size(); i++)
				text += (i == 0 ? " " : ", ") + set[i];
			text += ",";
		}
		return text;
	}

	// The error of an inference that doesn't resolve (Solvers != 1): the readings and the parameter proposed, or none
	inline std::string Message(const Inference& inference, const ParameterNames& names)
	{
		if (inference.Solvers == 0)
			return "fhhos4: " + Sizes(inference, names, true) + " fit no discretization of the library: P^l or Q^l in the cells and P^k or Q^k on the faces, in 2D or 3D, with l = k - 1, k or k + 1, the same in every cell and face. "
			       "These DoFs are read from the stored entries of A_TT and A_TF: the zeros inside the blocks of each cell must be stored. Check " + names.VectorComponents + ", or set " + names.Dimension + ", " + names.PolynomialSpace + ", " + names.FaceDegree + " and " + names.CellDegree + ".";
		std::string text = "fhhos4: " + Sizes(inference, names, true) + " fit discretizations that need different solvers:";
		for (const Discretization& d : inference.Readings)
			text += "\n  Dimension " + std::to_string(d.Dimension) + ", " + (Q(d) ? "Q" : "P") + ", face degree " + std::to_string(d.FaceDegree) + ", cell degree " + std::to_string(d.CellDegree);
		std::vector<std::string> values;
		for (int value : Values(inference.Readings, inference.Proposed))
			values.push_back(ValueText(inference.Proposed, value, names));
		return text + "\nSet " + Name(inference.Proposed, names) + " (" + Alternatives(values) + ").";
	}

	// What a resolved inference determined: the values of the parameters not set ("not needed": any value of the
	// readings gives the solver)
	inline std::string Description(const Inference& inference, const ParameterNames& names)
	{
		std::string text = "fhhos4: inferred from the blocks, " + Sizes(inference, names, false) + ":";
		const char* separator = " ";
		for (Parameter p : Parameters)
		{
			if (Value(inference.Set, p) != -1)
				continue;
			bool notNeeded = std::find(inference.NotNeeded.begin(), inference.NotNeeded.end(), p) != inference.NotNeeded.end();
			text += separator + Name(p, names) + " " + ValueText(p, Value(inference.Result, p), names) + (notNeeded ? " (not needed)" : "");
			separator = ", ";
		}
		return text;
	}
}
