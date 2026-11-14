// Tests of the discretization inferred from the blocks (src/Solver/DiscretizationInference.h, ADR-0023), without PETSc nor
// MPI: InferDiscretization() and CombineAnalyses() exhaustively, from the results of the analysis of the blocks (no
// matrix), and AnalyzeBlocks() on tiny matrices.
#include <gtest/gtest.h>
#include <algorithm>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include "Solver/DiscretizationInference.h"

using namespace DiscretizationInference;

namespace
{
	using Space = fhhos4::PolynomialSpace;
	constexpr fhhos4::PCoarsening PCoarsenings[] = { fhhos4::PCoarsening::Minus2, fhhos4::PCoarsening::Minus1, fhhos4::PCoarsening::ToZero };

	// The number of functions of P^k or Q^k in m variables (1 to 3), by their closed forms: independent of
	// PolynomialSpaceDimension(), both checked against the monomials counted (DimensionsOfTheSpaces)
	std::int64_t Functions(int m, std::int64_t k, bool q)
	{
		if (q)
			return m == 1 ? k + 1 : m == 2 ? (k + 1) * (k + 1) : (k + 1) * (k + 1) * (k + 1);
		return m == 1 ? k + 1 : m == 2 ? (k + 1) * (k + 2) / 2 : (k + 1) * (k + 2) * (k + 3) / 6;
	}

	// The DoFs per cell and per face of a discretization on nc components
	std::pair<int, int> BlockSizes(const Discretization& d, int nc)
	{
		bool q = d.PolynomialSpace == Space::Q;
		return { (int)(nc * Functions(d.Dimension, d.CellDegree, q)), (int)(nc * Functions(d.Dimension - 1, d.FaceDegree, q)) };
	}

	std::tuple<int, Space, int, int> Tuple(const Discretization& d)
	{
		return { d.Dimension, d.PolynomialSpace, d.FaceDegree, d.CellDegree };
	}

	std::string Text(const Discretization& d)
	{
		return "Dimension " + std::to_string(d.Dimension) + ", " + (d.PolynomialSpace == Space::Q ? "Q" : d.PolynomialSpace == Space::P ? "P" : "?") +
		       ", k = " + std::to_string(d.FaceDegree) + ", l = " + std::to_string(d.CellDegree);
	}

	// The configurations of the library up to the face degree maxDegree: 2D and 3D, P and Q, l = k - 1, k or k + 1
	std::vector<Discretization> Truths(int maxDegree)
	{
		std::vector<Discretization> truths;
		for (int dimension : { 2, 3 })
			for (Space space : { Space::P, Space::Q })
				for (int k = 0; k <= maxDegree; k++)
					for (int l = std::max(k - 1, 0); l <= k + 1; l++)
						truths.push_back({ dimension, space, k, l });
		return truths;
	}

	// The readings of these DoFs per cell and face by brute force: every dimension, space, face degree k up to the DoFs
	// per face, and cell degree k - 1 to k + 1
	std::vector<std::tuple<int, Space, int, int>> BruteForce(int cellBlockSize, int faceBlockSize, int nc)
	{
		std::vector<std::tuple<int, Space, int, int>> readings;
		for (int dimension : { 2, 3 })
			for (Space space : { Space::P, Space::Q })
				for (int k = 0; k <= faceBlockSize; k++)
					for (int l = std::max(k - 1, 0); l <= k + 1; l++)
						if (BlockSizes({ dimension, space, k, l }, nc) == std::pair(cellBlockSize, faceBlockSize))
							readings.push_back({ dimension, space, k, l });
		std::sort(readings.begin(), readings.end());
		return readings;
	}

	std::vector<std::tuple<int, Space, int, int>> Sorted(const std::vector<Discretization>& readings)
	{
		std::vector<std::tuple<int, Space, int, int>> tuples;
		for (const Discretization& d : readings)
			tuples.push_back(Tuple(d));
		std::sort(tuples.begin(), tuples.end());
		return tuples;
	}

	void SetValue(Discretization& d, Parameter p, int value)
	{
		switch (p)
		{
			case Parameter::Dimension:       d.Dimension = value; break;
			case Parameter::PolynomialSpace: d.PolynomialSpace = (Space)value; break;
			case Parameter::FaceDegree:      d.FaceDegree = value; break;
			case Parameter::CellDegree:      d.CellDegree = value; break;
		}
	}

	// The parameters of the truth in the subset (its bit i: Parameters[i]), the others not set
	Discretization Subset(const Discretization& truth, int subset)
	{
		Discretization set;
		for (int i = 0; i < 4; i++)
			if (subset & (1 << i))
				SetValue(set, Parameters[i], Value(truth, Parameters[i]));
		return set;
	}

	constexpr int DimensionAndDegrees = 0b1101; // Dimension, FaceDegree, CellDegree
}

// The closed forms of the test against the monomials counted, and PolynomialSpaceDimension() against them
TEST(DiscretizationInference, DimensionsOfTheSpaces)
{
	for (int m = 1; m <= 3; m++)
	{
		for (bool q : { false, true })
		{
			for (int k = 0; k <= 30; k++)
			{
				int monomials = 0;
				for (int a = 0; a <= k; a++)
					for (int b = 0; b <= (m >= 2 ? k : 0); b++)
						for (int c = 0; c <= (m >= 3 ? k : 0); c++)
							monomials += q || a + b + c <= k;
				EXPECT_EQ(Functions(m, k, q), monomials) << m << " variables, degree " << k << (q ? " Q" : " P");
				EXPECT_EQ(PolynomialSpaceDimension(m, k, q), monomials) << m << " variables, degree " << k << (q ? " Q" : " P");
			}
		}
	}
}

// Every configuration up to k = 10 (2D and 3D, P and Q, l = k - 1 to k + 1, 1 to 3 components), with the 3 p-coarsenings
// and the 16 subsets of the parameters set (18,432 cases): never no reading; resolved, the solver of the truth and the
// parameters set; ambiguous, the truth among the readings, the parameter proposed not set, which set to the truth's
// value resolves to the truth's solver in one step, Dimension proposed whenever it alone suffices; the readings of the
// DoFs per cell and face those of a brute force (a reading missed would be a silently wrong solver); Dimension and the
// degrees set without the space: resolved (from the sizes); VectorComponents forgotten (2 or 3 components read with 1):
// no reading
TEST(DiscretizationInference, Exhaustive)
{
	std::set<std::tuple<int, int, int>> sizesChecked;
	int cases = 0;
	for (const Discretization& truth : Truths(10))
	{
		for (int nc = 1; nc <= 3; nc++)
		{
			auto [cellBlockSize, faceBlockSize] = BlockSizes(truth, nc);
			std::string label = Text(truth) + ", " + std::to_string(nc) + " components";
			if (sizesChecked.insert({ cellBlockSize, faceBlockSize, nc }).second)
				EXPECT_EQ(Sorted(Readings(cellBlockSize, faceBlockSize, nc, Discretization())), BruteForce(cellBlockSize, faceBlockSize, nc)) << label;
			for (fhhos4::PCoarsening pCoarsening : PCoarsenings)
			{
				std::vector<std::int64_t> solver = SolverKey(truth, pCoarsening);
				for (int subset = 0; subset < 16; subset++, cases++)
				{
					Discretization set = Subset(truth, subset);
					std::string at = label + ", PCoarsening " + std::to_string((int)pCoarsening) + ", parameters set " + std::to_string(subset);
					Inference inference = InferDiscretization(cellBlockSize, faceBlockSize, nc, set, pCoarsening);
					ASSERT_GE(inference.Solvers, 1) << at;
					EXPECT_TRUE(std::any_of(inference.Readings.begin(), inference.Readings.end(), [&](const Discretization& d) { return Tuple(d) == Tuple(truth); })) << at;
					if (subset == DimensionAndDegrees)
						EXPECT_EQ(inference.Solvers, 1) << at;
					if (inference.Solvers == 1)
					{
						EXPECT_EQ(SolverKey(inference.Result, pCoarsening), solver) << at;
						for (Parameter p : Parameters)
							if (Value(set, p) != -1)
								EXPECT_EQ(Value(inference.Result, p), Value(set, p)) << at;
						continue;
					}
					Parameter proposed = inference.Proposed;
					EXPECT_EQ(Value(set, proposed), -1) << at;
					Discretization withProposed = set;
					SetValue(withProposed, proposed, Value(truth, proposed));
					Inference resolved = InferDiscretization(cellBlockSize, faceBlockSize, nc, withProposed, pCoarsening);
					EXPECT_EQ(resolved.Solvers, 1) << at;
					EXPECT_EQ(SolverKey(resolved.Result, pCoarsening), solver) << at;
					bool dimensionSuffices = Value(set, Parameter::Dimension) == -1 && Values(inference.Readings, Parameter::Dimension).size() > 1;
					for (int dimension : Values(inference.Readings, Parameter::Dimension))
					{
						Discretization withDimension = set;
						withDimension.Dimension = dimension;
						dimensionSuffices = dimensionSuffices && InferDiscretization(cellBlockSize, faceBlockSize, nc, withDimension, pCoarsening).Solvers == 1;
					}
					if (dimensionSuffices)
						EXPECT_TRUE(proposed == Parameter::Dimension) << at;
				}
				if (nc >= 2)
					EXPECT_EQ(InferDiscretization(cellBlockSize, faceBlockSize, 1, Subset(truth, DimensionAndDegrees), pCoarsening).Solvers, 0) << label << ": VectorComponents forgotten";
			}
		}
	}
	EXPECT_EQ(cases, 18432);
}

// The 112 common configurations (2D and 3D, P and Q, k <= 4, l = k - 1 to k + 1, scalar or with Dimension components),
// with only VectorComponents set: ambiguous, 2 with Minus2 (3D Q k=2 l=3, scalar and vector, which have the DoFs of 2D
// Q k=8 l=7), 10 with Minus1 (also 3D P k=1 l=1 and 2D Q k=2 l=1, 3D P k=1 l=2 and 2D P k=2 l=3), none with ToZero (the
// finest level the same, no other p-level)
TEST(DiscretizationInference, CommonConfigurations)
{
	std::map<fhhos4::PCoarsening, int> ambiguous;
	int configurations = 0;
	for (const Discretization& truth : Truths(4))
	{
		for (int nc : { 1, truth.Dimension })
		{
			configurations++;
			auto [cellBlockSize, faceBlockSize] = BlockSizes(truth, nc);
			for (fhhos4::PCoarsening pCoarsening : PCoarsenings)
				ambiguous[pCoarsening] += InferDiscretization(cellBlockSize, faceBlockSize, nc, Discretization(), pCoarsening).Solvers > 1;
		}
	}
	EXPECT_EQ(configurations, 112);
	EXPECT_EQ(ambiguous[fhhos4::PCoarsening::Minus2], 2);
	EXPECT_EQ(ambiguous[fhhos4::PCoarsening::Minus1], 10);
	EXPECT_EQ(ambiguous[fhhos4::PCoarsening::ToZero], 0);
	// Beyond k = 4: 2D Q k=6 l=5 and P k=6 l=7 (36 DoFs per cell and 7 per face): PolynomialSpace proposed
	Inference spaces = InferDiscretization(36, 7, 1, Discretization(), fhhos4::PCoarsening::Minus2);
	EXPECT_EQ(spaces.Solvers, 2);
	EXPECT_TRUE(spaces.Proposed == Parameter::PolynomialSpace);
}

// Refused (no reading): DoFs of no configuration, not divisible by the components, a parameter set contradicting them,
// CellDegree outside k - 1 to k + 1, sizes or components not positive
TEST(DiscretizationInference, InvalidInputs)
{
	fhhos4::PCoarsening minus2 = fhhos4::PCoarsening::Minus2;
	EXPECT_EQ(InferDiscretization(7, 3, 1, Discretization(), minus2).Solvers, 0);
	EXPECT_EQ(InferDiscretization(6, 3, 2, Discretization(), minus2).Solvers, 0); // 2D P k=2 l=1 on 2 components has 6 and 6
	EXPECT_EQ(InferDiscretization(3, 2, 1, Discretization(), minus2).Solvers, 1); // 2D P k=1 l=1
	EXPECT_EQ(InferDiscretization(3, 2, 1, Discretization{ 3 }, minus2).Solvers, 0);
	EXPECT_EQ(InferDiscretization(3, 2, 1, Discretization{ -1, Space::Q }, minus2).Solvers, 0);
	EXPECT_EQ(InferDiscretization(3, 2, 1, Discretization{ -1, Space::Inferred, 2 }, minus2).Solvers, 0);
	EXPECT_EQ(InferDiscretization(10, 2, 1, Discretization(), minus2).Solvers, 0); // P^3 in 2D, faces of degree 1
	EXPECT_EQ(InferDiscretization(10, 2, 1, Discretization{ -1, Space::Inferred, 1, 3 }, minus2).Solvers, 0);
	for (auto [cellBlockSize, faceBlockSize, nc] : { std::tuple{ 0, 1, 1 }, { 1, 0, 1 }, { 1, 1, 0 }, { -3, 2, 1 }, { 3, 2, -1 } })
		EXPECT_EQ(InferDiscretization(cellBlockSize, faceBlockSize, nc, Discretization(), minus2).Solvers, 0) << cellBlockSize << " " << faceBlockSize << " " << nc;
}

// The messages, exact: an ambiguity (with the members' names, and with the options of the PETSc preconditioner; the
// space proposed), no reading, and what a resolved inference prints (the parameters that don't change the solver: not
// needed)
TEST(DiscretizationInference, Messages)
{
	fhhos4::PCoarsening minus2 = fhhos4::PCoarsening::Minus2;
	Inference ambiguous = InferDiscretization(64, 9, 1, Discretization(), minus2);
	EXPECT_EQ(Message(ambiguous, MemberNames),
	          "fhhos4: 64 DoFs per cell and 9 per face (1 component) fit discretizations that need different solvers:\n"
	          "  Dimension 3, Q, face degree 2, cell degree 3\n"
	          "  Dimension 2, Q, face degree 8, cell degree 7\n"
	          "Set Dimension (2 or 3).");
	ParameterNames options = { "-pc_fhhos4_dimension", "-pc_fhhos4_polynomial_space", "-pc_fhhos4_face_degree", "-pc_fhhos4_cell_degree", "-pc_fhhos4_vector_components", "p", "q" };
	std::string message = Message(ambiguous, options);
	EXPECT_EQ(message.substr(message.rfind('\n') + 1), "Set -pc_fhhos4_dimension (2 or 3).");
	EXPECT_EQ(Message(InferDiscretization(36, 7, 1, Discretization(), minus2), options),
	          "fhhos4: 36 DoFs per cell and 7 per face (1 component) fit discretizations that need different solvers:\n"
	          "  Dimension 2, P, face degree 6, cell degree 7\n"
	          "  Dimension 2, Q, face degree 6, cell degree 5\n"
	          "Set -pc_fhhos4_polynomial_space (p or q).");
	EXPECT_EQ(Message(InferDiscretization(7, 3, 1, Discretization{ 3 }, minus2), MemberNames),
	          "fhhos4: 7 DoFs per cell and 3 per face (1 component), with Dimension 3, fit no discretization of the library: P^l or Q^l in the cells and P^k or Q^k on the faces, "
	          "in 2D or 3D, with l = k - 1, k or k + 1, the same in every cell and face. These DoFs are read from the stored entries of A_TT and A_TF: the zeros inside the blocks of "
	          "each cell must be stored. Check VectorComponents, or set Dimension, PolynomialSpace, FaceDegree and CellDegree.");
	EXPECT_EQ(Description(InferDiscretization(1, 1, 1, Discretization(), minus2), MemberNames),
	          "fhhos4: inferred from the blocks, 1 DoF per cell and 1 per face (1 component): Dimension 2 (not needed), PolynomialSpace P (not needed), FaceDegree 0, CellDegree 0");
	EXPECT_EQ(Description(InferDiscretization(64, 9, 1, Discretization{ 3 }, minus2), options),
	          "fhhos4: inferred from the blocks, 64 DoFs per cell and 9 per face (1 component): -pc_fhhos4_polynomial_space q, -pc_fhhos4_face_degree 2, -pc_fhhos4_cell_degree 3");
}

// The analyses of the processes combined: one process; processes without cells (no opinion), all but one; the GCDs of the
// processes (each a multiple of the DoFs per face) combined; refused: processes disagreeing on the cells, no face between
// two cells anywhere (also no cell at all), the message naming the parameters to set
TEST(DiscretizationInference, CombineAnalyses)
{
	using Sizes = std::pair<int, int>;
	EXPECT_EQ(CombineAnalyses({ { 6, 3, 1 } }, MemberNames), Sizes(6, 3));
	EXPECT_EQ(CombineAnalyses({ {}, { 6, 3, 1 }, { 6, 6, 0 } }, MemberNames), Sizes(6, 3));
	EXPECT_EQ(CombineAnalyses({ {}, {}, { 6, 6, 1 }, {} }, MemberNames), Sizes(6, 6));
	EXPECT_EQ(CombineAnalyses({ { 4, 6, 1 }, { 4, 9, 0 }, {} }, MemberNames), Sizes(4, 3));
	EXPECT_THROW(CombineAnalyses({ { 6, 3, 1 }, {}, { 3, 3, 1 } }, MemberNames), fhhos4::Error);
	for (std::vector<BlockAnalysis> processes : { std::vector<BlockAnalysis>{ { 6, 6, 0 } }, { { 6, 6, 0 }, { 6, 3, 0 } }, { {}, {} } })
	{
		try
		{
			CombineAnalyses(processes, MemberNames);
			ADD_FAILURE() << "no face between two cells accepted";
		}
		catch (const fhhos4::Error& e)
		{
			std::string expected = std::string("fhhos4: the blocks don't give the DoFs per face: ") + (processes.size() > 1 ? "no process has a face between two of its cells" : "no face is between two cells")
			                       + " (a column of A_TF in the rows of two cells). Set Dimension, PolynomialSpace, FaceDegree and CellDegree.";
			EXPECT_EQ(std::string(e.what()), expected);
		}
	}
}

namespace
{
	using Matrix = fhhos4::SparseMatrix;

	Matrix FromTriplets(int rows, int cols, const std::vector<Eigen::Triplet<double>>& triplets)
	{
		Matrix M(rows, cols);
		M.setFromTriplets(triplets.begin(), triplets.end());
		M.makeCompressed();
		return M;
	}

	// The blocks of a strip of cells: the cell c between the faces c and c + 1 (numbered consecutively, those of the
	// ends too), its block of A_TT dense (lower: lower-only), the values 1 on the diagonal and offDiagonal elsewhere
	// (stored, even if 0), A_TF dense
	std::pair<Matrix, Matrix> Strip(int cells, int cellDoFs, int faceDoFs, bool lower = false, double offDiagonal = 1)
	{
		std::vector<Eigen::Triplet<double>> TT, TF;
		for (int c = 0; c < cells; c++)
		{
			for (int i = 0; i < cellDoFs; i++)
			{
				for (int j = 0; j < (lower ? i + 1 : cellDoFs); j++)
					TT.emplace_back(c * cellDoFs + i, c * cellDoFs + j, i == j ? 1 : offDiagonal);
				for (int f : { c, c + 1 })
					for (int k = 0; k < faceDoFs; k++)
						TF.emplace_back(c * cellDoFs + i, f * faceDoFs + k, 1);
			}
		}
		return { FromTriplets(cells * cellDoFs, cells * cellDoFs, TT), FromTriplets(cells * cellDoFs, (cells + 1) * faceDoFs, TF) };
	}

	std::tuple<std::int64_t, std::int64_t, std::int64_t> Tuple(const BlockAnalysis& a)
	{
		return { a.CellBlockSize, a.FaceBoundaryGCD, a.InteriorFaceSeen };
	}
}

// The analysis of tiny blocks: a strip of cells whose faces are numbered consecutively (the GCD of the boundaries: the
// DoFs per face, once two cells share a face), A_TT full or lower-only, with exact zeros stored in its blocks (the
// pattern counts, not the values); one cell (its two faces read as one, no face between two cells); no cell; a
// vector-valued unknown whose components are uncoupled, the zeros between them not stored in A_TT (full: the row of
// phi_0 e_x reaches the last e_x of the cell); refused: cells of different sizes, and the same unknown in a lower-only
// A_TT (the first block then one DoF: the zeros must be stored)
TEST(DiscretizationInference, AnalyzeBlocks)
{
	using Analysis = std::tuple<std::int64_t, std::int64_t, std::int64_t>;
	for (bool lower : { false, true })
	{
		for (double offDiagonal : { 1.0, 0.0 })
		{
			auto [A_TT, A_TF] = Strip(4, 3, 2, lower, offDiagonal);
			EXPECT_EQ(Tuple(AnalyzeBlocks(A_TT, A_TF)), Analysis(3, 2, 1)) << "lower " << lower << ", off-diagonal " << offDiagonal;
		}
	}
	auto [oneTT, oneTF] = Strip(1, 3, 2);
	EXPECT_EQ(Tuple(AnalyzeBlocks(oneTT, oneTF)), Analysis(3, 4, 0));
	Matrix none(0, 0), noneTF(0, 5);
	EXPECT_EQ(Tuple(AnalyzeBlocks(none, noneTF)), Analysis(0, 0, 0));

	// 2 components of 3 scalar functions (P^1 in 2D), interleaved: the DoF 2 i + c is phi_i e_c, A_TT couples a
	// component to itself only; A_TF dense (faces of 2 x 2 DoFs)
	for (bool lower : { false, true })
	{
		std::vector<Eigen::Triplet<double>> TT;
		for (int c = 0; c < 3; c++)
			for (int i = 0; i < 6; i++)
				for (int j = i % 2; j < (lower ? i + 1 : 6); j += 2)
					TT.emplace_back(6 * c + i, 6 * c + j, i == j ? 2 : -0.5);
		Matrix A_TT = FromTriplets(18, 18, TT);
		Matrix A_TF = Strip(3, 6, 4).second;
		if (!lower)
			EXPECT_EQ(Tuple(AnalyzeBlocks(A_TT, A_TF)), Analysis(6, 4, 1));
		else
			EXPECT_THROW(AnalyzeBlocks(A_TT, A_TF), fhhos4::Error);
	}

	// Cells of 2 and 1 DoFs, of 1 and 2: refused; A_TT not of the rows of A_TF
	Matrix twoThenOne = FromTriplets(3, 3, { { 0, 0, 1 }, { 0, 1, 1 }, { 1, 0, 1 }, { 1, 1, 1 }, { 2, 2, 1 } });
	Matrix oneThenTwo = FromTriplets(3, 3, { { 0, 0, 1 }, { 1, 1, 1 }, { 1, 2, 1 }, { 2, 1, 1 }, { 2, 2, 1 } });
	Matrix A_TF = FromTriplets(3, 4, { { 0, 0, 1 }, { 1, 1, 1 }, { 2, 2, 1 } });
	EXPECT_THROW(AnalyzeBlocks(twoThenOne, A_TF), fhhos4::Error);
	EXPECT_THROW(AnalyzeBlocks(oneThenTwo, A_TF), fhhos4::Error);
	EXPECT_THROW(AnalyzeBlocks(twoThenOne, Matrix(2, 4)), fhhos4::Error);
}
