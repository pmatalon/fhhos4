// The Python module fhhos4 (nanobind): the API of the library (library/include/fhhos4.h), with the same names, on SciPy's
// sparse matrices (CSR, or any format that SciPy converts to it) and NumPy's vectors (float64). See README_library.md.
#include <memory>
#include <nanobind/nanobind.h>
#include <nanobind/eigen/dense.h>
#include <nanobind/eigen/sparse.h>
#include <nanobind/stl/string.h>
#include "fhhos4"

namespace nb = nanobind;
using namespace nb::literals;

using SparseMatrix = fhhos4::SparseMatrix;
using ConstVector = Eigen::Ref<const Eigen::VectorXd>;
using Vector = Eigen::Ref<Eigen::VectorXd>;

// fhhos4::Solver, which also owns the condensed matrix A of the setups from A: the library keeps a reference to it until
// the next setup, and the conversion of the SciPy matrix is a temporary. The new A replaces the previous one only after
// a successful setup: a setup that fails (on its checks) keeps the previous one, which still reads its A.
struct PySolver : fhhos4::Solver
{
	std::unique_ptr<SparseMatrix> A;

	template <typename Setup>
	void SetupCondensed(const SparseMatrix& A, Setup setup)
	{
		auto ownA = std::make_unique<SparseMatrix>(A);
		setup(*ownA);
		this->A = std::move(ownA); // the object stays where the library points to
	}
};

NB_MODULE(fhhos4, m)
{
	m.doc() = "fhhos4: the algebraic multigrid U-AMG for hybrid discretizations (HHO, HDG...). See README_library.md.";

	nb::exception<fhhos4::Error>(m, "Error", PyExc_RuntimeError);
	m.def("Version", &fhhos4::Version, "The version of the library, with the commit when built from a clone");

	nb::enum_<fhhos4::Krylov>(m, "Krylov")
		.value("Default", fhhos4::Krylov::Default)
		.value("FCG", fhhos4::Krylov::FCG)
		.value("GCR", fhhos4::Krylov::GCR)
		.value("BiCGSTAB", fhhos4::Krylov::BiCGSTAB)
		.value("None", fhhos4::Krylov::None);

	nb::enum_<fhhos4::Cycle>(m, "Cycle")
		.value("V", fhhos4::Cycle::V)
		.value("W", fhhos4::Cycle::W)
		.value("K", fhhos4::Cycle::K);

	nb::enum_<fhhos4::PCoarsening>(m, "PCoarsening")
		.value("Minus2", fhhos4::PCoarsening::Minus2)
		.value("Minus1", fhhos4::PCoarsening::Minus1)
		.value("DivideBy2", fhhos4::PCoarsening::DivideBy2)
		.value("ToZero", fhhos4::PCoarsening::ToZero);

	nb::enum_<fhhos4::Smoother>(m, "Smoother")
		.value("Default", fhhos4::Smoother::Default)
		.value("BlockGaussSeidel", fhhos4::Smoother::BlockGaussSeidel)
		.value("BackwardBlockGaussSeidel", fhhos4::Smoother::BackwardBlockGaussSeidel)
		.value("SymmetricBlockGaussSeidel", fhhos4::Smoother::SymmetricBlockGaussSeidel)
		.value("HybridBlockGaussSeidel", fhhos4::Smoother::HybridBlockGaussSeidel)
		.value("HybridBackwardBlockGaussSeidel", fhhos4::Smoother::HybridBackwardBlockGaussSeidel)
		.value("DampedBlockJacobi", fhhos4::Smoother::DampedBlockJacobi)
		.value("BlockChebyshevFirstKind", fhhos4::Smoother::BlockChebyshevFirstKind)
		.value("BlockChebyshevFourthKind", fhhos4::Smoother::BlockChebyshevFourthKind);

	nb::enum_<fhhos4::CoarseSolver>(m, "CoarseSolver")
		.value("Default", fhhos4::CoarseSolver::Default)
		.value("Cholesky", fhhos4::CoarseSolver::Cholesky)
		.value("LU", fhhos4::CoarseSolver::LU);

	nb::class_<fhhos4::Result>(m, "Result")
		.def_ro("Iterations", &fhhos4::Result::Iterations)
		.def_ro("RelativeResidual", &fhhos4::Result::RelativeResidual)
		.def_ro("Converged", &fhhos4::Result::Converged)
		.def("__repr__", [](const fhhos4::Result& r)
		{
			return "Result(Iterations=" + std::to_string(r.Iterations) + ", RelativeResidual=" + std::to_string(r.RelativeResidual)
				+ ", Converged=" + (r.Converged ? "True" : "False") + ")";
		});

	// The setups and solves release the GIL: they only read the converted matrices and the NumPy vectors
	auto noGIL = nb::call_guard<nb::gil_scoped_release>();

	nb::class_<PySolver>(m, "Solver", "fhhos4::Solver: set the parameters, then call a setup function, then Solve()")
		.def(nb::init<>())
		// Discretization
		.def_rw("Dimension", &fhhos4::Solver::Dimension)
		.def_rw("FaceDegree", &fhhos4::Solver::FaceDegree)
		.def_rw("CellDegree", &fhhos4::Solver::CellDegree)
		// Solver
		.def_rw("Krylov", &fhhos4::Solver::Krylov)
		.def_rw("Tolerance", &fhhos4::Solver::Tolerance)
		.def_rw("MaxIterations", &fhhos4::Solver::MaxIterations)
		// Multigrid
		.def_rw("Cycle", &fhhos4::Solver::Cycle)
		.def_rw("PCoarsening", &fhhos4::Solver::PCoarsening)
		.def_rw("PreSmoothingIterations", &fhhos4::Solver::PreSmoothingIterations)
		.def_rw("PostSmoothingIterations", &fhhos4::Solver::PostSmoothingIterations)
		.def_rw("PreSmoother", &fhhos4::Solver::PreSmoother)
		.def_rw("PostSmoother", &fhhos4::Solver::PostSmoother)
		.def_rw("CoarseningFactor", &fhhos4::Solver::CoarseningFactor)
		.def_rw("CoarseMatrixMaxSize", &fhhos4::Solver::CoarseMatrixMaxSize)
		.def_rw("CoarseSolver", &fhhos4::Solver::CoarseSolver)
		// Execution
		.def_rw("Threads", &fhhos4::Solver::Threads)
		.def_rw("Verbosity", &fhhos4::Solver::Verbosity)

		// Uncondensed system
		.def("SetupSymmetric", [](PySolver& s, const SparseMatrix& A_TT, const SparseMatrix& A_TF, const SparseMatrix& A_FF, const ConstVector& cellInterpOfOne, const ConstVector& faceInterpOfOne)
			{ s.SetupSymmetric(A_TT, A_TF, A_FF, cellInterpOfOne, faceInterpOfOne); },
			"A_TT"_a, "A_TF"_a, "A_FF"_a, "cellInterpOfOne"_a, "faceInterpOfOne"_a, noGIL)
		.def("SetupNonsymmetric", [](PySolver& s, const SparseMatrix& A_TT, const SparseMatrix& A_TF, const SparseMatrix& A_FT, const SparseMatrix& A_FF, const ConstVector& cellInterpOfOne, const ConstVector& faceInterpOfOne)
			{ s.SetupNonsymmetric(A_TT, A_TF, A_FT, A_FF, cellInterpOfOne, faceInterpOfOne); },
			"A_TT"_a, "A_TF"_a, "A_FT"_a, "A_FF"_a, "cellInterpOfOne"_a, "faceInterpOfOne"_a, noGIL)
		// x_T, x_F: NumPy vectors of float64, written in place (x_F: the initial guess on input)
		.def("Solve", [](PySolver& s, const ConstVector& b_T, const ConstVector& b_F, Vector x_T, Vector x_F)
			{ return s.Solve(b_T, b_F, x_T, x_F); },
			"b_T"_a, "b_F"_a, "x_T"_a.noconvert(), "x_F"_a.noconvert(), noGIL)

		// Condensed system
		.def("SetupSymmetricCondensed", [](PySolver& s, const SparseMatrix& A, const SparseMatrix& A_TT, const SparseMatrix& A_TF, const SparseMatrix& A_FF, const ConstVector& cellInterpOfOne, const ConstVector& faceInterpOfOne)
			{ s.SetupCondensed(A, [&](const SparseMatrix& ownA) { s.SetupSymmetricCondensed(ownA, A_TT, A_TF, A_FF, cellInterpOfOne, faceInterpOfOne); }); },
			"A"_a, "A_TT"_a, "A_TF"_a, "A_FF"_a, "cellInterpOfOne"_a, "faceInterpOfOne"_a, noGIL)
		.def("SetupSymmetricCondensed", [](PySolver& s, const SparseMatrix& A, const SparseMatrix& A_TT, const SparseMatrix& A_TF, const ConstVector& cellInterpOfOne, const ConstVector& faceInterpOfOne)
			{ s.SetupCondensed(A, [&](const SparseMatrix& ownA) { s.SetupSymmetricCondensed(ownA, A_TT, A_TF, cellInterpOfOne, faceInterpOfOne); }); },
			"A"_a, "A_TT"_a, "A_TF"_a, "cellInterpOfOne"_a, "faceInterpOfOne"_a, noGIL)
		.def("SetupNonsymmetricCondensed", [](PySolver& s, const SparseMatrix& A, const SparseMatrix& A_TT, const SparseMatrix& A_TF, const SparseMatrix& A_FT, const SparseMatrix& A_FF, const ConstVector& cellInterpOfOne, const ConstVector& faceInterpOfOne)
			{ s.SetupCondensed(A, [&](const SparseMatrix& ownA) { s.SetupNonsymmetricCondensed(ownA, A_TT, A_TF, A_FT, A_FF, cellInterpOfOne, faceInterpOfOne); }); },
			"A"_a, "A_TT"_a, "A_TF"_a, "A_FT"_a, "A_FF"_a, "cellInterpOfOne"_a, "faceInterpOfOne"_a, noGIL)
		.def("SetupNonsymmetricCondensed", [](PySolver& s, const SparseMatrix& A, const SparseMatrix& A_TT, const SparseMatrix& A_TF, const SparseMatrix& A_FT, const ConstVector& cellInterpOfOne, const ConstVector& faceInterpOfOne)
			{ s.SetupCondensed(A, [&](const SparseMatrix& ownA) { s.SetupNonsymmetricCondensed(ownA, A_TT, A_TF, A_FT, cellInterpOfOne, faceInterpOfOne); }); },
			"A"_a, "A_TT"_a, "A_TF"_a, "A_FT"_a, "cellInterpOfOne"_a, "faceInterpOfOne"_a, noGIL)
		// x: a NumPy vector of float64, written in place (the initial guess on input)
		.def("Solve", [](PySolver& s, const ConstVector& b, Vector x) { return s.Solve(b, x); },
			"b"_a, "x"_a.noconvert(), noGIL)
		// Solves A x = b from x = 0, and returns x
		.def("Solve", [](PySolver& s, const ConstVector& b) { return s.Solve(b); }, "b"_a, noGIL)

		// A copy of the condensed matrix (scipy.sparse.csr_matrix)
		.def("CondensedMatrix", [](const PySolver& s) { return SparseMatrix(s.CondensedMatrix()); })
		.def("IsSetUp", &fhhos4::Solver::IsSetUp)
		.def("NumberOfLevels", &fhhos4::Solver::NumberOfLevels);
}
