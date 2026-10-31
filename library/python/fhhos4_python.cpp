// The Python module fhhos4 (nanobind): the API of the library (library/include/fhhos4.h), with the same names, on NumPy
// arrays read in place, as the library reads its inputs: never copied nor converted (an array of another dtype, or not
// contiguous, is refused). The matrices (no Eigen in Python, the arrays of the library's implementation): a CSRView of
// NumPy arrays, a CSRMatrix (the result of the helpers that copy, as ToCSR()), or any object with the attributes
// indptr, indices, data, shape and format "csr" (SciPy's CSR matrices: SciPy itself is not needed); the vectors: NumPy
// vectors of float64. See README_library.md.
// Compiled with MPI (FHHOS4_MPI): Solver.Comm, from an mpi4py communicator, and the ghost faces (and cells) of the setups.
#include <cstdint>
#include <string>
#include <vector>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/eigen/dense.h>
#include <nanobind/stl/string.h>
#include "fhhos4"

namespace nb = nanobind;
using namespace nb::literals;

using ConstVector = Eigen::Ref<const Eigen::VectorXd>;
using Vector = Eigen::Ref<Eigen::VectorXd>;

namespace
{
	// The data of a 1-D C-contiguous array on the CPU, read in place. is64: an index array of 64 bits (int64, or uint64
	// read as int64: e.g. NGSolve's row pointers), 32 bits otherwise (int32).
	const void* IndexArray(nb::handle array, const char* name, bool& is64, size_t& size)
	{
		using Array = nb::ndarray<nb::ndim<1>, nb::c_contig, nb::device::cpu, nb::ro>;
		Array a;
		if (!nb::try_cast(array, a, false))
			throw nb::type_error((std::string("fhhos4: ") + name + " must be a contiguous 1-D NumPy array, read in place (not copied)").c_str());
		if (a.dtype() == nb::dtype<std::int64_t>() || a.dtype() == nb::dtype<std::uint64_t>())
			is64 = true;
		else if (a.dtype() == nb::dtype<std::int32_t>())
			is64 = false;
		else
			throw nb::type_error((std::string("fhhos4: ") + name + " must be an array of int32 or int64 (not converted)").c_str());
		size = a.shape(0);
		return a.data();
	}

	const double* ValueArray(nb::handle array, const char* name, size_t& size)
	{
		using Array = nb::ndarray<double, nb::ndim<1>, nb::c_contig, nb::device::cpu, nb::ro>;
		Array a;
		if (!nb::try_cast(array, a, false))
			throw nb::type_error((std::string("fhhos4: ") + name + " must be a contiguous 1-D NumPy array of float64, read in place (not converted)").c_str());
		size = a.shape(0);
		return a.data();
	}

	std::int64_t Last(const void* array, bool is64, size_t size)
	{
		return is64 ? ((const std::int64_t*)array)[size - 1] : ((const std::int32_t*)array)[size - 1];
	}

	// The arrays of a matrix in the CSR format, read in place, of any index types (32 or 64 bits, independently), and the
	// Python objects that own them: alive while the view is
	struct PyCSRView
	{
		std::int64_t Rows = 0, Cols = 0;
		const void* RowPointers = nullptr;
		const void* ColumnIndices = nullptr;
		const double* Values = nullptr;
		bool RowPointers64 = false, Columns64 = false;
		std::vector<nb::object> Owners;

		PyCSRView(std::int64_t rows, std::int64_t cols, nb::handle rowPointers, nb::handle columnIndices, nb::handle values)
			: Rows(rows), Cols(cols), Owners{ nb::borrow(rowPointers), nb::borrow(columnIndices), nb::borrow(values) }
		{
			size_t nRowPointers, nColumns, nValues;
			RowPointers = IndexArray(rowPointers, "the row pointers", RowPointers64, nRowPointers);
			ColumnIndices = IndexArray(columnIndices, "the column indices", Columns64, nColumns);
			Values = ValueArray(values, "the values", nValues);
			if (rows < 0 || cols < 0 || nRowPointers != (size_t)rows + 1)
				throw nb::value_error(("fhhos4: " + std::to_string(rows) + " rows: " + std::to_string(rows + 1) + " row pointers expected, got " + std::to_string(nRowPointers) + ".").c_str());
			std::int64_t nnz = Last(RowPointers, RowPointers64, nRowPointers);
			if (nnz < 0 || (size_t)nnz > nColumns || (size_t)nnz > nValues)
				throw nb::value_error(("fhhos4: " + std::to_string(nnz) + " non-zeros (the last row pointer), but " + std::to_string(nColumns) + " column indices and " + std::to_string(nValues) + " values.").c_str());
		}

		// The view of a CSRMatrix (its Python object: owner), as its arrays
		PyCSRView(const fhhos4::detail::CSRMatrix& M, nb::handle owner)
			: Owners{ nb::borrow(owner) }
		{
			fhhos4::detail::CSRView v = M.View();
			Rows = v.Rows;
			Cols = v.Cols;
			RowPointers = v.RowPointers;
			ColumnIndices = v.ColumnIndices;
			Values = v.Values;
			RowPointers64 = Columns64 = sizeof(fhhos4::SparseIndex) == 8;
		}

		template <typename RowPointer, typename Column>
		fhhos4::detail::BasicCSRView<RowPointer, Column> View() const
		{
			return { Rows, Cols, (const RowPointer*)RowPointers, (const Column*)ColumnIndices, Values };
		}

		// f(view), with the view of the index types of the arrays
		template <typename F>
		auto Visit(F f) const
		{
			if (RowPointers64)
				return Columns64 ? f(View<std::int64_t, std::int64_t>()) : f(View<std::int64_t, std::int32_t>());
			return Columns64 ? f(View<std::int32_t, std::int64_t>()) : f(View<std::int32_t, std::int32_t>());
		}

		// The view of the setups: the index type of the library
		fhhos4::detail::CSRView Native(const char* name) const
		{
			constexpr bool library64 = sizeof(fhhos4::SparseIndex) == 8;
			if (RowPointers64 != library64 || Columns64 != library64)
				throw nb::type_error((std::string("fhhos4: the indices of ") + name + " are not of the index type of the library (" + (library64 ? "int64" : "int32") + "): convert the matrix with fhhos4.ToCSR() (a copy)").c_str());
			return View<fhhos4::SparseIndex, fhhos4::SparseIndex>();
		}
	};

	std::string Format(nb::handle M)
	{
		return nb::hasattr(M, "format") ? nb::cast<std::string>(M.attr("format")) : std::string();
	}

	// The arrays of an object with the attributes indptr, indices, data and shape (e.g. SciPy's sparse matrices)
	PyCSRView FromAttributes(nb::handle M, bool transposed = false)
	{
		nb::object shape = M.attr("shape");
		std::int64_t rows = nb::cast<std::int64_t>(shape[transposed ? 1 : 0]), cols = nb::cast<std::int64_t>(shape[transposed ? 0 : 1]);
		return PyCSRView(rows, cols, M.attr("indptr"), M.attr("indices"), M.attr("data"));
	}

	bool HasArrays(nb::handle M)
	{
		return nb::hasattr(M, "indptr") && nb::hasattr(M, "indices") && nb::hasattr(M, "data") && nb::hasattr(M, "shape");
	}

	// The matrix M, read in place: a CSRView, a CSRMatrix, or an object with the CSR arrays (format "csr")
	PyCSRView ToView(nb::handle M, const char* name)
	{
		if (nb::isinstance<PyCSRView>(M))
			return nb::cast<const PyCSRView&>(M);
		if (nb::isinstance<fhhos4::detail::CSRMatrix>(M))
			return PyCSRView(nb::cast<const fhhos4::detail::CSRMatrix&>(M), M);
		if (HasArrays(M))
		{
			std::string format = Format(M);
			if (!format.empty() && format != "csr")
				throw nb::type_error((std::string("fhhos4: ") + name + " is a matrix in the format " + format + ", not CSR: convert it with fhhos4.ToCSR() (a copy), or, if it is a symmetric CSC matrix, view it with fhhos4.TransposedView() (no copy)").c_str());
			return FromAttributes(M);
		}
		throw nb::type_error((std::string("fhhos4: ") + name + " must be a fhhos4.CSRView, a fhhos4.CSRMatrix, or a CSR matrix (e.g. scipy.sparse.csr_array)").c_str());
	}

	// The DoF lists (fhhos4::DoFLists) of NumPy arrays of the same index type
	struct PyDoFLists
	{
		std::int64_t Count = 0;
		const void* Offsets = nullptr;
		const void* DoFs = nullptr;
		bool Is64 = false;
		std::vector<nb::object> Owners;

		PyDoFLists(std::int64_t count, nb::handle offsets, nb::handle dofs)
			: Count(count), Owners{ nb::borrow(offsets), nb::borrow(dofs) }
		{
			bool dofs64;
			size_t nOffsets, nDoFs;
			Offsets = IndexArray(offsets, "the offsets", Is64, nOffsets);
			DoFs = IndexArray(dofs, "the DoFs", dofs64, nDoFs);
			if (dofs64 != Is64)
				throw nb::type_error("fhhos4: the offsets and DoFs must have the same index type (int32 or int64)");
			if (count < 0 || nOffsets != (size_t)count + 1)
				throw nb::value_error(("fhhos4: " + std::to_string(count) + " lists: " + std::to_string(count + 1) + " offsets expected, got " + std::to_string(nOffsets) + ".").c_str());
			if (Last(Offsets, Is64, nOffsets) > (std::int64_t)nDoFs)
				throw nb::value_error("fhhos4: the last offset is beyond the DoFs.");
		}

		template <typename Index>
		fhhos4::DoFLists<Index> Lists() const
		{
			return { Count, (const Index*)Offsets, (const Index*)DoFs };
		}
	};

	// f(cells, faces), with the lists of their index type
	template <typename F>
	auto VisitLists(const PyDoFLists& cells, const PyDoFLists& faces, F f)
	{
		if (cells.Is64 != faces.Is64)
			throw nb::type_error("fhhos4: the DoF lists of the cells and faces must have the same index type");
		if (cells.Is64)
			return f(cells.Lists<std::int64_t>(), faces.Lists<std::int64_t>());
		return f(cells.Lists<std::int32_t>(), faces.Lists<std::int32_t>());
	}

	// fhhos4::Solver, which also keeps the views of the matrices that the solves read (A_FT of SetupNonSym()): their
	// arrays stay alive until the next setup. Replaced only after a successful
	// setup: a setup that fails (on its checks) keeps the previous one, which still reads them.
	struct PySolver : fhhos4::Solver
	{
		std::vector<PyCSRView> Kept;
		nb::object CommObject = nb::none(); // the Python object of Comm (MPI)

		// The setups of the library on the views of the arrays (fhhos4::Solver::Setup())
		void SetupViews(bool symmetric, const fhhos4::detail::CSRView& A_TT, const fhhos4::detail::CSRView& A_TF, const fhhos4::detail::CSRView* A_FT,
		                const fhhos4::detail::CSRView& A_FF, const ConstVector& cellInterpOfOne, const ConstVector& faceInterpOfOne,
		                const fhhos4::GhostNumbers& faceGhosts = fhhos4::GhostNumbers(), const fhhos4::GhostNumbers& cellGhosts = fhhos4::GhostNumbers())
		{
			Setup(symmetric, A_TT, A_TF, A_FT, A_FF, In(cellInterpOfOne), In(faceInterpOfOne), faceGhosts, cellGhosts);
		}

		fhhos4::detail::CSRView CondensedView() const
		{
			return CondensedMatrixView();
		}

		// setup(): the setup of the library, on the views, without the GIL. kept: the views the solves read.
		template <typename Setup>
		void Run(Setup setup, std::vector<PyCSRView> kept)
		{
			{
				nb::gil_scoped_release release;
				setup();
			}
			Kept = std::move(kept);
		}
	};

	// The global numbers of the ghost faces or cells (the setups on several processes): an int64 NumPy array, read in
	// place (only by the setup), or None (no ghosts)
	fhhos4::GhostNumbers ToGhostNumbers(nb::handle ghosts, const char* name)
	{
		if (ghosts.is_none())
			return fhhos4::GhostNumbers();
		bool is64;
		size_t size;
		const void* data = IndexArray(ghosts, name, is64, size);
		if (!is64)
			throw nb::type_error((std::string("fhhos4: ") + name + " must be an array of int64 (global numbers), not converted").c_str());
		return fhhos4::GhostNumbers(size, (const std::int64_t*)data);
	}

#if FHHOS4_MPI
	// The communicator of an mpi4py communicator (MPI.COMM_WORLD...), through its Fortran handle (py2f()), or of a
	// Fortran handle. MPI must be initialized (by mpi4py).
	MPI_Comm ToComm(nb::handle comm)
	{
		int initialized = 0;
		MPI_Initialized(&initialized);
		if (!initialized)
			throw fhhos4::Error("fhhos4: MPI is not initialized: import mpi4py (from mpi4py import MPI) before setting Solver.Comm.");
		nb::object handle = nb::hasattr(comm, "py2f") ? comm.attr("py2f")() : nb::borrow(comm);
		return MPI_Comm_f2c(nb::cast<MPI_Fint>(handle));
	}
#endif

	// The arrays of a CSRMatrix, as read-only NumPy arrays that keep it alive
	template <typename T>
	nb::ndarray<nb::numpy, const T, nb::ndim<1>> ArrayOf(const fhhos4::detail::CSRMatrix& M, const T* data, std::int64_t size)
	{
		size_t shape[1] = { (size_t)size };
		return nb::ndarray<nb::numpy, const T, nb::ndim<1>>(data, 1, shape, nb::find(M));
	}
}

NB_MODULE(fhhos4, m)
{
	m.doc() = "fhhos4: the algebraic multigrid U-AMG for hybrid discretizations (HHO, HDG...). See README_library.md.";

	nb::exception<fhhos4::Error>(m, "Error", PyExc_RuntimeError);
	m.def("Version", &fhhos4::Version, "The version of the library, with the commit when built from a clone");
	m.attr("WithMPI") = (bool)FHHOS4_MPI; // compiled with MPI: Solver.Comm (not set: one process)

	nb::enum_<fhhos4::Krylov>(m, "Krylov")
		.value("Default", fhhos4::Krylov::Default)
		.value("FCG", fhhos4::Krylov::FCG)
		.value("GCR", fhhos4::Krylov::GCR)
		.value("None", fhhos4::Krylov::None);

	nb::enum_<fhhos4::Cycle>(m, "Cycle")
		.value("V", fhhos4::Cycle::V)
		.value("W", fhhos4::Cycle::W)
		.value("K", fhhos4::Cycle::K);

	nb::enum_<fhhos4::PolynomialSpace>(m, "PolynomialSpace")
		.value("P", fhhos4::PolynomialSpace::P)
		.value("Q", fhhos4::PolynomialSpace::Q);
	nb::enum_<fhhos4::PCoarsening>(m, "PCoarsening")
		.value("Minus2", fhhos4::PCoarsening::Minus2)
		.value("Minus1", fhhos4::PCoarsening::Minus1)
		.value("DivideBy2", fhhos4::PCoarsening::DivideBy2)
		.value("ToZero", fhhos4::PCoarsening::ToZero);
	nb::enum_<fhhos4::Prolongation>(m, "Prolongation")
		.value("Reconstruction", fhhos4::Prolongation::Reconstruction)
		.value("HarmonicExtension", fhhos4::Prolongation::HarmonicExtension);

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
		.value("LU", fhhos4::CoarseSolver::LU)
		.value("MUMPS", fhhos4::CoarseSolver::MUMPS);

	nb::class_<fhhos4::Result>(m, "Result")
		.def_ro("Iterations", &fhhos4::Result::Iterations)
		.def_ro("RelativeResidual", &fhhos4::Result::RelativeResidual)
		.def_ro("Converged", &fhhos4::Result::Converged)
		.def("__repr__", [](const fhhos4::Result& r)
		{
			return "Result(Iterations=" + std::to_string(r.Iterations) + ", RelativeResidual=" + std::to_string(r.RelativeResidual)
				+ ", Converged=" + (r.Converged ? "True" : "False") + ")";
		});

	//----------------------------------------------- Matrices -----------------------------------------------//

	nb::class_<PyCSRView>(m, "CSRView", "A matrix in the CSR format: NumPy arrays read in place (int32 or int64 indices, float64 values)")
		.def(nb::init<std::int64_t, std::int64_t, nb::handle, nb::handle, nb::handle>(), "rows"_a, "cols"_a, "rowPointers"_a, "columnIndices"_a, "values"_a)
		.def_ro("Rows", &PyCSRView::Rows)
		.def_ro("Cols", &PyCSRView::Cols);

	nb::class_<fhhos4::detail::CSRMatrix>(m, "CSRMatrix", "A matrix in the CSR format of the library, owned by the library: the result of the helpers that copy")
		.def_prop_ro("Rows", [](const fhhos4::detail::CSRMatrix& M) { return M.View().Rows; })
		.def_prop_ro("Cols", [](const fhhos4::detail::CSRMatrix& M) { return M.View().Cols; })
		// Read-only NumPy arrays of the library's memory, which keep the matrix alive
		.def_prop_ro("RowPointers", [](const fhhos4::detail::CSRMatrix& M) { fhhos4::detail::CSRView v = M.View(); return ArrayOf(M, v.RowPointers, v.RowPointers ? v.Rows + 1 : 0); })
		.def_prop_ro("ColumnIndices", [](const fhhos4::detail::CSRMatrix& M) { fhhos4::detail::CSRView v = M.View(); return ArrayOf(M, v.ColumnIndices, v.NonZeros()); })
		.def_prop_ro("Values", [](const fhhos4::detail::CSRMatrix& M) { fhhos4::detail::CSRView v = M.View(); return ArrayOf(M, v.Values, v.NonZeros()); });

	nb::class_<PyDoFLists>(m, "DoFLists", "The DoFs of each cell (or face) in the numbering of a global matrix: NumPy arrays read in place")
		.def(nb::init<std::int64_t, nb::handle, nb::handle>(), "count"_a, "offsets"_a, "dofs"_a)
		.def_ro("Count", &PyDoFLists::Count);

	nb::class_<fhhos4::detail::HybridBlocks>(m, "HybridBlocks", "The blocks of a hybrid system, owned by the library")
		.def_prop_ro("A_TT", [](const fhhos4::detail::HybridBlocks& b) -> const fhhos4::detail::CSRMatrix& { return b.A_TT; }, nb::rv_policy::reference_internal)
		.def_prop_ro("A_TF", [](const fhhos4::detail::HybridBlocks& b) -> const fhhos4::detail::CSRMatrix& { return b.A_TF; }, nb::rv_policy::reference_internal)
		.def_prop_ro("A_FT", [](const fhhos4::detail::HybridBlocks& b) -> const fhhos4::detail::CSRMatrix& { return b.A_FT; }, nb::rv_policy::reference_internal)
		.def_prop_ro("A_FF", [](const fhhos4::detail::HybridBlocks& b) -> const fhhos4::detail::CSRMatrix& { return b.A_FF; }, nb::rv_policy::reference_internal);

	//------------------------------------------------ Helpers -----------------------------------------------//

	m.def("ExtractBlocksSym", [](nb::handle M, const PyDoFLists& cells, const PyDoFLists& faces)
		{
			PyCSRView view = ToView(M, "M");
			nb::gil_scoped_release release;
			return view.Visit([&](auto v) { return VisitLists(cells, faces, [&](auto c, auto f) { return fhhos4::detail::ExtractBlocksSym(v, c, f); }); });
		}, "M"_a, "cells"_a, "faces"_a, "A copy: the blocks of the uncondensed system, from its global matrix, for SetupSym()");

	m.def("ExtractBlocksNonSym", [](nb::handle M, const PyDoFLists& cells, const PyDoFLists& faces)
		{
			PyCSRView view = ToView(M, "M");
			nb::gil_scoped_release release;
			return view.Visit([&](auto v) { return VisitLists(cells, faces, [&](auto c, auto f) { return fhhos4::detail::ExtractBlocksNonSym(v, c, f); }); });
		}, "M"_a, "cells"_a, "faces"_a, "A copy: the four blocks of the uncondensed system, from its global matrix, for SetupNonSym()");

	m.def("Gather", [](const ConstVector& x, const PyDoFLists& cells, const PyDoFLists& faces, Vector x_T, Vector x_F)
		{
			VisitLists(cells, faces, [&](auto c, auto f) { fhhos4::Gather(x, c, f, x_T, x_F); });
		}, "x"_a.noconvert(), "cells"_a, "faces"_a, "x_T"_a.noconvert(), "x_F"_a.noconvert(), "A copy: the cells and faces DoFs of x, into x_T and x_F");

	m.def("Scatter", [](const ConstVector& x_T, const ConstVector& x_F, const PyDoFLists& cells, const PyDoFLists& faces, Vector x)
		{
			VisitLists(cells, faces, [&](auto c, auto f) { fhhos4::Scatter(x_T, x_F, c, f, x); });
		}, "x_T"_a.noconvert(), "x_F"_a.noconvert(), "cells"_a, "faces"_a, "x"_a.noconvert(), "A copy: x_T and x_F into x (its other DoFs unchanged)");

	m.def("ToCSR", [](nb::handle M)
		{
			if (HasArrays(M) && Format(M) == "csc") // the arrays of the transpose
			{
				PyCSRView transpose = FromAttributes(M, true);
				if (transpose.RowPointers64 != transpose.Columns64)
					throw nb::type_error("fhhos4.ToCSR(): the indices of a CSC matrix must have the same type");
				nb::gil_scoped_release release;
				if (transpose.RowPointers64)
					return fhhos4::detail::Transpose(transpose.View<std::int64_t, std::int64_t>());
				return fhhos4::detail::Transpose(transpose.View<std::int32_t, std::int32_t>());
			}
			PyCSRView view = ToView(M, "M");
			nb::gil_scoped_release release;
			return view.Visit([](auto v) { return fhhos4::detail::ToCSR(v); });
		}, "M"_a, "A copy: M in the CSR format of the library (from other index types, or from a CSC matrix)");

	m.def("TransposedView", [](nb::handle A)
		{
			if (!HasArrays(A) || Format(A) != "csc")
				throw nb::type_error("fhhos4.TransposedView(): for a CSC matrix (e.g. scipy.sparse.csc_array)");
			return FromAttributes(A, true);
		}, "A"_a, "No copy: the CSR view of A^T from the arrays of a CSC matrix A: A itself if A is symmetric");

	//------------------------------------------------ Solver ------------------------------------------------//

	// The setups and solves release the GIL: they only read the views and the NumPy vectors
	auto noGIL = nb::call_guard<nb::gil_scoped_release>();

	nb::class_<PySolver>(m, "Solver", "fhhos4::Solver: set the parameters, then call a setup function, then Solve()")
		.def(nb::init<>())
		// Discretization
		.def_rw("Dimension", &fhhos4::Solver::Dimension)
		.def_rw("FaceDegree", &fhhos4::Solver::FaceDegree)
		.def_rw("CellDegree", &fhhos4::Solver::CellDegree)
		.def_rw("PolynomialSpace", &fhhos4::Solver::PolynomialSpace)
		.def_rw("SingularMatrix", &fhhos4::Solver::SingularMatrix)
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
		.def_rw("Prolongation", &fhhos4::Solver::Prolongation)
		.def_rw("CoarseningFactor", &fhhos4::Solver::CoarseningFactor)
		.def_rw("CoarseMatrixMaxSize", &fhhos4::Solver::CoarseMatrixMaxSize)
		.def_rw("CoarseSolver", &fhhos4::Solver::CoarseSolver)
		.def_rw("MixedPrecision", &fhhos4::Solver::MixedPrecision)
		.def_rw("BlockCSRStorage", &fhhos4::Solver::BlockCSRStorage)
		// Execution
		.def_rw("Threads", &fhhos4::Solver::Threads)
		.def_rw("Verbosity", &fhhos4::Solver::Verbosity)
#if FHHOS4_MPI
		// The processes of the solver: an mpi4py communicator (MPI.COMM_SELF, MPI.COMM_WORLD...); not set: one process
		.def_prop_rw("Comm", [](const PySolver& s) { return s.CommObject; },
			[](PySolver& s, nb::handle comm) { s.Comm = ToComm(comm); s.CommObject = nb::borrow(comm); })
#endif

		// Setups from the blocks of the system. faceGhosts, cellGhosts: on several processes, the global numbers of the ghost
		// faces and cells (int64 arrays)
		.def("SetupSym", [](PySolver& s, nb::handle A_TT, nb::handle A_TF, nb::handle A_FF, const ConstVector& cellInterpOfOne, const ConstVector& faceInterpOfOne, nb::handle faceGhosts)
			{
				PyCSRView a_TT = ToView(A_TT, "A_TT"), a_TF = ToView(A_TF, "A_TF"), a_FF = ToView(A_FF, "A_FF");
				fhhos4::detail::CSRView tt = a_TT.Native("A_TT"), tf = a_TF.Native("A_TF"), ff = a_FF.Native("A_FF");
				fhhos4::GhostNumbers ghosts = ToGhostNumbers(faceGhosts, "faceGhosts");
				s.Run([&] { s.SetupViews(true, tt, tf, nullptr, ff, cellInterpOfOne, faceInterpOfOne, ghosts); }, {});
			},
			"A_TT"_a, "A_TF"_a, "A_FF"_a, "cellInterpOfOne"_a.noconvert(), "faceInterpOfOne"_a.noconvert(), "faceGhosts"_a = nb::none())
		.def("SetupNonSym", [](PySolver& s, nb::handle A_TT, nb::handle A_TF, nb::handle A_FT, nb::handle A_FF, const ConstVector& cellInterpOfOne, const ConstVector& faceInterpOfOne, nb::handle faceGhosts, nb::handle cellGhosts)
			{
				PyCSRView a_TT = ToView(A_TT, "A_TT"), a_TF = ToView(A_TF, "A_TF"), a_FT = ToView(A_FT, "A_FT"), a_FF = ToView(A_FF, "A_FF");
				fhhos4::detail::CSRView tt = a_TT.Native("A_TT"), tf = a_TF.Native("A_TF"), ft = a_FT.Native("A_FT"), ff = a_FF.Native("A_FF");
				fhhos4::GhostNumbers fGhosts = ToGhostNumbers(faceGhosts, "faceGhosts"), cGhosts = ToGhostNumbers(cellGhosts, "cellGhosts");
				s.Run([&] { s.SetupViews(false, tt, tf, &ft, ff, cellInterpOfOne, faceInterpOfOne, fGhosts, cGhosts); }, { a_FT });
			},
			"A_TT"_a, "A_TF"_a, "A_FT"_a, "A_FF"_a, "cellInterpOfOne"_a.noconvert(), "faceInterpOfOne"_a.noconvert(), "faceGhosts"_a = nb::none(), "cellGhosts"_a = nb::none())
		// x_T, x_F: NumPy vectors of float64, written in place (x_F: the initial guess on input)
		.def("Solve", [](PySolver& s, const ConstVector& b_T, const ConstVector& b_F, Vector x_T, Vector x_F)
			{ return s.Solve(b_T, b_F, x_T, x_F); },
			"b_T"_a.noconvert(), "b_F"_a.noconvert(), "x_T"_a.noconvert(), "x_F"_a.noconvert(), noGIL)

		// The condensed system A x = b (A: CondensedMatrix()), after a setup
		// x: a NumPy vector of float64, written in place (the initial guess on input)
		.def("Solve", [](PySolver& s, const ConstVector& b, Vector x) { return s.Solve(b, x); },
			"b"_a.noconvert(), "x"_a.noconvert(), noGIL)
		// Solves A x = b from x = 0, and returns x
		.def("Solve", [](PySolver& s, const ConstVector& b) { return s.Solve(b); }, "b"_a.noconvert(), noGIL)

		// A copy of the condensed matrix (fhhos4.CSRMatrix)
		.def("CondensedMatrix", [](const PySolver& s) { return fhhos4::detail::ToCSR(s.CondensedView()); })
		.def("IsSetUp", &fhhos4::Solver::IsSetUp)
		.def("NumberOfLevels", &fhhos4::Solver::NumberOfLevels);
}
