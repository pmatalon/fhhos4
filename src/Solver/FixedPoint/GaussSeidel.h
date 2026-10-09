#pragma once
#include "BlockSOR.h"
#include "DiagonalBlockPositions.h"
using namespace std;

class GaussSeidel : public IterativeSolver
{
protected:
	Direction _direction;
	bool _hybrid = false; // see HybridSweep()
	vector<SparseMatrixIndex> _diagPos; // position of the diagonal coefficient of each row in the arrays of the matrix
public:
	GaussSeidel() : GaussSeidel(Direction::Forward) {}

	GaussSeidel(Direction direction, bool hybrid = false)
	{
		this->_direction = direction;
		this->_hybrid = hybrid;
	}

	virtual void Serialize(ostream& os) const override
	{
		if (_hybrid)
			os << "Hybrid ";
		if (_direction == Direction::Forward)
			os << "Gauss-Seidel (forward)";
		else if (_direction == Direction::Backward)
			os << "Gauss-Seidel (backward)";
		else if (_direction == Direction::Symmetric)
			os << "Symmetric Gauss-Seidel";
		else if (_direction == Direction::AlternatingForwardFirst)
			os << "Gauss-Seidel (alternating, forward first)";
		else if (_direction == Direction::AlternatingBackwardFirst)
			os << "Gauss-Seidel (alternating, backward first)";
		else
			os << "Gauss-Seidel (UNDEFINED)";
	}

	//------------------------------------------------//
	// Big assumption: the matrix must be symmetric!  //
	//------------------------------------------------//

	void Setup(const SparseMatrix& A) override
	{
		IterativeSolver::Setup(A);

		// The rows are sorted by column index: the diagonal coefficient separates the strictly lower and upper parts
		DiagonalBlockPositions diagonal;
		diagonal.Setup(A, 1);
		for (BigNumber i = 0; i < A.rows(); ++i)
		{
			if (diagonal.End[i] != diagonal.Begin[i] + 1)
				Utils::FatalError("Gauss-Seidel: no diagonal coefficient in row " + to_string(i) + ".");
		}
		_diagPos = std::move(diagonal.Begin);

		this->SetupComputationalWork = 0;
	}

	// The residual is deduced from the products computed in the sweep (half of the matrix): not in the hybrid sweeps,
	// whose rows don't all read the same values of x
	bool CanOptimizeResidualComputation() override
	{
		return !_hybrid;
	}

private:
	IterationResult ExecuteOneIteration(const Vector& b, Vector& x, bool& xEquals0, bool computeResidual, bool computeAx, const IterationResult& oldResult) override
	{
		IterationResult result(oldResult);
		assert(!_hybrid || (!computeResidual && !computeAx)); // see CanOptimizeResidualComputation()

		const SparseMatrix& A = *this->Matrix;

		if (!computeResidual && !computeAx)
		{
			if (_direction == Direction::Forward)
				ForwardSweep(b, x, xEquals0, result);
			else if (_direction == Direction::Backward)
				BackwardSweep(b, x, xEquals0, result);
			else if (_direction == Direction::Symmetric)
			{
				ForwardSweep(b, x, xEquals0, result);
				BackwardSweep(b, x, xEquals0, result);
			}
			else if (_direction == Direction::AlternatingForwardFirst || _direction == Direction::AlternatingBackwardFirst)
			{
				int modulo = _direction == Direction::AlternatingForwardFirst ? 0 : 1;
				if (this->IterationCount % 2 == modulo)
					ForwardSweep(b, x, xEquals0, result);
				else
					BackwardSweep(b, x, xEquals0, result);
			}
			else
				Utils::FatalError("direction not managed");
		}
		else
		{
			if (_direction == Direction::Forward)
				ForwardSweepAndComputeResidualOrAx(b, x, xEquals0, computeAx, result);
			else if (_direction == Direction::Backward)
				BackwardSweepAndComputeResidualOrAx(b, x, xEquals0, computeResidual, computeAx, result);
			else if (_direction == Direction::Symmetric)
			{
				ForwardSweep(b, x, xEquals0, result);
				BackwardSweepAndComputeResidualOrAx(b, x, xEquals0, computeResidual, computeAx, result);
			}
			else if (_direction == Direction::AlternatingForwardFirst || _direction == Direction::AlternatingBackwardFirst)
			{
				int modulo = _direction == Direction::AlternatingForwardFirst ? 0 : 1;
				if (this->IterationCount % 2 == modulo)
					ForwardSweepAndComputeResidualOrAx(b, x, xEquals0, computeAx, result);
				else
					BackwardSweepAndComputeResidualOrAx(b, x, xEquals0, computeResidual, computeAx, result);
			}
			else
				Utils::FatalError("direction not managed");
		}

		result.SetX(x);
		return result;
	}

	void ForwardSweep(const Vector& b, Vector& x, bool& xEquals0, IterationResult& result)
	{
		const SparseMatrix& A = *this->Matrix;
		// x(new) = (L+D)^{-1} * (b-Ux)
		if (!xEquals0)
			result.AddWorkInFlops(Cost::DAXPY_StrictTri(A));
		if (_hybrid)
			HybridSweep(b, x, xEquals0, false);
		else
			ForwardSweepKernel(b, x, xEquals0, nullptr);
		                                                                          result.AddWorkInFlops(Cost::SpFWElimination(A));
		xEquals0 = false;
	}

	void BackwardSweep(const Vector& b, Vector& x, bool& xEquals0, IterationResult& result)
	{
		const SparseMatrix& A = *this->Matrix;
		// x(new) = (D+U)^{-1} * (b-Lx)
		if (!xEquals0)
			result.AddWorkInFlops(Cost::DAXPY_StrictTri(A));
		if (_hybrid)
			HybridSweep(b, x, xEquals0, true);
		else
			BackwardSweepKernel(b, x, xEquals0, nullptr);
		                                                                          result.AddWorkInFlops(Cost::SpBWSubstitution(A));
		xEquals0 = false;
	}

	// The residual is always computed (it is cheaper than Ax), and Ax is deduced from it if requested.
	void ForwardSweepAndComputeResidualOrAx(const Vector& b, Vector& x, bool& xEquals0, bool computeAx, IterationResult& result)
	{
		const SparseMatrix& A = *this->Matrix;

		if (!xEquals0)
		{
			// Sweep: x(new) = (L+D)^{-1} * (b-Ux)
			Vector Ux(A.rows());
			ForwardSweepKernel(b, x, false, &Ux);                 result.AddWorkInFlops(Cost::SpMatVec(NNZ::StrictTriPart(A)));
			                                                      result.AddWorkInFlops(Cost::AddVec(b) + Cost::SpFWElimination(A));
			xEquals0 = false;

			// Residual: Ux(old) - Ux(new)
			result.Residual = Vector(A.rows());
			StrictUpperProduct(x, &Ux, result.Residual);         result.AddWorkInFlops(Cost::DAXPY_StrictTri(A));
		}
		else
		{
			// Sweep
			ForwardSweepKernel(b, x, true, nullptr);              result.AddWorkInFlops(Cost::SpFWElimination(A));
			xEquals0 = false;

			// Residual: r = -Ux
			result.Residual = Vector(A.rows());
			StrictUpperProduct(x, nullptr, result.Residual);      result.AddWorkInFlops(Cost::DAXPY_StrictTri(A));
		}

		if (computeAx)
		{
			result.Ax = b - result.Residual;                      result.AddWorkInFlops(Cost::AddVec(b));
		}
	}

	void BackwardSweepAndComputeResidualOrAx(const Vector& b, Vector& x, bool& xEquals0, bool computeResidual, bool computeAx, IterationResult& result)
	{
		const SparseMatrix& A = *this->Matrix;

		if (!xEquals0)
		{
			// Sweep: x(new) = (D+U)^{-1} * (b-Lx)
			Vector b_Lx(A.rows());
			BackwardSweepKernel(b, x, false, &b_Lx);              result.AddWorkInFlops(Cost::DAXPY_StrictTri(A));
			                                                      result.AddWorkInFlops(Cost::SpBWSubstitution(A));
			xEquals0 = false;

			if (computeAx || computeResidual)
			{
				result.Ax = Vector(A.rows());
				StrictLowerProduct(x, &b_Lx, result.Ax);          result.AddWorkInFlops(Cost::DAXPY_StrictTri(A));
				if (computeResidual)
				{
					result.Residual = b - result.Ax;              result.AddWorkInFlops(Cost::AddVec(b));
				}
			}
		}
		else
		{
			// Sweep: x(new) = (D+U)^{-1} * b
			BackwardSweepKernel(b, x, true, nullptr);             result.AddWorkInFlops(Cost::SpBWSubstitution(A));
			xEquals0 = false;

			if (computeAx)
			{
				// Ax = b + Lx
				result.Ax = Vector(A.rows());
				StrictLowerProduct(x, &b, result.Ax);             result.AddWorkInFlops(Cost::DAXPY_StrictTri(A));
				if (computeResidual)
				{
					result.Residual = b - result.Ax;              result.AddWorkInFlops(Cost::AddVec(b));
				}
			}
			else if (computeResidual)
			{
				// Residual: Lx(old) - Lx(new) = -Lx
				result.Residual = Vector(A.rows());
				StrictLowerProduct(x, nullptr, result.Residual);  result.AddWorkInFlops(Cost::DAXPY_StrictTri(A));
			}
		}
	}

	// The kernels below compute the same operations, in the same order, as the Eigen expressions they replace
	// (A.triangularView<...>().solve(b), A.triangularView<...>() * x): same results (up to the sign of zeros). But a
	// sweep reads the rows of the matrix once (Eigen's b - U*x followed by the solve reads them twice), and finds their
	// diagonal directly.
	// - triangular solve: x_i = (b_i - sum_j a_ij x_j) / a_ii, the products subtracted one by one, in the order of the row;
	// - sparse matrix-vector product: for each row, two partial sums (even and odd positions) added at the end.

	// Sum of the a_ij*x_j over the positions [begin, end) of the arrays of the matrix, as computed by Eigen's sparse
	// matrix-vector product
	static double RowProduct(const double* val, const SparseMatrixIndex* col, const double* x, SparseMatrixIndex begin, SparseMatrixIndex end)
	{
		double sumEven = 0;
		double sumOdd = 0;
		SparseMatrixIndex p = begin;
		for (; p + 1 < end; p += 2)
		{
			sumEven += val[p] * x[col[p]];
			sumOdd += val[p + 1] * x[col[p + 1]];
		}
		if (p < end)
			sumEven += val[p] * x[col[p]];
		return sumEven + sumOdd;
	}

	// x = (L+D)^{-1} * (b - Ux): Eigen's L_plus_D.solve(b - U*x), or L_plus_D.solve(b) if xEquals0.
	// If Ux is given (and !xEquals0), Ux = U*x (old x).
	void ForwardSweepKernel(const Vector& b, Vector& x, bool xEquals0, Vector* Ux)
	{
		BigNumber n = this->Matrix->rows();
		if (xEquals0)
			x.resize(n);
		ForwardSweepKernel(b, x, xEquals0, Ux, 0, n);
	}

	// The same on the rows [begin, end) only
	void ForwardSweepKernel(const Vector& b, Vector& x, bool xEquals0, Vector* Ux, BigNumber begin, BigNumber end)
	{
		const SparseMatrix& A = *this->Matrix;
		const SparseMatrixIndex* outer = A.outerIndexPtr();
		const SparseMatrixIndex* col = A.innerIndexPtr();
		const double* val = A.valuePtr();
		double* xp = x.data();

		for (BigNumber i = begin; i < end; ++i)
		{
			SparseMatrixIndex d = _diagPos[i];
			double tmp = b[i];
			if (!xEquals0)
			{
				double Ux_i = RowProduct(val, col, xp, d + 1, outer[i + 1]); // x_j, j > i: old values
				if (Ux)
					(*Ux)[i] = Ux_i;
				tmp -= Ux_i;
			}
			for (SparseMatrixIndex p = outer[i]; p < d; ++p) // x_j, j < i: new values
				tmp -= val[p] * xp[col[p]];
			xp[i] = tmp / val[d];
		}
	}

	// x = (D+U)^{-1} * (b - Lx): Eigen's D_plus_U.solve(b - L*x), or D_plus_U.solve(b) if xEquals0.
	// If b_Lx is given (and !xEquals0), b_Lx = b - L*x (old x).
	void BackwardSweepKernel(const Vector& b, Vector& x, bool xEquals0, Vector* b_Lx)
	{
		BigNumber n = this->Matrix->rows();
		if (xEquals0)
			x.resize(n);
		BackwardSweepKernel(b, x, xEquals0, b_Lx, 0, n);
	}

	// The same on the rows [begin, end) only
	void BackwardSweepKernel(const Vector& b, Vector& x, bool xEquals0, Vector* b_Lx, BigNumber begin, BigNumber end)
	{
		const SparseMatrix& A = *this->Matrix;
		const SparseMatrixIndex* outer = A.outerIndexPtr();
		const SparseMatrixIndex* col = A.innerIndexPtr();
		const double* val = A.valuePtr();
		double* xp = x.data();

		for (BigNumber i = end; i-- > begin; )
		{
			SparseMatrixIndex d = _diagPos[i];
			double tmp = b[i];
			if (!xEquals0)
			{
				tmp -= RowProduct(val, col, xp, outer[i], d); // x_j, j < i: old values
				if (b_Lx)
					(*b_Lx)[i] = tmp;
			}
			for (SparseMatrixIndex p = d + 1; p < outer[i + 1]; ++p) // x_j, j > i: new values
				tmp -= val[p] * xp[col[p]];
			xp[i] = tmp / val[d];
		}
	}

	// Hybrid sweep: each thread sweeps its own contiguous chunk of rows (Gauss-Seidel), and reads the values of x of the
	// other chunks as they are when it reads them (Jacobi between the chunks, up to the timing of the threads: x is read
	// by a thread while another one writes it). The result depends on the number of threads and varies from one run to
	// the next, as with BlockGaussSeidel's hybrid sweeps. From x = 0, the products by the values not computed yet in the
	// sweep order are skipped: in the other chunks, those are then the old values (0), whatever the timing.
	void HybridSweep(const Vector& b, Vector& x, bool xEquals0, bool backward)
	{
		const SparseMatrix& A = *this->Matrix;
		BigNumber n = A.rows();
		if (xEquals0 && x.size() != n)
			x = Vector::Zero(n); // read by the other chunks
		#pragma omp parallel if (A.nonZeros() > 20000)
		{
			auto [begin, end] = Parallelism::ThreadChunk(n);
			if (backward)
				BackwardSweepKernel(b, x, xEquals0, nullptr, begin, end);
			else
				ForwardSweepKernel(b, x, xEquals0, nullptr, begin, end);
		}
	}

	// result = v - U*x, or result = -U*x if v is null.
	// The rows are independent: parallel like Eigen's product, with the same result.
	void StrictUpperProduct(const Vector& x, const Vector* v, Vector& result) const
	{
		const SparseMatrix& A = *this->Matrix;
		const SparseMatrixIndex* outer = A.outerIndexPtr();
		const SparseMatrixIndex* col = A.innerIndexPtr();
		const double* val = A.valuePtr();
		const double* xp = x.data();
		BigNumber n = A.rows();
		if (v)
		{
			#pragma omp parallel for if (A.nonZeros() > 20000)
			for (BigNumber i = 0; i < n; ++i)
				result[i] = (*v)[i] - RowProduct(val, col, xp, _diagPos[i] + 1, outer[i + 1]);
		}
		else
		{
			#pragma omp parallel for if (A.nonZeros() > 20000)
			for (BigNumber i = 0; i < n; ++i)
				result[i] = -RowProduct(val, col, xp, _diagPos[i] + 1, outer[i + 1]);
		}
	}

	// result = v + L*x, or result = -L*x if v is null.
	void StrictLowerProduct(const Vector& x, const Vector* v, Vector& result) const
	{
		const SparseMatrix& A = *this->Matrix;
		const SparseMatrixIndex* outer = A.outerIndexPtr();
		const SparseMatrixIndex* col = A.innerIndexPtr();
		const double* val = A.valuePtr();
		const double* xp = x.data();
		BigNumber n = A.rows();
		if (v)
		{
			#pragma omp parallel for if (A.nonZeros() > 20000)
			for (BigNumber i = 0; i < n; ++i)
				result[i] = (*v)[i] + RowProduct(val, col, xp, outer[i], _diagPos[i]);
		}
		else
		{
			#pragma omp parallel for if (A.nonZeros() > 20000)
			for (BigNumber i = 0; i < n; ++i)
				result[i] = -RowProduct(val, col, xp, outer[i], _diagPos[i]);
		}
	}
};
