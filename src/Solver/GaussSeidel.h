#pragma once
#include "BlockSOR.h"
using namespace std;

class GaussSeidel : public IterativeSolver
{
protected:
	Direction _direction;
	RowMajorSparseMatrix _rowMajorA;
	bool _useEigen = true;
public:
	GaussSeidel() : GaussSeidel(Direction::Forward) {}

	GaussSeidel(Direction direction)
	{
		this->_direction = direction;
	}

	virtual void Serialize(ostream& os) const override
	{
		if (_direction == Direction::Symmetric)
			os << "Symmetric ";
		os << "Gauss-Seidel";

		if (_direction != Direction::Symmetric)
			os << " (direction=" << (_direction == Direction::Forward ? "forward" : "backward") << ")";
	}

	//------------------------------------------------//
	// Big assumption: the matrix must be symmetric!  //
	//------------------------------------------------//

	void Setup(const SparseMatrix& A) override
	{
		IterativeSolver::Setup(A);
		if (!A.IsRowMajor)
			this->_rowMajorA = A;

		this->SetupComputationalWork = 0;
	}

	bool CanOptimizeResidualComputation() override
	{
		return true;
	}

private:
	IterationResult ExecuteOneIteration(const Vector& b, Vector& x, bool& xEquals0, bool computeResidual, bool computeAx, const IterationResult& oldResult) override
	{
		IterationResult result(oldResult);

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
			else
				Utils::FatalError("direction not managed");
		}
		else
		{
			if (_direction == Direction::Forward)
				ForwardSweepAndComputeResidual(b, x, xEquals0, result);
			else if (_direction == Direction::Backward)
				BackwardSweepAndComputeResidualOrAx(b, x, xEquals0, computeResidual, computeAx, result);
			else if (_direction == Direction::Symmetric)
			{
				ForwardSweep(b, x, xEquals0, result);
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
		if (_useEigen)
		{
			if (!xEquals0)
			{
				Vector v = b - A.triangularView<Eigen::StrictlyUpper>() * x;      result.AddCost(Cost::DAXPY_StrictTri(A));
				x = A.triangularView<Eigen::Lower>().solve(v);                    result.AddCost(Cost::SpFWElimination(A));
			}
			else
			{
				x = A.triangularView<Eigen::Lower>().solve(b);                    result.AddCost(Cost::SpFWElimination(A));
			}
		}
		else
		{
			for (BigNumber i = 0; i < A.rows(); ++i)
				ProcessRow(i, b, x);
			                                                                      result.AddCost(2 * A.nonZeros() + A.rows());
		}
		xEquals0 = false;
	}

	void BackwardSweep(const Vector& b, Vector& x, bool& xEquals0, IterationResult& result)
	{
		const SparseMatrix& A = *this->Matrix;
		if (_useEigen)
		{
			if (!xEquals0)
			{
				Vector v = b - A.triangularView<Eigen::StrictlyLower>() * x;      result.AddCost(Cost::DAXPY_StrictTri(A));
				x = A.triangularView<Eigen::Upper>().solve(v);                    result.AddCost(Cost::SpBWSubstitution(A));
			}
			else
			{
				x = A.triangularView<Eigen::Upper>().solve(b);                    result.AddCost(Cost::SpBWSubstitution(A));
			}
		}
		else
		{
			for (BigNumber i = 0; i < A.rows(); ++i)
				ProcessRow(A.rows() - i - 1, b, x);
			                                                                      result.AddCost(2 * A.nonZeros() + A.rows());
		}
		xEquals0 = false;
	}

	void ForwardSweepAndComputeResidual(const Vector& b, Vector& x, bool& xEquals0, IterationResult& result)
	{
		const SparseMatrix& A = *this->Matrix;
		auto U = A.triangularView<Eigen::StrictlyUpper>();
		auto L_plus_D = A.triangularView<Eigen::Lower>();

		if (!xEquals0)
		{
			// Sweep: x(new) = (L+D)^{-1} * (b-Ux)
			Vector Ux = U * x;                                    result.AddCost(Cost::SpMatVec(NNZ::StrictTriPart(A)));
			x = L_plus_D.solve(b - Ux);                           result.AddCost(Cost::AddVec(b) + Cost::SpFWElimination(A));
			xEquals0 = false;

			// Residual: Ux(old) - Ux(new)
			result.Residual = Ux - U * x;                         result.AddCost(Cost::DAXPY_StrictTri(A));
		}
		else
		{
			// Sweep
			x = L_plus_D.solve(b);                                result.AddCost(Cost::SpFWElimination(A));
			xEquals0 = false;

			// Residual: r = -Ux
			result.Residual = -U * x;                                          result.AddCost(Cost::DAXPY_StrictTri(A));
		}
	}

	void BackwardSweepAndComputeResidualOrAx(const Vector& b, Vector& x, bool& xEquals0, bool computeResidual, bool computeAx, IterationResult& result)
	{
		const SparseMatrix& A = *this->Matrix;
		auto L = A.triangularView<Eigen::StrictlyLower>();
		auto D_plus_U = A.triangularView<Eigen::Upper>();

		if (!xEquals0)
		{
			// Sweep: x(new) = (D+U)^{-1} * (b-Lx)
			Vector b_Lx = b - L*x;                                result.AddCost(Cost::DAXPY_StrictTri(A));
			x = D_plus_U.solve(b_Lx);                             result.AddCost(Cost::SpBWSubstitution(A));
			xEquals0 = false;

			if (computeAx || computeResidual)
			{
				result.Ax = b_Lx + L * x;                         result.AddCost(Cost::DAXPY_StrictTri(A));
				if (computeResidual)
				{
					result.Residual = b - result.Ax;              result.AddCost(Cost::AddVec(b));
				}
			}
		}
		else
		{
			// Sweep: x(new) = (D+U)^{-1} * b
			x = D_plus_U.solve(b);                                result.AddCost(Cost::SpBWSubstitution(A));
			xEquals0 = false;

			if (computeAx)
			{
				result.Ax = b + L * x;                            result.AddCost(Cost::DAXPY_StrictTri(A));
				if (computeResidual)
				{
					result.Residual = b - result.Ax;              result.AddCost(Cost::AddVec(b));
				}
			}
			else if (computeResidual)
			{
				// Residual: Lx(old) - Lx(new)
				result.Residual = -L * x;                         result.AddCost(Cost::DAXPY_StrictTri(A));
			}
		}
	}

	// Same as ProcessBlockRow, but optimized for blockSize = 1
	void ProcessRow(BigNumber currentRow, const Vector& b, Vector& x)
	{
		BigNumber i = currentRow;
		const SparseMatrix& A = *this->Matrix;
		double tmp_x = b(i);
		double a_ii = 0;

		// RowMajor --> the following line iterates over the non-zeros of the i-th row.
		for (RowMajorSparseMatrix::InnerIterator it(A, i); it; ++it)
		{
			BigNumber j = it.col();
			double a_ij = it.value();
			if (i == j) // Di
				a_ii = a_ij;
			else // Li and Ui
				tmp_x -= a_ij * x(j);
		}

		x(i) = tmp_x / a_ii;
	}
};