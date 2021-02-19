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

private:
	IterationResult ExecuteOneIteration(const Vector& b, Vector& x, bool& xEquals0, const IterationResult& oldResult) override
	{
		IterationResult result(oldResult);

		const SparseMatrix& A = *this->Matrix;

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

		return result;
	}

	pair<IterationResult, Vector> ExecuteOneIterationAndComputeResidual(const Vector& b, Vector& x, bool& xEquals0, const IterationResult& oldResult) override
	{
		pair<IterationResult, Vector> p;
		auto&[result, r] = p;

		result = IterationResult(oldResult);

		const SparseMatrix& A = *this->Matrix;

		if (_direction == Direction::Forward)
			r = ForwardSweepAndComputeResidual(b, x, xEquals0, result);
		else if (_direction == Direction::Backward)
			r = BackwardSweepAndComputeResidual(b, x, xEquals0, result);
		else if (_direction == Direction::Symmetric)
		{
			ForwardSweep(b, x, xEquals0, result);
			r = BackwardSweepAndComputeResidual(b, x, xEquals0, result);
		}
		else
			Utils::FatalError("direction not managed");

		result.SetX(x);
		return p;
	}

	void ForwardSweep(const Vector& b, Vector& x, bool& xEquals0, IterationResult& result)
	{
		const SparseMatrix& A = *this->Matrix;
		if (_useEigen)
		{
			if (!xEquals0)
			{
				Vector v = b - A.triangularView<Eigen::StrictlyUpper>() * x;      result.AddCost(Cost::DAXPY(NNZ::StrictTriPart(A), A.rows()));
				x = A.triangularView<Eigen::Lower>().solve(v);                    result.AddCost(Cost::SpFWElimination(NNZ::TriPart(A)));
			}
			else
			{
				x = A.triangularView<Eigen::Lower>().solve(b);                    result.AddCost(Cost::SpFWElimination(NNZ::TriPart(A)));
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
				Vector v = b - A.triangularView<Eigen::StrictlyLower>() * x;      result.AddCost(Cost::DAXPY(NNZ::StrictTriPart(A), A.rows()));
				x = A.triangularView<Eigen::Upper>().solve(v);                    result.AddCost(Cost::SpBWSubstitution(NNZ::TriPart(A)));
			}
			else
			{
				x = A.triangularView<Eigen::Upper>().solve(b);                    result.AddCost(Cost::SpBWSubstitution(NNZ::TriPart(A)));
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

	Vector ForwardSweepAndComputeResidual(const Vector& b, Vector& x, bool& xEquals0, IterationResult& result)
	{
		const SparseMatrix& A = *this->Matrix;
		auto U = A.triangularView<Eigen::StrictlyUpper>();
		auto L_plus_D = A.triangularView<Eigen::Lower>();

		Vector r;
		if (!xEquals0)
		{
			// Sweep: x(new) = (L+D)^{-1} * (b-Ux)
			Vector Ux = U * x;                                    result.AddCost(Cost::SpMatVec(NNZ::StrictTriPart(A)));
			x = L_plus_D.solve(b - Ux);                           result.AddCost(Cost::AddVec(b) + Cost::SpFWElimination(NNZ::TriPart(A)));
			xEquals0 = false;

			// Residual: Ux(old) - Ux(new)
			r = Ux - U * x;                                       result.AddCost(Cost::DAXPY(NNZ::StrictTriPart(A), A.rows()));
		}
		else
		{
			// Sweep
			x = L_plus_D.solve(b);                                result.AddCost(Cost::SpFWElimination(NNZ::TriPart(A)));
			xEquals0 = false;

			// Residual
			r = - U * x;                                          result.AddCost(Cost::DAXPY(NNZ::StrictTriPart(A), A.rows()));
		}
		return r;
	}

	Vector BackwardSweepAndComputeResidual(const Vector& b, Vector& x, bool& xEquals0, IterationResult& result)
	{
		const SparseMatrix& A = *this->Matrix;
		auto L = A.triangularView<Eigen::StrictlyLower>();
		auto D_plus_U = A.triangularView<Eigen::Upper>();

		Vector r;
		if (!xEquals0)
		{
			// Sweep: x(new) = (D+U)^{-1} * (b-Lx)
			Vector Lx = L * x;                                    result.AddCost(Cost::SpMatVec(NNZ::StrictTriPart(A)));
			x = D_plus_U.solve(b - Lx);                           result.AddCost(Cost::AddVec(b) + Cost::SpBWSubstitution(NNZ::TriPart(A)));
			xEquals0 = false;

			// Residual: Lx(old) - Lx(new)
			r = Lx - L * x;                                       result.AddCost(Cost::DAXPY(NNZ::StrictTriPart(A), A.rows()));
		}
		else
		{
			// Sweep
			x = D_plus_U.solve(b);                                result.AddCost(Cost::SpBWSubstitution(NNZ::TriPart(A)));
			xEquals0 = false;

			// Residual: Lx(old) - Lx(new)
			r = - L * x;                                          result.AddCost(Cost::DAXPY(NNZ::StrictTriPart(A), A.rows()));
		}
		return r;
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