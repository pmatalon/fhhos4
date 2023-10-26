#pragma once
#include "BlockSOR.h"
using namespace std;

class GaussSeidel : public IterativeSolver
{
protected:
	Direction _direction;
	RowMajorSparseMatrix _rowMajorA;
	bool _useEigen = true;
	Vector invD; // inverses of the diagonal elements
public:
	GaussSeidel() : GaussSeidel(Direction::Forward) {}

	GaussSeidel(Direction direction)
	{
		this->_direction = direction;
	}

	virtual void Serialize(ostream& os) const override
	{
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
		if (!A.IsRowMajor)
			this->_rowMajorA = A;

		/*this->invD = Vector(A.rows());
		for (BigNumber i = 0; i < A.rows(); ++i)
			this->invD[i] = 1.0 / A.coeff(i, i);*/

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
				ForwardSweepAndComputeResidual(b, x, xEquals0, result);
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
					ForwardSweepAndComputeResidual(b, x, xEquals0, result);
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
		if (_useEigen)
		{
			if (!xEquals0)
			{
				Vector v = b - A.triangularView<Eigen::StrictlyUpper>() * x;      result.AddWorkInFlops(Cost::DAXPY_StrictTri(A));
				x = A.triangularView<Eigen::Lower>().solve(v);                    result.AddWorkInFlops(Cost::SpFWElimination(A));
			}
			else
			{
				x = A.triangularView<Eigen::Lower>().solve(b);                    result.AddWorkInFlops(Cost::SpFWElimination(A));
			}
		}
		else
		{
			for (BigNumber i = 0; i < A.rows(); ++i)
				ProcessRow(i, b, x);
			                                                                      result.AddWorkInFlops(2 * A.nonZeros() + A.rows());
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
				Vector v = b - A.triangularView<Eigen::StrictlyLower>() * x;      result.AddWorkInFlops(Cost::DAXPY_StrictTri(A));
				x = A.triangularView<Eigen::Upper>().solve(v);                    result.AddWorkInFlops(Cost::SpBWSubstitution(A));
			}
			else
			{
				x = A.triangularView<Eigen::Upper>().solve(b);                    result.AddWorkInFlops(Cost::SpBWSubstitution(A));
			}
		}
		else
		{
			for (BigNumber i = 0; i < A.rows(); ++i)
				ProcessRow(A.rows() - i - 1, b, x);
			                                                                      result.AddWorkInFlops(2 * A.nonZeros() + A.rows());
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
			Vector Ux = U * x;                                    result.AddWorkInFlops(Cost::SpMatVec(NNZ::StrictTriPart(A)));
			x = L_plus_D.solve(b - Ux);                           result.AddWorkInFlops(Cost::AddVec(b) + Cost::SpFWElimination(A));
			xEquals0 = false;

			// Residual: Ux(old) - Ux(new)
			result.Residual = Ux - U * x;                         result.AddWorkInFlops(Cost::DAXPY_StrictTri(A));
		}
		else
		{
			// Sweep
			x = L_plus_D.solve(b);                                result.AddWorkInFlops(Cost::SpFWElimination(A));
			//ForwardElimination(x, b);                            result.AddWorkInFlops(Cost::SpFWElimination(A));
			xEquals0 = false;

			// Residual: r = -Ux
			result.Residual = -U * x;                                          result.AddWorkInFlops(Cost::DAXPY_StrictTri(A));
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
			Vector b_Lx = b - L*x;                                result.AddWorkInFlops(Cost::DAXPY_StrictTri(A));
			x = D_plus_U.solve(b_Lx);                             result.AddWorkInFlops(Cost::SpBWSubstitution(A));
			//BackwardSubstitution(x, b_Lx);                        result.AddWorkInFlops(Cost::SpBWSubstitution(A));
			xEquals0 = false;

			if (computeAx || computeResidual)
			{
				result.Ax = b_Lx + L * x;                         result.AddWorkInFlops(Cost::DAXPY_StrictTri(A));
				if (computeResidual)
				{
					result.Residual = b - result.Ax;              result.AddWorkInFlops(Cost::AddVec(b));
				}
			}
		}
		else
		{
			// Sweep: x(new) = (D+U)^{-1} * b
			x = D_plus_U.solve(b);                                result.AddWorkInFlops(Cost::SpBWSubstitution(A));
			xEquals0 = false;

			if (computeAx)
			{
				result.Ax = b + L * x;                            result.AddWorkInFlops(Cost::DAXPY_StrictTri(A));
				if (computeResidual)
				{
					result.Residual = b - result.Ax;              result.AddWorkInFlops(Cost::AddVec(b));
				}
			}
			else if (computeResidual)
			{
				// Residual: Lx(old) - Lx(new)
				result.Residual = -L * x;                         result.AddWorkInFlops(Cost::DAXPY_StrictTri(A));
			}
		}
	}

	// Solves (L+D)x=b, where L is the strictly lower triangular part of A and D the diagonal
	void ForwardElimination(Vector& x, const Vector& b)
	{
		const SparseMatrix& A = *this->Matrix;
		//Vector x(A.rows());

		for (BigNumber i = 0; i < A.rows(); ++i)
		{
			double tmp_x = b(i);

			// RowMajor --> the following line iterates over the non-zeros of the i-th row.
			for (RowMajorSparseMatrix::InnerIterator it(A, i); it; ++it)
			{
				BigNumber j = it.col(); 
				if (i == j) // Ui or Di
					break;
				assert(i > j);
				// Li
				double a_ij = it.value();
				tmp_x -= a_ij * x(j);
			}

			x(i) = tmp_x * this->invD[i];
		}
		//return x;
	}

	// Solves (D+U)x=b, where U is the strictly upper triangular part of A and D the diagonal
	void BackwardSubstitution(Vector& x, const Vector& b)
	{
		const SparseMatrix& A = *this->Matrix;
		//Vector x(A.rows());

		for (BigNumber i = A.rows()-1; i >= 0; --i)
		{
			double tmp_x = b(i);

			// RowMajor --> the following line iterates over the non-zeros of the i-th row.
			for (RowMajorSparseMatrix::InnerIterator it(A, i); it; ++it)
			{
				BigNumber j = it.col();
				if (i < j) // U
				{
					double a_ij = it.value();
					tmp_x -= a_ij * x(j);
				}
			}

			x(i) = tmp_x * this->invD[i];
		}
		//return x;
	}

	
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