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
	IterationResult ExecuteOneIteration(const Vector& b, Vector& x, const IterationResult& oldResult) override
	{
		IterationResult result(oldResult);

		const SparseMatrix& A = *this->Matrix;

		if (_direction == Direction::Forward)
			ForwardSweep(b, x);
		else if (_direction == Direction::Backward)
			BackwardSweep(b, x);
		else if (_direction == Direction::Symmetric)
		{
			ForwardSweep(b, x);
			BackwardSweep(b, x);
		}
		else
			Utils::FatalError("direction not managed");

		result.SetX(x);
		double sweepCost = 2 * A.nonZeros() + A.rows(); // 1 DAXPY
		if (_direction == Direction::Symmetric)
			sweepCost *= 2;
		result.AddCost(sweepCost);
		return result;
	}

	pair<IterationResult, Vector> ExecuteOneIterationAndComputeResidual(const Vector& b, Vector& x, const IterationResult& oldResult) override
	{
		pair<IterationResult, Vector> p;
		auto&[result, r] = p;

		result = IterationResult(oldResult);

		const SparseMatrix& A = *this->Matrix;

		if (_direction == Direction::Forward)
			r = ForwardSweepAndComputeResidual(b, x);
		else if (_direction == Direction::Backward)
			r = BackwardSweepAndComputeResidual(b, x);
		else if (_direction == Direction::Symmetric)
		{
			ForwardSweep(b, x);
			r = BackwardSweepAndComputeResidual(b, x);
		}
		else
			Utils::FatalError("direction not managed");

		result.SetX(x);

		double sweepCost = 2 * A.nonZeros() + A.rows();
		if (_direction == Direction::Symmetric)
			sweepCost *= 2;
		double residualCost = Cost::DAXPY(NNZ::StrictTriPart(A), A.rows());
		result.AddCost(sweepCost + residualCost);

		return p;
	}

	void ForwardSweep(const Vector& b, Vector& x)
	{
		const SparseMatrix& A = *this->Matrix;
		if (_useEigen)
		{
			Vector v = b - A.triangularView<Eigen::StrictlyUpper>() * x;
			x = A.triangularView<Eigen::Lower>().solve(v);
		}
		else
		{
			for (BigNumber i = 0; i < A.rows(); ++i)
				ProcessRow(i, b, x);
		}
	}

	Vector ForwardSweepAndComputeResidual(const Vector& b, Vector& x)
	{
		const SparseMatrix& A = *this->Matrix;
		auto U = A.triangularView<Eigen::StrictlyUpper>();
		auto L_plus_D = A.triangularView<Eigen::Lower>();

		// Sweep: x(new) = (L+D)^{-1} * (b-Ux)
		Vector Ux = U * x;                        // Cost::MatVec(NNZ::StrictTriPart(A), A.rows());
		x = L_plus_D.solve(b - Ux);               // Cost::AddVec(b) + Cost::SpFWElimination(NNZ::TriPart(A))

		// Residual: Ux(old) - Ux(new)
		Vector r = Ux - U * x;                    // Cost::DAXPY(NNZ::StrictTriPart(A), A.rows())
		return r;
	}

	void BackwardSweep(const Vector& b, Vector& x)
	{
		const SparseMatrix& A = *this->Matrix;
		if (_useEigen)
		{
			Vector v = b - A.triangularView<Eigen::StrictlyLower>() * x;
			x = A.triangularView<Eigen::Upper>().solve(v);
		}
		else
		{
			for (BigNumber i = 0; i < A.rows(); ++i)
				ProcessRow(A.rows() - i - 1, b, x);
		}
	}

	Vector BackwardSweepAndComputeResidual(const Vector& b, Vector& x)
	{
		const SparseMatrix& A = *this->Matrix;
		auto L = A.triangularView<Eigen::StrictlyLower>();
		auto D_plus_U = A.triangularView<Eigen::Upper>();

		// Sweep: x(new) = (D+U)^{-1} * (b-Lx)
		Vector Lx = L * x;
		x = D_plus_U.solve(b - Lx);

		// Residual: Lx(old) - Lx(new)
		Vector r = Lx - L * x;
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