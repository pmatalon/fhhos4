#pragma once
#include "BlockSOR.h"
using namespace std;

class GaussSeidel : public IterativeSolver
{
protected:
	Direction _direction;

	RowMajorSparseMatrix _rowMajorA;
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

		auto nb = A.rows();
		if (_direction == Direction::Forward)
		{
			for (BigNumber i = 0; i < nb; ++i)
				ProcessRow(i, b, x);
		}
		else if (_direction == Direction::Backward)
		{
			for (BigNumber i = 0; i < nb; ++i)
				ProcessRow(nb - i - 1, b, x);
		}
		else // Symmetric
		{
			// Forward
			for (BigNumber i = 0; i < nb; ++i)
				ProcessRow(i, b, x);
			// Backward
			for (BigNumber i = 0; i < nb; ++i)
				ProcessRow(nb - i - 1, b, x);
		}

		result.SetX(x);
		double sweepCost = 2 * A.nonZeros() + nb;
		result.AddCost(_direction == Direction::Symmetric ? 2 * sweepCost : sweepCost);
		return result;
	}

	// Same as ProcessBlockRow, but optimized for blockSize = 1
	void ProcessRow(BigNumber currentRow, const Vector& b, Vector& x)
	{
		BigNumber i = currentRow;
		const SparseMatrix& A = *this->Matrix;
		double tmp_x = b(i);
		double a_ii = 0;

		// RowMajor --> the following line iterates over the non-zeros of the i-th row.
		for (RowMajorSparseMatrix::InnerIterator it(A.IsRowMajor ? A : _rowMajorA, i); it; ++it)
		{
			auto j = it.col();
			auto a_ij = it.value();
			if (i == j) // Di
				a_ii = a_ij;
			else // Li and Ui
				tmp_x -= a_ij * x(j);
		}

		x(i) = tmp_x / a_ii;
	}
};