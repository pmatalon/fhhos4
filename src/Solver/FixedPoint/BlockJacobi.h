#pragma once
#include "../IterativeSolver.h"
#include "BlockDiagonalSolver.h"
using namespace std;

class BlockJacobi : public IterativeSolver
{
protected:
	int _blockSize;
	double _omega;

	BlockDiagonalSolver diagBlockSolver; // solves the systems with the diagonal blocks
	RowMajorSparseMatrix _rowMajorA;
public:

	BlockJacobi(int blockSize, double omega = 1)
	{
		assert(omega > 0 && omega < 2);
		this->_blockSize = blockSize;
		this->_omega = omega;
	}

	virtual void Serialize(ostream& os) const override
	{
		if (_blockSize == 1)
		{
			os << "Jacobi";
			if (_omega != 1)
			{
				os << " (omega=";
				if (_omega == 2.0/3.0)
					os << "2/3";
				else
					os << _omega;
				os << ")";
			}
		}
		else
		{
			os << "Block Jacobi (blockSize=" << _blockSize;
			if (_omega != 1)
			{
				os << ", omega=";
				if (_omega == 2.0 / 3.0)
					os << "2/3";
				else
				os <<  _omega;
			}
			os << ")";
		}
	}

	//------------------------------------------------//
	// Big assumption: the matrix must be symmetric!  //
	//------------------------------------------------//

	void Setup(const SparseMatrix& A) override
	{
		IterativeSolver::Setup(A);
		if (!A.IsRowMajor)
			this->_rowMajorA = A;

		auto nb = A.rows() / _blockSize;
		this->diagBlockSolver.Setup(A, _blockSize);

		this->SetupComputationalWork = nb * 2.0/3.0*pow(_blockSize, 3)*1e-6;
	}

private:
	IterationResult ExecuteOneIteration(const Vector& b, Vector& xOld, bool& xEquals0, bool computeResidual, bool computeAx, const IterationResult& oldResult) override
	{
		IterationResult result(oldResult);
		assert(!computeResidual && !computeAx);

		const SparseMatrix& A = *this->Matrix;

		auto nb = A.rows() / _blockSize;

		Vector xNew(xOld.rows());

		#pragma omp parallel
		{
			Vector tmp_x(_blockSize); // work vector, allocated once per thread
			#pragma omp for
			for (BigNumber i = 0; i < nb; i++)
				ProcessBlockRow(i, b, xOld, xNew, tmp_x);
		}

		xOld.swap(xNew);
		result.SetX(xOld);
		result.AddWorkInFlops(2 * A.nonZeros() + nb * pow(_blockSize, 2));
		return result;
	}

protected:
	inline void ProcessBlockRow(BigNumber currentBlockRow, const Vector& b, const Vector& xOld, Vector& xNew, Vector& tmp_x)
	{
		const SparseMatrix& A = *this->Matrix;

		// BlockRow i: [ --- Li --- | Di | --- Ui --- ]

		tmp_x = _omega * b.segment(currentBlockRow * _blockSize, _blockSize);

		for (int k = 0; k < _blockSize; k++)
		{
			BigNumber iBlock = currentBlockRow;
			BigNumber i = iBlock * _blockSize + k;
			// RowMajor --> the following line iterates over the non-zeros of the i-th row.
			for (RowMajorSparseMatrix::InnerIterator it(A.IsRowMajor ? A : _rowMajorA, i); it; ++it)
			{
				auto j = it.col();
				auto jBlock = j / this->_blockSize;
				auto a_ij = it.value();
				if (iBlock == jBlock) // Di
					tmp_x(k) += (1 - _omega) * a_ij * xOld(j);
				else // Li and Ui
					tmp_x(k) += -_omega * a_ij * xOld(j);
			}
		}

		auto xNew_i = xNew.segment(currentBlockRow * _blockSize, _blockSize);
		this->diagBlockSolver.Solve(currentBlockRow, tmp_x, xNew_i);
	}

public:
	SparseMatrix IterationMatrix()
	{
		const SparseMatrix& A = *this->Matrix;
		// J = I - omega * D^-1 * A
		NonZeroCoefficients coeffs(A.nonZeros());

		DenseMatrix oneMinusOmegaIdentity = (1 - _omega)*DenseMatrix::Identity(_blockSize, _blockSize);

		for (int iBlock = 0; iBlock < A.rows() / _blockSize; iBlock++)
		{
			vector<BigNumber> jBlocksTreated;
			BigNumber rowBlock = iBlock * _blockSize;

			for (int k = 0; k < _blockSize; k++)
			{
				BigNumber row = iBlock * _blockSize + k;
				for (RowMajorSparseMatrix::InnerIterator it(A, row); it; ++it)
				{
					auto j = it.col();
					auto jBlock = j / this->_blockSize;

					if (find(jBlocksTreated.begin(), jBlocksTreated.end(), jBlock) != jBlocksTreated.end())
						continue;

					if (iBlock == jBlock) // Di
						coeffs.Add(rowBlock, rowBlock, oneMinusOmegaIdentity);
					else
					{
						BigNumber colBlock = jBlock * _blockSize;
						DenseMatrix A_ij = A.block(rowBlock, colBlock, _blockSize, _blockSize);
						DenseMatrix jacobiBlock(_blockSize, _blockSize);
						this->diagBlockSolver.Solve(iBlock, A_ij, jacobiBlock);
						jacobiBlock *= -_omega;
						coeffs.Add(rowBlock, colBlock, jacobiBlock);
					}

					jBlocksTreated.push_back(jBlock);
				}
			}
		}
		SparseMatrix J(A.rows(), A.cols());
		coeffs.Fill(J);
		return J;
	}
};