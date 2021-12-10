#pragma once
#include "AggregAMG.h"
using namespace std;

class HighOrderAggregLevel : public Level
{
private:
	int _blockSize;
	double _omega;
private:
	SparseMatrix Q_T;
	SparseMatrix Q_F;

public:
	HighOrderAggregLevel(int number, int blockSize, double omega)
		: Level(number)
	{
		this->_blockSize = blockSize;
		this->_omega = omega;
	}

	int PolynomialDegree() override
	{
		return 0; // TODO for p-multigrid
	}

	void CoarsenMesh(H_CoarsStgy coarseningStgy, FaceCoarseningStrategy faceCoarseningStgy, double coarseningFactor, bool& noCoarserMeshProvided, bool& coarsestPossibleMeshReached) override
	{
		const SparseMatrix& A = *this->OperatorMatrix;

		BlockJacobi blockJacobi(_blockSize, _omega);
		blockJacobi.Setup(A);
		SparseMatrix J = blockJacobi.IterationMatrix();
		//SparseMatrix J = BlockJacobiMatrix(A, _omega);

		NonZeroCoefficients coeffs(A.rows() / _blockSize);
		for (BigNumber i = 0; i < A.rows(); i+=_blockSize)
			coeffs.Add(i, i/_blockSize, 1);
		SparseMatrix P0(A.rows(), A.rows() / _blockSize);
		coeffs.Fill(P0);

		this->P = J * P0;
	}

/*private:
	SparseMatrix BlockJacobiMatrix(const SparseMatrix& A, double omega)
	{
		NonZeroCoefficients coeffs(A.nonZeros());

		for (int iBlock = 0; iBlock < A.rows()/_blockSize; iBlock++)
		{
			DenseMatrix A_ii = A.block(iBlock, iBlock, _blockSize, _blockSize);
			for (int k = 0; k < _blockSize; k++)
			{
				BigNumber i = iBlock * _blockSize + k;
				for (RowMajorSparseMatrix::InnerIterator it(A, i); it; ++it)
				{
					if (it.row() == it.col())
						coeffs.Add(it.row(), it.col(), 1 - omega);
					else
						coeffs.Add(it.row(), it.col(), -omega * it.value() / a_ii);
				}
			}
		}
		SparseMatrix J(A.rows(), A.cols());
		coeffs.Fill(J);
		return J;
	}*/

public:
	void OnStartSetup() override
	{
		cout << "\t\tMesh                : " << this->NUnknowns() / _blockSize << " elements";
		if (!this->IsFinestLevel())
		{
			double nFine = this->FinerLevel->NUnknowns();
			double nCoarse = this->NUnknowns();
			cout << ", coarsening factor = " << (nFine/nCoarse);
		}
		cout << endl;
	}

	void SetupProlongation() override
	{
	}

	void SetupRestriction() override
	{
		double scalingFactor = 1.0;
		R = (scalingFactor * P.transpose()).eval();
	}
};

class HighOrderAggregAMG : public Multigrid
{
private:
	int _blockSize;
	double _strongCouplingThreshold;
	double _omega;
public:

	HighOrderAggregAMG(int blockSize, double strongCouplingThreshold, double omega, int nLevels = 0)
		: Multigrid(nLevels)
	{
		this->_blockSize = blockSize;
		this->_strongCouplingThreshold = strongCouplingThreshold;
		this->_omega = omega;
		this->BlockSizeForBlockSmoothers = blockSize;
		this->UseGalerkinOperator = true;
	}

	void BeginSerialize(ostream& os) const override
	{
		os << "HighOrderAggregAMG" << endl;
	}

	void EndSerialize(ostream& os) const override
	{
	}

	Vector Solve(const Vector& b, string initialGuessCode) override
	{
		if (initialGuessCode.compare("smooth") == 0)
			Utils::Warning("Smooth initial guess unmanaged in AMG.");
		return Multigrid::Solve(b, initialGuessCode);
	}
	
protected:
	Level* CreateFineLevel() const override
	{
		return new HighOrderAggregLevel(0, _blockSize, _omega);
	}

	Level* CreateCoarseLevel(Level* fineLevel, CoarseningType coarseningType, int coarseDegree) override
	{
		//AggregLevel* coarseLevel = new AggregLevel(fineLevel->Number + 1, 1, _strongCouplingThreshold);
		//return coarseLevel;
		HighOrderAggregLevel* coarseLevel = new HighOrderAggregLevel(fineLevel->Number + 1, 1, 0);
		return coarseLevel;
	}
};