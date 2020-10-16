#pragma once
#include "GeometricShape.h"

template <int Dim>
class ReferenceShape : public GeometricShape<Dim>
{
private:
	DenseMatrix _massMatrix1;
	FunctionalBasis<Dim>* _massMatrix1Basis = nullptr;

	DenseMatrix _massMatrix2;
	FunctionalBasis<Dim>* _massMatrix2Basis = nullptr;

	DenseMatrix _cellReconstructMassMatrix;
	FunctionalBasis<Dim>* _cellReconstructMassMatrixCellBasis;
	FunctionalBasis<Dim>* _cellReconstructMassMatrixReconstructBasis;

public:
	ReferenceShape() {}

	double InRadius() const override
	{
		assert(false);
	}

	virtual vector<RefPoint> QuadraturePoints() const = 0;
	virtual vector<RefPoint> QuadraturePoints(int polynomialDegree) const = 0;

	virtual double Integral(RefFunction func) const = 0;
	virtual double Integral(RefFunction func, int polynomialDegree) const = 0;

	double Integral(DomFunction globalFunction) const override
	{
		assert(false);
	}
	double Integral(DomFunction globalFunction, int polynomialDegree) const override
	{
		assert(false);
	}

	virtual void Serialize(ostream& os) const
	{
		assert(false);
	}

	//--------//
	//   DG   //
	//--------//

	double MassTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		return this->_massMatrix1(phi1->LocalNumber, phi2->LocalNumber);
	}

	//---------//
	//   HHO   //
	//---------//

	DenseMatrix MassMatrix(FunctionalBasis<Dim>* basis)
	{
		if (basis == _massMatrix1Basis)
			return _massMatrix1;
		else if (basis == _massMatrix2Basis)
			return _massMatrix2;
		Utils::Warning("The mass matrix for this basis should be computed once and stored for this reference element.");
		return this->ComputeAndReturnMassMatrix(basis);
	}
	DenseMatrix CellReconstructMassMatrix(FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim>* reconstructBasis)
	{
		if (cellBasis == _cellReconstructMassMatrixCellBasis && reconstructBasis == _cellReconstructMassMatrixReconstructBasis)
			return _cellReconstructMassMatrix;
		return this->ComputeAndReturnMassMatrix(cellBasis, reconstructBasis);
	}
	void ComputeAndStoreMassMatrix(FunctionalBasis<Dim>* basis)
	{
		if (!_massMatrix1Basis)
			ComputeAndStoreMassMatrix1(basis);
		else if (!_massMatrix2Basis)
			ComputeAndStoreMassMatrix2(basis);
		else
			assert(false);
	}

private:
	void ComputeAndStoreMassMatrix1(FunctionalBasis<Dim>* basis)
	{
		if (_massMatrix1.rows() == 0)
		{
			_massMatrix1 = this->ComputeAndReturnMassMatrix(basis);
			_massMatrix1Basis = basis;
		}
	}
	void ComputeAndStoreMassMatrix2(FunctionalBasis<Dim>* basis)
	{
		if (_massMatrix2.rows() == 0)
		{
			_massMatrix2 = this->ComputeAndReturnMassMatrix(basis);
			_massMatrix2Basis = basis;
		}
	}

public:
	void ComputeAndStoreCellReconstructMassMatrix(FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim>* reconstructBasis)
	{
		if (_cellReconstructMassMatrix.rows() == 0)
		{
			_cellReconstructMassMatrix = this->ComputeAndReturnMassMatrix(cellBasis, reconstructBasis);
			_cellReconstructMassMatrixCellBasis = cellBasis;
			_cellReconstructMassMatrixReconstructBasis = reconstructBasis;
		}
	}

};