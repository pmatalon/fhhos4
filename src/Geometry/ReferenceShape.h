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

	double MassTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) const
	{
		return this->_massMatrix1(phi1->LocalNumber, phi2->LocalNumber);
	}

	//---------//
	//   HHO   //
	//---------//

	const DenseMatrix& StoredMassMatrix(FunctionalBasis<Dim>* basis) const
	{
		if (basis == _massMatrix1Basis)
			return _massMatrix1;
		else if (basis == _massMatrix2Basis)
			return _massMatrix2;
		assert(false && "The mass matrix for this basis should be computed once and stored for this reference element.");
	}

	const DenseMatrix& StoredCellReconstructMassMatrix() const
	{
		assert(_cellReconstructMassMatrix.rows() > 0 && "The cell-reconstruct mass matrix should be computed once and stored for this reference element.");
		return _cellReconstructMassMatrix;
	}

	void ComputeAndStoreMassMatrix(FunctionalBasis<Dim>* basis)
	{
		if (!_massMatrix1Basis)
		{
			_massMatrix1 = this->ComputeAndReturnMassMatrix(basis);
			_massMatrix1Basis = basis;
		}
		else if (!_massMatrix2Basis)
		{
			_massMatrix2 = this->ComputeAndReturnMassMatrix(basis);
			_massMatrix2Basis = basis;
		}
		else
			assert(false && "Only 2 mass matrices can be stored.");
	}
	
	void ComputeAndStoreCellReconstructMassMatrix(FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim>* reconstructBasis)
	{
		assert(_cellReconstructMassMatrix.rows() == 0 && "The cellReconstructMassMatrix is already computed.");
		_cellReconstructMassMatrix = this->ComputeAndReturnMassMatrix(cellBasis, reconstructBasis);
	}

};