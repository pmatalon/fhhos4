#pragma once
#include "GeometricShape.h"

template <int Dim>
class ReferenceShape : public GeometricShape<Dim>
{
protected:
	// For DG
	DenseMatrix _massMatrix;

	// For HHO
	DenseMatrix _cellMassMatrix;
	FunctionalBasis<Dim>* _cellMassMatrixBasis;

	DenseMatrix _reconstructMassMatrix;
	FunctionalBasis<Dim>* _reconstructMassMatrixBasis;

	DenseMatrix _cellReconstructMassMatrix;
	FunctionalBasis<Dim>* _cellReconstructMassMatrixCellBasis;
	FunctionalBasis<Dim>* _cellReconstructMassMatrixReconstructBasis;

	DenseMatrix _faceMassMatrix;
	FunctionalBasis<Dim>* _faceMassMatrixBasis;

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
		return this->_massMatrix(phi1->LocalNumber, phi2->LocalNumber);
	}

	void ComputeAndStoreMassMatrix(FunctionalBasis<Dim>* basis)
	{
		if (_massMatrix.rows() == 0)
			_massMatrix = this->ComputeAndReturnMassMatrix(basis);
	}

	//---------//
	//   HHO   //
	//---------//

	DenseMatrix CellMassMatrix(FunctionalBasis<Dim>* basis)
	{
		if (basis == _cellMassMatrixBasis)
			return _cellMassMatrix;
		return this->ComputeAndReturnMassMatrix(basis);
	}
	DenseMatrix CellReconstructMassMatrix(FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim>* reconstructBasis)
	{
		if (cellBasis == _cellReconstructMassMatrixCellBasis && reconstructBasis == _cellReconstructMassMatrixReconstructBasis)
			return _cellReconstructMassMatrix;
		return this->ComputeAndReturnMassMatrix(cellBasis, reconstructBasis);
	}
	DenseMatrix ReconstructMassMatrix(FunctionalBasis<Dim>* basis)
	{
		if (basis == _reconstructMassMatrixBasis)
			return _reconstructMassMatrix;
		return this->ComputeAndReturnMassMatrix(basis);
	}
	DenseMatrix FaceMassMatrix(FunctionalBasis<Dim>* basis)
	{
		if (basis == _faceMassMatrixBasis)
			return _faceMassMatrix;
		return this->ComputeAndReturnMassMatrix(basis);
	}


	void ComputeAndStoreCellMassMatrix(FunctionalBasis<Dim>* basis)
	{
		if (_cellMassMatrix.rows() == 0)
		{
			_cellMassMatrix = this->ComputeAndReturnMassMatrix(basis);
			_cellMassMatrixBasis = basis;
		}
	}
	void ComputeAndStoreReconstructMassMatrix(FunctionalBasis<Dim>* basis)
	{
		if (_reconstructMassMatrix.rows() == 0)
		{
			_reconstructMassMatrix = this->ComputeAndReturnMassMatrix(basis);
			_reconstructMassMatrixBasis = basis;
		}
	}
	void ComputeAndStoreCellReconstructMassMatrix(FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim>* reconstructBasis)
	{
		if (_cellReconstructMassMatrix.rows() == 0)
		{
			_cellReconstructMassMatrix = this->ComputeAndReturnMassMatrix(cellBasis, reconstructBasis);
			_cellReconstructMassMatrixCellBasis = cellBasis;
			_cellReconstructMassMatrixReconstructBasis = reconstructBasis;
		}
	}
	void ComputeAndStoreFaceMassMatrix(FunctionalBasis<Dim>* basis)
	{
		if (_faceMassMatrix.rows() == 0)
		{
			_faceMassMatrix = this->ComputeAndReturnMassMatrix(basis);
			_faceMassMatrixBasis = basis;
		}
	}

};