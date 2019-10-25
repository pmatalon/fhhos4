#pragma once
#include <Eigen/Dense>
#include "../FunctionalBasis/FunctionalBasis.h"

template <int Dim>
class ReferenceElement
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
	ReferenceElement() {}

	virtual double ComputeIntegral(RefFunction func) const = 0;
	virtual double ComputeIntegral(RefFunction func, int polynomialDegree) const = 0;

	virtual double ComputeIntegral(BasisFunction<Dim>* phi) const
	{
		RefFunction func = [phi](RefPoint p) {
			return phi->Eval(p);
		};
		return ComputeIntegral(func, phi->GetDegree());
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
		return ComputeAndReturnMassMatrix(basis);
	}
	DenseMatrix CellReconstructMassMatrix(FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim>* reconstructBasis)
	{
		if (cellBasis == _cellReconstructMassMatrixCellBasis && reconstructBasis == _cellReconstructMassMatrixReconstructBasis)
			return _cellReconstructMassMatrix;
		return ComputeAndReturnMassMatrix(cellBasis, reconstructBasis);
	}
	DenseMatrix ReconstructMassMatrix(FunctionalBasis<Dim>* basis)
	{
		if (basis == _reconstructMassMatrixBasis)
			return _reconstructMassMatrix;
		return ComputeAndReturnMassMatrix(basis);
	}
	DenseMatrix FaceMassMatrix(FunctionalBasis<Dim>* basis)
	{
		if (basis == _faceMassMatrixBasis)
			return _faceMassMatrix;
		return ComputeAndReturnMassMatrix(basis);
	}


	void ComputeAndStoreCellMassMatrix(FunctionalBasis<Dim>* basis)
	{
		if (_cellMassMatrix.rows() == 0)
		{
			_cellMassMatrix = ComputeAndReturnMassMatrix(basis);
			_cellMassMatrixBasis = basis;
		}
	}
	void ComputeAndStoreReconstructMassMatrix(FunctionalBasis<Dim>* basis)
	{
		if (_reconstructMassMatrix.rows() == 0)
		{
			_reconstructMassMatrix = ComputeAndReturnMassMatrix(basis);
			_reconstructMassMatrixBasis = basis;
		}
	}
	void ComputeAndStoreCellReconstructMassMatrix(FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim>* reconstructBasis)
	{
		if (_cellReconstructMassMatrix.rows() == 0)
		{
			_cellReconstructMassMatrix = ComputeAndReturnMassMatrix(cellBasis, reconstructBasis);
			_cellReconstructMassMatrixCellBasis = cellBasis;
			_cellReconstructMassMatrixReconstructBasis = reconstructBasis;
		}
	}
	void ComputeAndStoreFaceMassMatrix(FunctionalBasis<Dim>* basis)
	{
		if (_faceMassMatrix.rows() == 0)
		{
			_faceMassMatrix = ComputeAndReturnMassMatrix(basis);
			_faceMassMatrixBasis = basis;
		}
	}


protected:
	DenseMatrix ComputeAndReturnMassMatrix(FunctionalBasis<Dim>* basis)
	{
		DenseMatrix M = DenseMatrix(basis->Size(), basis->Size());
		for (BasisFunction<Dim>* phi1 : basis->LocalFunctions)
		{
			for (BasisFunction<Dim>* phi2 : basis->LocalFunctions)
			{
				if (phi2->LocalNumber > phi1->LocalNumber)
					break;
				double term = ComputeMassTerm(phi1, phi2);
				M(phi1->LocalNumber, phi2->LocalNumber) = term;
				M(phi2->LocalNumber, phi1->LocalNumber) = term;
			}
		}
		return M;
	}

	DenseMatrix ComputeAndReturnMassMatrix(FunctionalBasis<Dim>* basis1, FunctionalBasis<Dim>* basis2)
	{
		DenseMatrix M(basis1->LocalFunctions.size(), basis2->LocalFunctions.size());
		for (BasisFunction<Dim>* phi1 : basis1->LocalFunctions)
		{
			for (BasisFunction<Dim>* phi2 : basis2->LocalFunctions)
			{
				double term = ComputeMassTerm(phi1, phi2);
				M(phi1->LocalNumber, phi2->LocalNumber) = term;
			}
		}
		return M;
	}

public:
	double ComputeMassTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		RefFunction functionToIntegrate = [phi1, phi2](RefPoint p) {
			return phi1->Eval(p)*phi2->Eval(p);
		};

		int polynomialDegree = phi1->GetDegree() + phi2->GetDegree();
		return ComputeIntegral(functionToIntegrate, polynomialDegree);
	}

};