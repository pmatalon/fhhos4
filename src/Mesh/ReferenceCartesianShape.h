#pragma once
#include "../FunctionalBasis/FunctionalBasis.h"

template <int Dim>
class ReferenceCartesianShape
{
private:
	// For DG
	DenseMatrix _stiffnessMatrix;
	DenseMatrix _massMatrix;

	// For HHO
	DenseMatrix _cellMassMatrix;
	FunctionalBasis<Dim>* _cellMassMatrixBasis;

	DenseMatrix _reconstructMassMatrix;
	FunctionalBasis<Dim>* _reconstructMassMatrixBasis;

	DenseMatrix _cellStiffnessMatrix;

	DenseMatrix _reconstructK1StiffnessMatrix;
	Tensor<Dim>* _K1;
	DenseMatrix _reconstructK2StiffnessMatrix;
	Tensor<Dim>* _K2;

	DenseMatrix _faceMassMatrix;
	FunctionalBasis<Dim>* _faceMassMatrixBasis;

	DenseMatrix _cellReconstructMassMatrix;
	FunctionalBasis<Dim>* _cellReconstructMassMatrixCellBasis;
	FunctionalBasis<Dim>* _cellReconstructMassMatrixReconstructBasis;

public:
	ReferenceCartesianShape() {}

	double ComputeIntegral(BasisFunction<Dim>* phi)
	{
		return Utils::Integral(phi);
	}

	double ComputeIntegral(function<double(RefPoint)> func, int polynomialDegree)
	{
		return Utils::Integral<Dim>(func, polynomialDegree);
	}

	double ComputeIntegral(function<double(RefPoint)> func)
	{
		return Utils::Integral<Dim>(func);
	}

	//--------//
	//   DG   //
	//--------//

	double MassTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		return this->_massMatrix(phi1->LocalNumber, phi2->LocalNumber);
	}
	double StiffnessTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		return this->_stiffnessMatrix(phi1->LocalNumber, phi2->LocalNumber);
	}

	void ComputeAndStoreMassMatrix(FunctionalBasis<Dim>* basis)
	{
		if (_massMatrix.rows() == 0)
			_massMatrix = ComputeAndReturnMassMatrix(basis);
	}
	void ComputeAndStoreStiffnessMatrix(FunctionalBasis<Dim>* basis)
	{
		if (_stiffnessMatrix.rows() == 0)
			_stiffnessMatrix = ComputeAndReturnStiffnessMatrix(basis);
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
	DenseMatrix CellReconstructMassMatrix(FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim>* reconstructBasis)
	{
		if (cellBasis == _cellReconstructMassMatrixCellBasis && reconstructBasis == _cellReconstructMassMatrixReconstructBasis)
			return _cellReconstructMassMatrix;
		return ComputeAndReturnMassMatrix(cellBasis, reconstructBasis);
	}

	double ReconstructKStiffnessTerm(Tensor<Dim>* K, BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		if (K == _K1)
			return this->_reconstructK1StiffnessMatrix(phi1->LocalNumber, phi2->LocalNumber);
		else if (K == _K2)
			return this->_reconstructK2StiffnessMatrix(phi1->LocalNumber, phi2->LocalNumber);
		assert(false);
	}

	void ComputeAndStoreCellMassMatrix(FunctionalBasis<Dim>* basis)
	{
		if (_cellMassMatrix.rows() == 0)
		{
			_cellMassMatrix = ComputeAndReturnMassMatrix(basis);
			_cellMassMatrixBasis = basis;
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
	void ComputeAndStoreReconstructMassMatrix(FunctionalBasis<Dim>* basis)
	{
		if (_reconstructMassMatrix.rows() == 0)
		{
			_reconstructMassMatrix = ComputeAndReturnMassMatrix(basis);
			_reconstructMassMatrixBasis = basis;
		}
	}
	void ComputeAndStoreCellStiffnessMatrix(FunctionalBasis<Dim>* basis)
	{
		if (_cellStiffnessMatrix.rows() == 0)
			_cellStiffnessMatrix = ComputeAndReturnStiffnessMatrix(basis);
	}
	void ComputeAndStoreReconstructK1StiffnessMatrix(Tensor<Dim>* K, FunctionalBasis<Dim>* basis)
	{
		if (_reconstructK1StiffnessMatrix.rows() == 0)
		{
			_reconstructK1StiffnessMatrix = ComputeAndReturnKStiffnessMatrix(K, basis);
			_K1 = K;
		}
	}
	void ComputeAndStoreReconstructK2StiffnessMatrix(Tensor<Dim>* K, FunctionalBasis<Dim>* basis)
	{
		if (_reconstructK2StiffnessMatrix.rows() == 0)
		{
			_reconstructK2StiffnessMatrix = ComputeAndReturnKStiffnessMatrix(K, basis);
			_K2 = K;
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

private:
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

	DenseMatrix ComputeAndReturnStiffnessMatrix(FunctionalBasis<Dim>* basis)
	{
		DenseMatrix stiffnessMatrix = DenseMatrix(basis->Size(), basis->Size());
		for (BasisFunction<Dim>* phi1 : basis->LocalFunctions)
		{
			for (BasisFunction<Dim>* phi2 : basis->LocalFunctions)
				stiffnessMatrix(phi1->LocalNumber, phi2->LocalNumber) = ComputeIntegralGradGrad(phi1, phi2);
		}
		return stiffnessMatrix;
	}

	DenseMatrix ComputeAndReturnKStiffnessMatrix(Tensor<Dim>* K, FunctionalBasis<Dim>* basis)
	{
		DenseMatrix stiffnessMatrix = DenseMatrix(basis->Size(), basis->Size());
		for (BasisFunction<Dim>* phi1 : basis->LocalFunctions)
		{
			for (BasisFunction<Dim>* phi2 : basis->LocalFunctions)
				stiffnessMatrix(phi1->LocalNumber, phi2->LocalNumber) = ComputeIntegralKGradGrad(K, phi1, phi2);
		}
		return stiffnessMatrix;
	}
public:
	double ComputeMassTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		function<double(RefPoint)> functionToIntegrate = [phi1, phi2](RefPoint p) {
			return phi1->Eval(p)*phi2->Eval(p);
		};

		int polynomialDegree = phi1->GetDegree() + phi2->GetDegree();
		return Utils::Integral<Dim>(functionToIntegrate, polynomialDegree);
	}

	double ComputeIntegralGradGrad(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		if (phi1->GetDegree() == 0 || phi2->GetDegree() == 0)
			return 0;

		function<double(RefPoint)> functionToIntegrate = [phi1, phi2](RefPoint p) {
			return phi1->Grad(p).dot(phi2->Grad(p));
		};

		int polynomialDegree = max(0, phi1->GetDegree() + phi2->GetDegree() - 2);
		return Utils::Integral<Dim>(functionToIntegrate, polynomialDegree);
	}

	double ComputeIntegralKGradGrad(Tensor<Dim>* K, BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		if (phi1->GetDegree() == 0 || phi2->GetDegree() == 0)
			return 0;

		function<double(RefPoint)> functionToIntegrate = [K, phi1, phi2](RefPoint p) {
			return (K * phi1->Grad(p)).dot(phi2->Grad(p));
		};

		int polynomialDegree = max(0, phi1->GetDegree() + phi2->GetDegree() - 2);
		return Utils::Integral<Dim>(functionToIntegrate, polynomialDegree);
	}
};