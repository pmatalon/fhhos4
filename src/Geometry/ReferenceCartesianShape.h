#pragma once
#include "ReferenceShape.h"
#include "../Utils/GaussLegendre.h"

template <int Dim>
class ReferenceCartesianShape : public ReferenceShape<Dim>
{
private:
	// For DG
	DenseMatrix _stiffnessMatrix;

	// For HHO
	DenseMatrix _cellStiffnessMatrix;

	DenseMatrix _reconstructK1StiffnessMatrix;
	Tensor<Dim>* _K1;
	DenseMatrix _reconstructK2StiffnessMatrix;
	Tensor<Dim>* _K2;

public:
	ReferenceCartesianShape() : ReferenceShape<Dim>() {}

	virtual double Diameter() const override
	{
		return 2;
	}
	inline double Measure() const override
	{
		return pow(2, Dim);
	}
	virtual DomPoint Center() const override
	{
		return DomPoint(0);
	}

	vector<RefPoint> QuadraturePoints() const
	{
		GaussLegendre* gs = GaussLegendre::Get(GaussLegendre::MAX_POINTS);
		return gs->QuadraturePoints<Dim>();
	}

	double Integral(RefFunction func, int polynomialDegree) const override
	{
		if (Dim == 0)
			return func(0);
		int nPoints = GaussLegendre::NumberOfRequiredQuadraturePoint(polynomialDegree);
		return Integral(nPoints, func);
	}
	double Integral(RefFunction func) const override
	{
		return Integral(GaussLegendre::MAX_POINTS, func);
	}
	double Integral(int nPoints, RefFunction func) const
	{
		if (Dim == 0)
			return func(0);
		GaussLegendre* gs = GaussLegendre::Get(nPoints);
		return gs->QuadratureDim<Dim>(func);
	}

	//--------//
	//   DG   //
	//--------//

	double StiffnessTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		return this->_stiffnessMatrix(phi1->LocalNumber, phi2->LocalNumber);
	}

	void ComputeAndStoreStiffnessMatrix(FunctionalBasis<Dim>* basis)
	{
		if (_stiffnessMatrix.rows() == 0)
			_stiffnessMatrix = ComputeAndReturnStiffnessMatrix(basis);
	}

	//---------//
	//   HHO   //
	//---------//

	double ReconstructKStiffnessTerm(Tensor<Dim>* K, BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		if (K == _K1)
			return this->_reconstructK1StiffnessMatrix(phi1->LocalNumber, phi2->LocalNumber);
		else if (K == _K2)
			return this->_reconstructK2StiffnessMatrix(phi1->LocalNumber, phi2->LocalNumber);
		assert(false);
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

private:
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

	double ComputeIntegralGradGrad(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		if (phi1->GetDegree() == 0 || phi2->GetDegree() == 0)
			return 0;

		RefFunction functionToIntegrate = [phi1, phi2](RefPoint p) {
			return phi1->Grad(p).dot(phi2->Grad(p));
		};

		int polynomialDegree = max(0, phi1->GetDegree() + phi2->GetDegree() - 2);
		return Integral(functionToIntegrate, polynomialDegree);
	}

	double ComputeIntegralKGradGrad(Tensor<Dim>* K, BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		if (phi1->GetDegree() == 0 || phi2->GetDegree() == 0)
			return 0;

		RefFunction functionToIntegrate = [K, phi1, phi2](RefPoint p) {
			return (K * phi1->Grad(p)).dot(phi2->Grad(p));
		};

		int polynomialDegree = max(0, phi1->GetDegree() + phi2->GetDegree() - 2);
		return Integral(functionToIntegrate, polynomialDegree);
	}
};