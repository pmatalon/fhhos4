#pragma once
#include "ReferenceShape.h"
#include "../QuadratureRules/GaussLegendre/GaussLegendre.h"
#include "../TestCases/Diffusion/Tensor.h"

template <int Dim>
class ReferenceCartesianShape : public ReferenceShape<Dim>
{
private:
	// For DG
	DenseMatrix _stiffnessMatrix;

	// For HHO
	map<const Tensor<Dim>*, DenseMatrix> _reconstructStiffnessMatrices;

public:
	ReferenceCartesianShape() : ReferenceShape<Dim>() {}

	string Name() const override { return "Reference CartesianShape<" + to_string(Dim) + ">"; }

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

	vector<RefPoint> QuadraturePoints() const override
	{
		GaussLegendre* gs = GaussLegendre::Get(GaussLegendre::MAX_POINTS);
		return gs->QuadraturePoints<Dim>();
	}

	vector<RefPoint> QuadraturePoints(int polynomialDegree) const override
	{
		GaussLegendre* gs = GaussLegendre::Get(polynomialDegree);
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

	double ReconstructStiffnessTerm(const Tensor<Dim>& K, BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		auto it = _reconstructStiffnessMatrices.find(&K);
		if (it != _reconstructStiffnessMatrices.end())
			return it->second(phi1->LocalNumber, phi2->LocalNumber);
		else
			return this->ComputeIntegralKGradGrad(K, phi1, phi2);
	}

	void ComputeAndStoreReconstructStiffnessMatrix(const Tensor<Dim>& K, FunctionalBasis<Dim>* basis)
	{
		if (_reconstructStiffnessMatrices.find(&K) == _reconstructStiffnessMatrices.end())
			_reconstructStiffnessMatrices[&K] = ComputeAndReturnKStiffnessMatrix(K, basis);
	}

private:
	DenseMatrix ComputeAndReturnStiffnessMatrix(FunctionalBasis<Dim>* basis)
	{
		DenseMatrix stiffnessMatrix = DenseMatrix(basis->Size(), basis->Size());
		for (BasisFunction<Dim>* phi1 : basis->LocalFunctions())
		{
			for (BasisFunction<Dim>* phi2 : basis->LocalFunctions())
				stiffnessMatrix(phi1->LocalNumber, phi2->LocalNumber) = ComputeIntegralGradGrad(phi1, phi2);
		}
		return stiffnessMatrix;
	}

	DenseMatrix ComputeAndReturnKStiffnessMatrix(const Tensor<Dim>& K, FunctionalBasis<Dim>* basis)
	{
		DenseMatrix stiffnessMatrix = DenseMatrix(basis->Size(), basis->Size());
		for (BasisFunction<Dim>* phi1 : basis->LocalFunctions())
		{
			for (BasisFunction<Dim>* phi2 : basis->LocalFunctions())
				stiffnessMatrix(phi1->LocalNumber, phi2->LocalNumber) = ComputeIntegralKGradGrad(K, phi1, phi2);
		}
		return stiffnessMatrix;
	}

public:

	double ComputeIntegralGradGrad(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		if (phi1->GetDegree() == 0 || phi2->GetDegree() == 0)
			return 0;

		RefFunction functionToIntegrate = [phi1, phi2](const RefPoint& p) {
			return phi1->Grad(p).dot(phi2->Grad(p));
		};

		int polynomialDegree = max(0, phi1->GetDegree() + phi2->GetDegree() - 2);
		return Integral(functionToIntegrate, polynomialDegree);
	}

	double ComputeIntegralKGradGrad(const Tensor<Dim>& K, BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		if (phi1->GetDegree() == 0 || phi2->GetDegree() == 0)
			return 0;

		RefFunction functionToIntegrate = [&K, phi1, phi2](const RefPoint& p) {
			return (K * phi1->Grad(p)).dot(phi2->Grad(p));
		};

		int polynomialDegree = max(0, phi1->GetDegree() + phi2->GetDegree() - 2);
		return Integral(functionToIntegrate, polynomialDegree);
	}
};