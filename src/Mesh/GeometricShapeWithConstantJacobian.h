#pragma once
#include "GeometricShapeWithReferenceShape.h"

template <int Dim>
class GeometricShapeWithConstantJacobian : public GeometricShapeWithReferenceShape<Dim>
{
public:
	GeometricShapeWithConstantJacobian() : GeometricShapeWithReferenceShape<Dim>() {}

	//-----------------------//
	//   Virtual functions   //
	//-----------------------//

	virtual DimMatrix<Dim> InverseJacobianTranspose() const = 0;
	virtual double DetJacobian() const = 0;

	//-----------------------//

	inline DimMatrix<Dim> InverseJacobianTranspose(RefPoint p) const override
	{
		return InverseJacobianTranspose();
	}
	inline double DetJacobian(RefPoint p) const override
	{
		return DetJacobian();
	}
	inline int DetJacobianDegree() const override
	{
		return 0;
	}

	//-------------------//
	//     Integrals     //
	//-------------------//

	inline double Integral(RefFunction func) const override
	{
		return DetJacobian() * this->RefShape()->Integral(func);
	}

	inline double Integral(RefFunction func, int polynomialDegree) const override
	{
		return DetJacobian() * this->RefShape()->Integral(func, polynomialDegree);
	}

	//----------------------------//
	//     Specific integrals     //
	//----------------------------//

	double ComputeIntegralGradGrad(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) const override
	{
		if (phi1->GetDegree() == 0 || phi2->GetDegree() == 0)
			return 0;

		DimMatrix<Dim> invJ = InverseJacobianTranspose();

		RefFunction functionToIntegrate = [phi1, phi2, invJ](RefPoint p) {
			DimVector<Dim> gradPhi1 = invJ * phi1->Grad(p);
			DimVector<Dim> gradPhi2 = invJ * phi2->Grad(p);
			return gradPhi1.dot(gradPhi2);
		};

		int polynomialDegree = max(0, phi1->GetDegree() + phi2->GetDegree() - 2);
		return Integral(functionToIntegrate, polynomialDegree);
	}

	double ComputeIntegralKGradGrad(Tensor<Dim>* K, BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) const override
	{
		if (phi1->GetDegree() == 0 || phi2->GetDegree() == 0)
			return 0;

		DimMatrix<Dim> invJ = InverseJacobianTranspose();

		RefFunction functionToIntegrate = [K, phi1, phi2, invJ](RefPoint p) {
			DimVector<Dim> gradPhi1 = invJ * phi1->Grad(p);
			DimVector<Dim> gradPhi2 = invJ * phi2->Grad(p);
			return (K * gradPhi1).dot(gradPhi2);
		};

		int polynomialDegree = max(0, phi1->GetDegree() + phi2->GetDegree() - 2);
		return Integral(functionToIntegrate, polynomialDegree);
	}

	//----------------------------//
	//             DG             //
	//----------------------------//

	inline double MassTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) override
	{
		return DetJacobian() * this->RefShape()->MassTerm(phi1, phi2);
	}

	//-----------------------------//
	//             HHO             //
	//-----------------------------//

	inline DenseMatrix FaceMassMatrix(FunctionalBasis<Dim>* basis) override
	{
		return DetJacobian() * this->RefShape()->FaceMassMatrix(basis);
	}

	inline DenseMatrix CellMassMatrix(FunctionalBasis<Dim>* basis) override
	{
		return DetJacobian() * this->RefShape()->CellMassMatrix(basis);
	}

	inline DenseMatrix CellReconstructMassMatrix(FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim>* reconstructBasis) override
	{
		return DetJacobian() * this->RefShape()->CellReconstructMassMatrix(cellBasis, reconstructBasis);
	}
};