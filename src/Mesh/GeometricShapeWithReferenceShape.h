#pragma once
#include "ReferenceShape.h"
#include "../Problem/Tensor.h"

template <int Dim>
class GeometricShapeWithReferenceShape : public GeometricShape<Dim>
{
public:
	GeometricShapeWithReferenceShape() : GeometricShape<Dim>() {}

	//-----------------------//
	//   Virtual functions   //
	//-----------------------//

	virtual ReferenceShape<Dim>* RefShape() const = 0;

	// Transformation to reference element
	virtual DomPoint ConvertToDomain(RefPoint refPoint) const = 0;
	virtual RefPoint ConvertToReference(DomPoint domainPoint) const = 0;

	virtual DimMatrix<Dim> InverseJacobianTranspose(RefPoint p) const = 0;
	virtual double DetJacobian(RefPoint p) const = 0;
	virtual int DetJacobianDegree() const = 0;

	//-------------------//
	//     Integrals     //
	//-------------------//

	double Integral(BasisFunction<Dim>* phi) const
	{
		return GeometricShape<Dim>::Integral(phi);
	}

	virtual double Integral(RefFunction f) const override
	{
		RefFunction func = [this, f](RefPoint p) {
			return DetJacobian(p) * f(p);
		};
		return RefShape()->Integral(func);
	}

	virtual double Integral(RefFunction f, int polynomialDegree) const override
	{
		RefFunction func = [this, f](RefPoint p) {
			return DetJacobian(p) * f(p);
		};
		return RefShape()->Integral(func, polynomialDegree + DetJacobianDegree());
	}

	virtual double Integral(DomFunction globalFunction) const
	{
		RefFunction refFunction = [this, globalFunction](RefPoint refElementPoint) {
			DomPoint domainPoint = this->ConvertToDomain(refElementPoint);
			return globalFunction(domainPoint);
		};

		return Integral(refFunction);
	}

	virtual double Integral(DomFunction globalFunction, int polynomialDegree) const
	{
		RefFunction refFunction = [this, globalFunction](RefPoint refElementPoint) {
			DomPoint domainPoint = this->ConvertToDomain(refElementPoint);
			return globalFunction(domainPoint);
		};

		return Integral(refFunction, polynomialDegree);
	}

	//----------------------------//
	//     Specific integrals     //
	//----------------------------//

	virtual double ComputeIntegralGradGrad(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) const
	{
		if (phi1->GetDegree() == 0 || phi2->GetDegree() == 0)
			return 0;

		RefFunction functionToIntegrate = [this, phi1, phi2](RefPoint p) {
			DimMatrix<Dim> invJ = InverseJacobianTranspose(p);
			DimVector<Dim> gradPhi1 = invJ * phi1->Grad(p);
			DimVector<Dim> gradPhi2 = invJ * phi2->Grad(p);
			return gradPhi1.dot(gradPhi2);
		};

		int polynomialDegree = max(0, phi1->GetDegree() + phi2->GetDegree() - 2);
		return Integral(functionToIntegrate, polynomialDegree);
	}

	virtual double ComputeIntegralKGradGrad(Tensor<Dim>* K, BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) const
	{
		if (phi1->GetDegree() == 0 || phi2->GetDegree() == 0)
			return 0;

		RefFunction functionToIntegrate = [this, K, phi1, phi2](RefPoint p) {
			DimMatrix<Dim> invJ = InverseJacobianTranspose(p);
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

	virtual double MassTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		return this->ComputeMassTerm(phi1, phi2);
	}

	//-----------------------------//
	//             HHO             //
	//-----------------------------//

	virtual DenseMatrix FaceMassMatrix(FunctionalBasis<Dim>* basis)
	{
		return this->ComputeAndReturnMassMatrix(basis);
	}

	virtual DenseMatrix CellMassMatrix(FunctionalBasis<Dim>* basis)
	{
		return this->ComputeAndReturnMassMatrix(basis);
	}

	virtual DenseMatrix CellReconstructMassMatrix(FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim>* reconstructBasis)
	{
		return this->ComputeAndReturnMassMatrix(cellBasis, reconstructBasis);
	}
};