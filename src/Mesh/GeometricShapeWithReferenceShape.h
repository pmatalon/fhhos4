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
	virtual DimMatrix<Dim> InverseJacobianTranspose() const = 0;
	virtual double DetJacobian() const = 0;

	//-------------------//
	//     Integrals     //
	//-------------------//

	double Integral(BasisFunction<Dim>* phi) const
	{
		return GeometricShape<Dim>::Integral(phi);
	}

	double Integral(RefFunction func) const override
	{
		return DetJacobian() * RefShape()->Integral(func);
	}

	double Integral(RefFunction func, int polynomialDegree) const override
	{
		return DetJacobian() * RefShape()->Integral(func, polynomialDegree);
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

	double ComputeIntegralGradGrad(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) const
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

	double ComputeIntegralKGradGrad(Tensor<Dim>* K, BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) const
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

	inline double MassTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		return DetJacobian() * RefShape()->MassTerm(phi1, phi2);
	}

	//-----------------------------//
	//             HHO             //
	//-----------------------------//

	inline DenseMatrix FaceMassMatrix(FunctionalBasis<Dim>* basis)
	{
		return DetJacobian() * RefShape()->FaceMassMatrix(basis);
	}

	inline DenseMatrix CellMassMatrix(FunctionalBasis<Dim>* basis)
	{
		return DetJacobian() * RefShape()->CellMassMatrix(basis);
	}

	inline DenseMatrix CellReconstructMassMatrix(FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim>* reconstructBasis)
	{
		return DetJacobian() * RefShape()->CellReconstructMassMatrix(cellBasis, reconstructBasis);
	}
};