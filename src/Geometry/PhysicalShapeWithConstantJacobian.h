#pragma once
#include "PhysicalShape.h"
#include "../FunctionalBasis/Orthogonal/OrthogonalBasisOnCstJacShape.h"

template <int Dim>
class PhysicalShapeWithConstantJacobian : public PhysicalShape<Dim>
{
public:
	PhysicalShapeWithConstantJacobian() : PhysicalShape<Dim>() {}

	//-----------------------//
	//   Virtual functions   //
	//-----------------------//

	virtual DimMatrix<Dim> InverseJacobianTranspose() const = 0;
	virtual double DetJacobian() const = 0;

	//-----------------------//

	inline DimMatrix<Dim> InverseJacobianTranspose(const RefPoint& p) const override
	{
		return InverseJacobianTranspose(); // does not depend on p
	}
	inline double DetJacobian(const RefPoint& p) const override
	{
		return DetJacobian(); // does not depend on p
	}
	inline int DetJacobianDegree() const override
	{
		return 0; // DetJacobian is a constant, so degree 0
	}

	//-------------------//
	//     Integrals     //
	//-------------------//
	// As the Jacobian is constant, DetJacobian() can be put out of the integral.

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
	// As the Jacobian is independent of p, InverseJacobianTranspose() can be retrieved once out of the integral.
	// But that's a very small optimization, not sure it fastens anything...

	double ComputeIntegralGradGrad(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) const override
	{
		if (phi1->GetDegree() == 0 || phi2->GetDegree() == 0)
			return 0;

		DimMatrix<Dim> invJ = InverseJacobianTranspose();

		RefFunction functionToIntegrate = [phi1, phi2, &invJ](const RefPoint& p) {
			DimVector<Dim> gradPhi1 = invJ * phi1->Grad(p);
			DimVector<Dim> gradPhi2 = invJ * phi2->Grad(p);
			return gradPhi1.dot(gradPhi2);
		};

		int polynomialDegree = max(0, phi1->GetDegree() + phi2->GetDegree() - 2);
		return Integral(functionToIntegrate, polynomialDegree);
	}

	double ComputeIntegralKGradGrad(const Tensor<Dim>& K, BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) const override
	{
		if (phi1->GetDegree() == 0 || phi2->GetDegree() == 0)
			return 0;

		DimMatrix<Dim> invJ = InverseJacobianTranspose();

		RefFunction functionToIntegrate = [&K, phi1, phi2, &invJ](const RefPoint& p) {
			DimVector<Dim> gradPhi1 = invJ * phi1->Grad(p);
			DimVector<Dim> gradPhi2 = invJ * phi2->Grad(p);
			return (K * gradPhi1).dot(gradPhi2);
		};

		int polynomialDegree = max(0, phi1->GetDegree() + phi2->GetDegree() - 2);
		return Integral(functionToIntegrate, polynomialDegree);
	}

	// See http://arturo.imati.cnr.it/~marini/didattica/Metodi-engl/Intro2FEM.pdf (page 30)
	DenseMatrix IntegralKGradGradMatrix(const Tensor<Dim>& K, FunctionalBasis<Dim>* basis) const override
	{
		FunctionalBasis<Dim>* refShapeBasis = basis;
		OrthogonalBasisOnCstJacShape<Dim>* orthogBasis = dynamic_cast<OrthogonalBasisOnCstJacShape<Dim>*>(basis);
		if (orthogBasis)
			refShapeBasis = orthogBasis->RefShapeBasis;

		double coeff = DetJacobian();
		if (orthogBasis && orthogBasis->IsNormalized())
			coeff = 1; // = coeff / sqrt(DetJacobian())^2

		DimMatrix<Dim> C = coeff * InverseJacobianTranspose().transpose() * K.TensorMatrix * InverseJacobianTranspose();
		const StiffnessMatrices& refStiff = this->RefShape()->StoredStiffnessMatrices(refShapeBasis);
		int t = 0;
		int u = 1;
		int v = 2;
		if (Dim == 1)
			return C(t, t)*refStiff.tt;
		else if (Dim == 2)
			return C(t, t)*refStiff.tt + C(u, u)*refStiff.uu + C(t, u)*(refStiff.tu + refStiff.tu.transpose());
		else if (Dim == 3)
			return C(t, t)*refStiff.tt + C(u, u)*refStiff.uu + C(v, v)*refStiff.vv + C(t, u)*(refStiff.tu + refStiff.tu.transpose()) + C(t, v)*(refStiff.tv + refStiff.tv.transpose()) + C(u, v)*(refStiff.uv + refStiff.uv.transpose());
		assert(false);
		return DenseMatrix();
	}


	//----------------------------//
	//             DG             //
	//----------------------------//

	inline double MassTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) const override
	{
		return DetJacobian() * this->RefShape()->MassTerm(phi1, phi2);
	}

	//-----------------------------//
	//             HHO             //
	//-----------------------------//

	DenseMatrix MassMatrix(FunctionalBasis<Dim>* basis) const override
	{
		return DetJacobian() * this->RefShape()->StoredMassMatrix(basis);
	}

	DenseMatrix CellReconstructMassMatrix(FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim>* reconstructBasis) const override
	{
		return DetJacobian() * this->RefShape()->StoredCellReconstructMassMatrix(cellBasis, reconstructBasis);
	}

	Vector Integral(FunctionalBasis<Dim>* basis) const override
	{
		FunctionalBasis<Dim>* refShapeBasis = basis;
		OrthogonalBasisOnCstJacShape<Dim>* orthogBasis = dynamic_cast<OrthogonalBasisOnCstJacShape<Dim>*>(basis);
		if (orthogBasis)
			refShapeBasis = orthogBasis->RefShapeBasis;

		double coeff = DetJacobian();
		if (orthogBasis && orthogBasis->IsNormalized())
			coeff = sqrt(DetJacobian()); // coeff / sqrt(DetJacobian())
		return coeff * this->RefShape()->StoredIntegralVector(refShapeBasis);
	}

	//---------------------------------------------------------------------//
	// This is f***ing useless, it should be automatic due to inheritance! //
	// But without that it doesn't compile for some reason :-(             //
	//---------------------------------------------------------------------//

	double Integral(DomFunction globalFunction) const override
	{
		return PhysicalShape<Dim>::Integral(globalFunction);
	}
	double Integral(DomFunction globalFunction, int polynomialDegree) const override
	{
		return PhysicalShape<Dim>::Integral(globalFunction, polynomialDegree);
	}
};