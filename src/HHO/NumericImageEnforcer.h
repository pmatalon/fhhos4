#pragma once
#include "Diffusion_HHO.h"


// --- Numeric Image Enforcer in the case of a pure Neumann problem (to ensure b \in Im(A))
// 
// The fundamental theorem of linear algebra states that Im(A) is the orthogonal of Ker(A^T).
// Since A is symmetric, Im(A) is the orthogonal of Ker(A).
// 
// The constant function 1 is in the kernel of the bilinear form, so its vector of coefficients, denoted 'one', is in Ker(A).
// Let b be a vector. Projecting b onto Im(A) then corresponds to orthogonalizing b to _one:
//              b <- b - (b^T*one / one^T*one) * one.        (eq.1)
// 
// Implementation:
// 
// Let (.|.) denote the L2-inner product.
// Let ipw1 := [(phi_i|1)]_i  (stands for "inner prods with one").
// In a setup function, we compute
//                one := M^-1 * [(phi_i|1)]_i,
//              _nOne := one / sqrt(one^T*one),
// where M^-1 is the face mass matrix and (1|1).
// (eq.1) becomes 
//              b <- b - b^T*_nOne * _nOne.                      (eq.2)


class NumericImageEnforcer
{
protected:
	Vector _nOne;
public:
	void Setup()
	{
		Vector one = CoefficientsOfFunctionOne();
		double normOne = sqrt(one.dot(one));
		_nOne = std::move(one);
		_nOne /= normOne;
	}
	double OrthogonalityFactor(const Vector& b)
	{
		return b.dot(_nOne);
	}
	bool Check(const Vector& b)
	{
		return abs(OrthogonalityFactor(b)) < Utils::NumericalZero;
	}
	double CheckKernel(const SparseMatrix& A)
	{
		return (A * _nOne).norm();
	}
	void ProjectOntoImage(Vector& b)
	{
		b -= OrthogonalityFactor(b) * _nOne;
	}
protected:
	virtual Vector CoefficientsOfFunctionOne() = 0;
};



//------------------------------//
// Polynomial over the skeleton //
//------------------------------//

template<int Dim>
class NumericImageEnforcerFromFaceCoeffs : public NumericImageEnforcer
{
private:
	Diffusion_HHO<Dim>* _discretePb;
public:
	NumericImageEnforcerFromFaceCoeffs() {}
	NumericImageEnforcerFromFaceCoeffs(Diffusion_HHO<Dim>* discretePb) : _discretePb(discretePb) {}
private:
	Vector CoefficientsOfFunctionOne() override
	{
		return _discretePb->ProjectOnFaceDiscreteSpace(Utils::ConstantFunctionOne);
	}
	/*Vector BasisInnerProdWith1() override
	{
		return _discretePb->InnerProdWithFaceBasis(Utils::ConstantFunctionOne);
	}
	Vector SolveMassMatrix(const Vector& x) override
	{
		return _discretePb->SolveFaceMassMatrix(x);
	}
	double OneScalOne() override
	{
		return _discretePb->_mesh->SkeletonMeasure();
	}*/
};
