#pragma once
#include "Diffusion_HHO.h"


// --- Zero-mean condition (to enforce unicity of the solution)
// 
// 
// Let v be a function represented by the vector of coefficients x.
// Applying the zero mean condition is equivalent to orthogonalizing v to 1 (the constant function 1):
//              v <- v - (v|1)/(1|1)*1
// Let '_one' be the vector of coefficients corresponding to the constant function 1. 
// Then the discrete conterpart gives
//              x <- x - x^T*[(phi_i|1)]_i / (1|1) * _one       (eq.1)
// Let ipw1 := [(phi_i|1)]_i  (stands for "inner prods with one")
// In a setup function, we compute and store the vectors
//              _one   := M^-1 * ipw1,
//              _nipw1 := ipw1 / (1|1),
// where M^-1 is the face mass matrix and (1|1).
// (eq.1) becomes 
//              x <- x - x^T*_nipw1 * _one                      (eq.2)
// which is implemented in the Enforce() function.


class ZeroMeanEnforcer
{
protected:
	Vector _nipw1;
	Vector _one;
public:
	void Setup()
	{
		Vector ipw1 = BasisInnerProdWith1();
		_one = SolveMassMatrix(ipw1);
		_nipw1 = std::move(ipw1);
		_nipw1 /= OneScalOne();
	}
	virtual void Enforce(Vector& x)
	{
		assert(_nipw1.rows() > 0 && "Setup() must be called before using Enforce()");
		assert(x.rows() == _nipw1.rows());
		x -= x.dot(_nipw1) * _one;
	};
protected:
	// returns [(phi_i|1)]_i
	virtual Vector BasisInnerProdWith1() = 0;
	// returns M^-1 * x
	virtual Vector SolveMassMatrix(const Vector& x) = 0;
	// returns (1|1)
	virtual double OneScalOne() = 0;
};

//------------------------------------------//
// Reconstructed polynomial over the domain //
//------------------------------------------//

template<int Dim>
class ZeroMeanEnforcerFromReconstructCoeffs : public ZeroMeanEnforcer
{
private:
	Diffusion_HHO<Dim>* _discretePb;
public:
	ZeroMeanEnforcerFromReconstructCoeffs() {}
	ZeroMeanEnforcerFromReconstructCoeffs(Diffusion_HHO<Dim>* discretePb) : _discretePb(discretePb) {}
private:
	Vector BasisInnerProdWith1() override
	{
		return _discretePb->InnerProdWithReconstructBasis(Utils::ConstantFunctionOne);
	}
	Vector SolveMassMatrix(const Vector& x) override
	{
		return _discretePb->SolveReconstructMassMatrix(x);
	}
	double OneScalOne() override
	{
		return _discretePb->_mesh->Measure();
	}
};


//------------------------------//
// Polynomial over the skeleton //
//------------------------------//

template<int Dim>
class ZeroMeanEnforcerFromFaceCoeffs : public ZeroMeanEnforcer
{
private:
	Diffusion_HHO<Dim>* _discretePb;
public:
	ZeroMeanEnforcerFromFaceCoeffs() {}
	ZeroMeanEnforcerFromFaceCoeffs(Diffusion_HHO<Dim>* discretePb) : _discretePb(discretePb) {}
private:
	Vector BasisInnerProdWith1() override
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
	}
};

//------------------------------//
// Polynomial over the boundary //
//------------------------------//

template<int Dim>
class ZeroMeanEnforcerFromBoundaryFaceCoeffs : public ZeroMeanEnforcer
{
private:
	Diffusion_HHO<Dim>* _discretePb;
public:
	ZeroMeanEnforcerFromBoundaryFaceCoeffs() {}
	ZeroMeanEnforcerFromBoundaryFaceCoeffs(Diffusion_HHO<Dim>* discretePb) : _discretePb(discretePb) {}
private:
	Vector BasisInnerProdWith1() override
	{
		return _discretePb->InnerProdWithBoundaryFaceBasis(Utils::ConstantFunctionOne);
	}
	Vector SolveMassMatrix(const Vector& x) override
	{
		return _discretePb->SolveBoundaryFaceMassMatrix(x);
	}
	double OneScalOne() override
	{
		return _discretePb->_mesh->BoundaryMeasure();
	}
};