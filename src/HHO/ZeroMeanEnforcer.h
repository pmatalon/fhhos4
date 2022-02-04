#pragma once
#include "Diffusion_HHO.h"
#include "HigherOrderBoundary.h"
#include "DiscreteSpaces/IDiscreteSpace.h"


// --- Zero-mean condition (to enforce unicity of the solution)
// 
// Let (.|.) denote the L2-inner product.
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
private:
	IDiscreteSpace* _discreteSpace = nullptr;
	Vector _nipw1;
	Vector _one;

public:
	ZeroMeanEnforcer() {}

	ZeroMeanEnforcer(IDiscreteSpace* discreteSpace)
	{
		_discreteSpace = discreteSpace;
	}

	void Setup()
	{
		Vector ipw1 = _discreteSpace->InnerProdWithBasis(Utils::ConstantFunctionOne);
		_one = _discreteSpace->SolveMassMatrix(ipw1);
		_nipw1 = std::move(ipw1);
		_nipw1 /= _discreteSpace->Measure(); // (1|1) = measure
	}

	double OrthogonalityFactor(const Vector& x)
	{
		return x.dot(_nipw1);
	}

	void Enforce(Vector& x)
	{
		assert(_nipw1.rows() > 0 && "Setup() must be called before using Enforce()");
		assert(x.rows() == _nipw1.rows());
		x -= OrthogonalityFactor(x) * _one;
	}

	bool Check(const Vector& x)
	{
		return abs(OrthogonalityFactor(x)) < Utils::Eps;
	}
};