#pragma once
#include "../../../Utils/Utils.h"


class IDiscreteSpace
{
public:
	// Returns the number of DoFs in the space
	virtual BigNumber Dimension() = 0;

	// Returns the measure of the support, which equals to (1|1)
	virtual double Measure() = 0;

	// Returns [(phi_i|func)]_i
	virtual Vector InnerProdWithBasis(DomFunction func) = 0;

	// Returns M * x
	virtual Vector ApplyMassMatrix(const Vector& x) = 0;

	// Returns M^-1 * x
	virtual Vector SolveMassMatrix(const Vector& x) = 0;

	// Returns the vector of coefficients representing 'func' in the basis.
	// Equivalent to: SolveMassMatrix(InnerProdWithBasis(func))
	virtual Vector Project(DomFunction func) = 0;

	// Returns v1^T * M * v2
	virtual double L2InnerProd(const Vector& v1, const Vector& v2) = 0;

	virtual double Integral(const Vector& coeffs) = 0;

	virtual double Integral(DomFunction func) = 0;

	double MeanValue(const Vector& coeffs)
	{
		return Integral(coeffs) / Measure();
	}

	double MeanValue(DomFunction func)
	{
		return Integral(func) / Measure();
	}

	double L2Norm(const Vector& v)
	{
		return sqrt(L2InnerProd(v, v));
	}

	double RelativeL2Norm(const Vector& error, const Vector& discreteExactSolution)
	{
		return L2Norm(error) / L2Norm(discreteExactSolution);
	}

	double RelativeL2Error(const Vector& approximation, const Vector& discreteExactSolution)
	{
		return L2Norm(discreteExactSolution - approximation) / L2Norm(discreteExactSolution);
	}

	bool Contains(const Vector& v)
	{
		return v.rows() == Dimension();
	}

	Vector ZeroVector()
	{
		return Vector::Zero(Dimension());
	}
};
