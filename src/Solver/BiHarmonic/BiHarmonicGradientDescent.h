#pragma once
#include "../IterativeSolver.h"
#include "../../Discretizations/HHO/BiHarmonicMixedForm.h"
using namespace std;

class BiHarmonicGradientDescent : public IterativeSolver
{
private:
	BiHarmonicMixedForm* _biHarPb;

public:
	BiHarmonicGradientDescent(BiHarmonicMixedForm* biHarPb)
	{
		_biHarPb = biHarPb;
	}

	virtual void Serialize(ostream& os) const override
	{
		os << "Gradient descent";
	}

	void Solve(const Vector& b, Vector& theta, bool xEquals0) override
	{
		this->SolvingComputationalWork = 0;

		IterationResult result = CreateFirstIterationResult(b, theta);

		Vector r = xEquals0 ? b : b - A(theta);

		//--------------------//
		//  Gradient descent  //
		//--------------------//

		Vector r_old = r;

		result.SetResidualNorm(sqrt(L2InnerProdOnBoundary(r, r)));

		this->IterationCount = 0;

		if (this->PrintIterationResults)
			cout << result << endl;

		double step = 1;


		while (!StoppingCriteriaReached(result))
		{
			result = IterationResult(result);
			
			// Update theta in the opposite direction of the gradient (= r)
			theta += step * r;

			r = b - A(theta);

			// Compute the step
			Vector r_minus_r_old = r - r_old;
			step = abs(step * L2InnerProdOnBoundary(r_old, r_minus_r_old)) / L2InnerProdOnBoundary(r_minus_r_old, r_minus_r_old);

			r_old = r;
			
			//------------------------------------

			this->IterationCount++;

			result.SetX(theta);
			result.SetResidualNorm(sqrt(L2InnerProdOnBoundary(r, r)));

			if (this->PrintIterationResults)
				cout << result << endl;
		}

		if (this->PrintIterationResults)
			cout << endl;

		this->SolvingComputationalWork = result.SolvingComputationalWork();
	}

private:
	double L2InnerProdOnBoundary(const Vector& v1, const Vector& v2)
	{
		//return _biHarPb->L2InnerProdOnBoundary(v1, v2);
		return v1.dot(v2);
	}

	Vector A(const Vector& x)
	{
		Vector delta = _biHarPb->Solve1stDiffProblem_Homogeneous(x);
		return -_biHarPb->Solve2ndDiffProblem_Homogeneous(delta);
	}
};