#pragma once
#include "../IterativeSolver.h"
#include "../../HHO/BiHarmonicMixedForm_HHO.h"
using namespace std;

template <int Dim>
class BiHarmonicGradientDescent : public IterativeSolver
{
private:
	BiHarmonicMixedForm_HHO<Dim>& _biHarPb;
	double _step;

public:
	BiHarmonicGradientDescent(BiHarmonicMixedForm_HHO<Dim>& biHarPb, double step) :
		_biHarPb(biHarPb)
	{
		_step = step;
	}

	virtual void Serialize(ostream& os) const override
	{
		os << "Conjugate Gradient for the bi-harmonic problem";
	}

	Vector Solve() override
	{
		this->SolvingComputationalWork = 0;

		// Find initial theta verifying the compatibility condition
		//               (source|1) + <theta|1> = 0
		Vector theta = _biHarPb.FindCompatibleTheta();

		IterationResult result = CreateFirstIterationResult(Vector::Zero(theta.rows()), theta);

		// Solve 1st problem (f=source, Neum=theta --> lambda s.t. (lambda|1)=0)
		Vector lambda = _biHarPb.Solve1stDiffProblem(theta);

		// Solve 2nd problem (f=lamda, Neum=0 --> r s.t. <r|1>=0) //
		Vector r = _biHarPb.Solve2ndDiffProblem(lambda, true);

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
			theta -= step * r;

			// Solve 1st problem (f=source, Neum=theta --> lambda s.t. (lambda|1)=0)
			lambda = _biHarPb.Solve1stDiffProblem(theta);

			// Solve 2nd problem (f=lamda, Neum=0 --> r s.t. <r|1>=0) //
			r = _biHarPb.Solve2ndDiffProblem(lambda, true);

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

		return theta;
	}

private:
	double L2InnerProdOnBoundary(const Vector& v1, const Vector& v2)
	{
		return _biHarPb.DiffPb().L2InnerProdOnBoundary(v1, v2);
		//return v1.dot(v2);
	}
};