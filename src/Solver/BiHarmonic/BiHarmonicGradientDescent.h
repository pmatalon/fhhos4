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

		result.SetResidualNorm(r.norm());

		this->IterationCount = 0;

		if (this->PrintIterationResults)
			cout << result << endl;

		while (!StoppingCriteriaReached(result))
		{
			result = IterationResult(result);

			// Solve 1st problem (f=source, Neum=theta --> lambda s.t. (lambda|1)=0)
			lambda = _biHarPb.Solve1stDiffProblem(theta);

			// Solve 2nd problem (f=lamda, Neum=0 --> r s.t. <r|1>=0) //
			r = _biHarPb.Solve2ndDiffProblem(lambda, true);

			// Compute the step
			//step = abs((theta - theta_old).dot(r - r_old))

			// Update theta in the opposite direction of the gradient (= r)
			theta -= _step * r;
			
			//------------------------------------

			this->IterationCount++;

			result.SetX(theta);
			result.SetResidualNorm(r.norm());

			if (this->PrintIterationResults)
				cout << result << endl;
		}

		if (this->PrintIterationResults)
			cout << endl;

		this->SolvingComputationalWork = result.SolvingComputationalWork();

		return theta;
	}

};