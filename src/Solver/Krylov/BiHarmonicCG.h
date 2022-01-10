#pragma once
#include "../IterativeSolver.h"
#include "../../HHO/BiHarmonicMixedForm_HHO.h"
using namespace std;

template <int Dim>
class BiHarmonicCG : public IterativeSolver
{
private:
	BiHarmonicMixedForm_HHO<Dim>& _biHarPb;

public:
	BiHarmonicCG(BiHarmonicMixedForm_HHO<Dim>& biHarPb) :
		_biHarPb(biHarPb)
	{}

	virtual void Serialize(ostream& os) const override
	{
		os << "Conjugate Gradient for the bi-harmonic problem";
	}

	Vector Solve()
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
		// Conjugate Gradient //
		//--------------------//

		result.SetResidualNorm(r.norm());               result.AddWorkInFlops(Cost::Norm(r));

		// Direction of research
		Vector p = r;

		double r_dot_r = r.dot(r);

		this->IterationCount = 0;

		if (this->PrintIterationResults)
			cout << result << endl;

		while (!StoppingCriteriaReached(result))
		{
			result = IterationResult(result);

			// Solve 1st diffusion problem (f=0, Neum=p)
			// compatibility condition: <p|1> = 0
			//_integralOnBoundary.Enforce(p);
			Vector delta = _biHarPb.Solve1stDiffProblemWithZeroSource(p);

			// Solve 2nd diffusion problem (f=delta, Neum=0)
			// compatibility condition: (delta|1) = 0
			//_zeroMeanForReconstruct.Enforce(delta);
			Vector gamma = _biHarPb.Solve2ndDiffProblem(delta, true);


			// Step for theta in the direction of research
			double rho = p.dot(r) / gamma.dot(p);
			// Move theta in the direction of research
			theta += rho * p;

			double r_dot_r_old = r_dot_r; // save the dot product before overwriting r and p

			// Update residual
			r -= rho * gamma;

			r_dot_r = r.dot(r);

			// Step for the direction of research
			double q = r_dot_r / r_dot_r_old;
			// Update the direction of research
			p = r + q * p;

			//------------------------------------

			this->IterationCount++;

			result.SetX(theta);
			result.SetResidualNorm(r.norm());                   result.AddWorkInFlops(Cost::Norm(r));

			if (this->PrintIterationResults)
				cout << result << endl;
		}

		if (this->PrintIterationResults)
			cout << endl;

		this->SolvingComputationalWork = result.SolvingComputationalWork();

		return theta;
	}

};