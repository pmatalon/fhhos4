#pragma once
#include "../IterativeSolver.h"
#include "../../Discretizations/HHO/BiHarmonicMixedForm.h"
using namespace std;

class BiHarmonicCG : public IterativeSolver
{
private:
	BiHarmonicMixedForm* _biHarPb;

public:
	Preconditioner* Precond;
	int Restart = 0;

	BiHarmonicCG(BiHarmonicMixedForm* biHarPb, int restart = 0)
	{
		_biHarPb = biHarPb;
		Restart = restart;
	}

	virtual void Serialize(ostream& os) const override
	{
		os << "Conjugate Gradient for the biharmonic problem";
		if (Restart > 0)
			os << " (restart = " << Restart << ")";
	}

	void Solve(const Vector& b, Vector& theta, bool xEquals0) override
	{
		this->SolvingComputationalWork = 0;

		//--------------------//
		// Conjugate Gradient //
		//--------------------//

		IterationResult result = CreateFirstIterationResult(b, theta);

		Vector r = xEquals0 ? b : b-A(theta);
		
		double r_dot_r = L2InnerProdOnBoundary(r, r);

		result.SetResidualNorm(sqrt(r_dot_r));

		Vector z = Precond->Apply(r);

		// Direction of research
		Vector p = z;

		double r_dot_z = r.dot(z);

		this->IterationCount = 0;

		if (this->PrintIterationResults)
			cout << result << endl;

		while (!StoppingCriteriaReached(result))
		{
			result = IterationResult(result);

			Vector Ap = A(p);

			// Step for theta in the direction of research
			double alpha = r_dot_z / L2InnerProdOnBoundary(p, Ap);

			// Move theta in the direction of research
			theta += alpha * p;

			double r_dot_z_old = r_dot_z; // save the dot product before overwriting r

			// Restart?
			if (this->IterationCount > 0 && this->Restart > 0 && this->IterationCount % this->Restart == 0)
			{
				// Restart algorithm
				r = b - A(theta);
				z = Precond->Apply(r);
				r_dot_z = L2InnerProdOnBoundary(r, z);

				p = z;
			}
			else
			{
				// Update residual
				r -= alpha * Ap;
				z = Precond->Apply(r);

				r_dot_z = L2InnerProdOnBoundary(r, z);

				// Step for the direction of research
				double beta = r_dot_z / r_dot_z_old;
				// Update the direction of research
				p = z + beta * p;
			}


			//------------------------------------
			/*
			// Recompute explicitly the residual, by computing the solution u
			Vector lambda_final = _biHarPb->Solve1stDiffProblemWithFSource(theta);
			Vector u_boundary_real = _biHarPb->Solve2ndDiffProblem(lambda_final, true);

			Vector r_u_boundary = r - u_boundary_real;
			//cout << "r = " << endl << r << endl;
			//cout << "u_boundary = " << endl << u_boundary << endl;
			//cout << "r - u_boundary = " << endl << r_u_boundary << endl;
			cout << "||u_boundary|| = " << sqrt(L2InnerProdOnBoundary(u_boundary_real, u_boundary_real)) << "   ";
			cout << "||r - u_boundary|| = " << sqrt(L2InnerProdOnBoundary(r_u_boundary, r_u_boundary)) << "   ";
			cout << "||r|| = " << sqrt(r_dot_r_old) << endl;
			*/
			//------------------------------------

			this->IterationCount++;

			r_dot_r = L2InnerProdOnBoundary(r, r);

			result.SetX(theta);
			result.SetResidualNorm(sqrt(r_dot_r));

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

	/*Vector A0(const Vector& x)
	{
		Vector delta = _biHarPb->Solve1stDiffProblem(x);
		return -_biHarPb->Solve2ndDiffProblem(delta, true);
		//return _biHarPb->DiffPb().A * x;
	}*/

	Vector A(const Vector& x)
	{
		Vector delta = _biHarPb->Solve1stDiffProblemWithZeroSource(x);
		return -_biHarPb->Solve2ndDiffProblem(delta, true);
		//return _biHarPb->DiffPb().A * x;
	}

};