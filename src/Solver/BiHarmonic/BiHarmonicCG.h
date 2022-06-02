#pragma once
#include "../IterativeSolver.h"
#include "../../Discretizations/HHO/BiHarmonicMixedForm.h"
using namespace std;

class BiHarmonicCG : public IterativeSolver
{
private:
	BiHarmonicMixedForm* _biHarPb;
	double _diffSolverToleranceStep = 1e-3;
	int _restartPeriod = 0;

public:
	Preconditioner* Precond = nullptr;

	BiHarmonicCG(BiHarmonicMixedForm* biHarPb, double diffSolverToleranceStep = 1e-3, int restartPeriod = 0)
	{
		_biHarPb = biHarPb;
		_diffSolverToleranceStep = diffSolverToleranceStep;
		_restartPeriod = restartPeriod;
	}

	void Serialize(ostream& os) const override
	{
		os << "Conjugate Gradient" << endl;
		os << "        Laplacian solver tolerance: ";
		if (_diffSolverToleranceStep > 0)
			os << "dynamic (step = " << std::scientific << std::setprecision(1) << _diffSolverToleranceStep << ")";
		else
			os << "constant (" << std::scientific << std::setprecision(1) << this->Tolerance << ")";
		if (_restartPeriod > 0)
			os << endl << "        Restart period = " << _restartPeriod;

		if (Precond)
			os << endl << "Preconditioner: " << *Precond;
	}

	void Setup(const SparseMatrix& A) override
	{
		IterativeSolver::Setup(A);
		if (_diffSolverToleranceStep > 0)
			_biHarPb->SetDiffSolverTolerance(max(_diffSolverToleranceStep, this->Tolerance));
		else
			_biHarPb->SetDiffSolverTolerance(this->Tolerance);
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

			// Update residual
			r -= alpha * Ap;
			z = Precond->Apply(r);

			r_dot_z = L2InnerProdOnBoundary(r, z);

			// Step for the direction of research
			double beta = r_dot_z / r_dot_z_old;
			// Update the direction of research
			p = z + beta * p;

			//-------------------------------------------

			this->IterationCount++;

			r_dot_r = L2InnerProdOnBoundary(r, r);

			result.SetX(theta);
			result.SetResidualNorm(sqrt(r_dot_r));

			// Update dynamic tolerance for the diffusion solver
			bool toleranceUpdated = false;
			if (_diffSolverToleranceStep > 0 && result.NormalizedResidualNorm < _biHarPb->DiffSolverTolerance() && result.NormalizedResidualNorm > this->Tolerance)
			{
				double newTolerance = max(_biHarPb->DiffSolverTolerance() * _diffSolverToleranceStep, this->Tolerance);
				_biHarPb->SetDiffSolverTolerance(newTolerance);
				toleranceUpdated = true;
			}

			// Restart?
			if (toleranceUpdated || (this->_restartPeriod > 0 && this->IterationCount % this->_restartPeriod == 0))
			{
				result.AddAtTheEndOfTheLine("R");
				if (toleranceUpdated)
				{
					stringstream ss;
					ss << "(new tol = " << std::setprecision(1) << std::scientific << _biHarPb->DiffSolverTolerance() << ")";
					result.AddAtTheEndOfTheLine(ss.str());
				}

				// Restart algorithm
				r = b - A(theta);
				z = Precond->Apply(r);
				r_dot_z = L2InnerProdOnBoundary(r, z);

				p = z;
			}

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