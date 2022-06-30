#pragma once
#include "../IterativeSolver.h"
#include "../../Discretizations/HHO/BiHarmonicMixedForm.h"
using namespace std;

class BiHarmonicCG : public IterativeSolver
{
private:
	BiHarmonicMixedForm* _biHarPb;
	ToleranceStrategy _toleranceStgy = ToleranceStrategy::DynamicVariableStep;
	double _diffSolverStartingTol = 1e-3;
	double _diffSolverToleranceStep = 1e-3;

	int _restartPeriod = 0;

public:
	Preconditioner* Precond = nullptr;

	BiHarmonicCG(BiHarmonicMixedForm* biHarPb, ToleranceStrategy toleranceStgy, double diffSolverToleranceStep = 1e-3, int restartPeriod = 0)
	{
		_biHarPb = biHarPb;
		_toleranceStgy = toleranceStgy;
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
		if (_toleranceStgy == ToleranceStrategy::Fixed)
			_biHarPb->SetDiffSolverTolerance(this->Tolerance);
		else
			_biHarPb->SetDiffSolverTolerance(max(_diffSolverStartingTol, this->Tolerance));
	}

	void Solve(const Vector& b, Vector& theta, bool xEquals0) override
	{
		this->SolvingComputationalWork = 0;

		//--------------------//
		// Conjugate Gradient //
		//--------------------//

		IterationResult result = CreateFirstIterationResult(b, theta);

		Vector r = xEquals0 ? b : b-A(theta);
		
		double r_dot_r = r.dot(r);

		result.SetResidualNorm(sqrt(r_dot_r));

		Vector z = Precond->Apply(r);

		// Direction of research
		Vector p = z;

		double r_dot_z = r.dot(z);

		this->IterationCount = 0;

		if (this->PrintIterationResults)
			cout << result << endl;

		bool restartNow = false;

		while (!StoppingCriteriaReached(result))
		{
			result = IterationResult(result);

			if (restartNow)
			{
				// Restart algorithm
				r = b - A(theta);
				z = Precond->Apply(r);
				p = z;
				r_dot_z = r.dot(z);

				restartNow = false;
			}

			Vector Ap = A(p);

			Timer cgWithoutATimer;
			cgWithoutATimer.Start();

			// Step for theta in the direction of research
			double alpha = r_dot_z / p.dot(Ap);

			// Move theta in the direction of research
			theta += alpha * p;

			double r_dot_z_old = r_dot_z; // save the dot product before overwriting r

			// Update residual
			r -= alpha * Ap;
			z = Precond->Apply(r);

			r_dot_z = r.dot(z);

			// Step for the direction of research
			double beta = r_dot_z / r_dot_z_old;
			// Update the direction of research
			p = z + beta * p;

			//-------------------------------------------

			this->IterationCount++;

			r_dot_r = r.dot(r);

			cgWithoutATimer.Stop();

			result.SetX(theta);
			result.SetResidualNorm(sqrt(r_dot_r));

			// Update dynamic tolerance for the diffusion solver
			bool toleranceUpdated = false;
			if (_toleranceStgy != ToleranceStrategy::Fixed && _biHarPb->IterativeDiffSolver() && result.NormalizedResidualNorm < _biHarPb->DiffSolverTolerance() && result.NormalizedResidualNorm > this->Tolerance)
			{
				UpdateTolerance(result, cgWithoutATimer);
				toleranceUpdated = true;
			}

			// Should restart?
			if (toleranceUpdated || (this->_restartPeriod > 0 && this->IterationCount % this->_restartPeriod == 0))
			{
				restartNow = true;
				result.AddAtTheEndOfTheLine("R");
				if (toleranceUpdated)
				{
					stringstream ss;
					ss << "(new tol = " << std::setprecision(1) << std::scientific << _biHarPb->DiffSolverTolerance() << ")";
					result.AddAtTheEndOfTheLine(ss.str());
				}
			}

			if (this->PrintIterationResults)
				cout << result << endl;
		}

		if (this->PrintIterationResults)
			cout << endl;

		this->SolvingComputationalWork = result.SolvingComputationalWork();
	}

private:
	Vector A(const Vector& x)
	{
		Vector delta = _biHarPb->Solve1stDiffProblemWithZeroSource(x);
		return -_biHarPb->Solve2ndDiffProblem(delta, true);
	}

	void UpdateTolerance(const IterationResult& cgResult, const Timer& cgWithoutDiffTimer)
	{
		double newTolerance;
		if (_toleranceStgy == ToleranceStrategy::DynamicFixedStep)
		{
			newTolerance = _biHarPb->DiffSolverTolerance() * _diffSolverToleranceStep;
		}
		else if (_toleranceStgy == ToleranceStrategy::DynamicVariableStep)
		{
			IterationResult& diffResult = _biHarPb->IterativeDiffSolver()->LastIterationResult;
			int cost1IterDiff = diffResult.IterationTimer().CPU().InMilliseconds;
			int cost1IterCGNoDiff = cgWithoutDiffTimer.CPU().InMilliseconds;
			double rhoCG = cgResult.AsymptoticConvRate;
			double coeff = 1;//Utils::ProgramArgs.Actions.DoubleParam1;
			double worsenRhoCG = log(1 + pow(rhoCG, coeff) * (exp(1) - 1));
			//cout << "rhoCG=" << std::setprecision(2) << rhoCG << ", worsenRhoCG=" << worsenRhoCG << endl;
			double lnRhoDiff = log(diffResult.AsymptoticConvRate);
			double lnRhoCG = log(worsenRhoCG);
			double lnTol = log(this->Tolerance);
			double lnRes = log(_biHarPb->DiffSolverTolerance());//log(cgResult.NormalizedResidualNorm);
			assert(cost1IterDiff > 0 && cost1IterCGNoDiff >= 0 && lnRhoDiff < 0 && lnRhoCG < 0);

			/*double a = 2 * cost1IterDiff / (lnRhoCG * lnRhoDiff);
			double b_low = cost1IterCGNoDiff / lnRhoCG + 2 * cost1IterDiff / lnRhoDiff * (1 - log(cgResult.NormalizedResidualNorm) / lnRhoCG);
			double b_high = b_low + 2 * cost1IterDiff * (1 / lnRhoCG + 1 / lnRhoDiff);

			double epsLow = exp(-b_low / (2 * a));
			double epsHigh = exp(-b_high / (2 * a));

			newTolerance = cgResult.NormalizedResidualNorm * epsLow;
			newTolerance = _biHarPb->DiffSolverTolerance() * epsLow;*/
			//double CDiff = 2 * cost1IterDiff / lnRhoDiff;

			/*double a = CDiff / lnRhoCG;
			double b = -CDiff * ((1 - lnRes/lnRhoCG) + lnTol/lnRhoCG);
			double c = CDiff * lnTol * (1 + lnTol / lnRhoCG);
			c -= ((lnTol-lnRes)/lnRhoCG + 1) * CDiff * lnTol;*/
			double a = 1;
			double b = lnRhoCG - lnRes - lnTol;
			double c = lnTol * lnRes;

			double det = b * b - 4 * a * c;
			if (det <= 0)
				newTolerance = this->Tolerance;
			else
			{
				double eps = exp(-b / (2 * a));
				newTolerance = eps;
			}
		}
		newTolerance = max(newTolerance, this->Tolerance);
		_biHarPb->SetDiffSolverTolerance(newTolerance);
	}

};