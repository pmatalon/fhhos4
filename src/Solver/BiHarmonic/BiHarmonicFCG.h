#pragma once
#include "../IterativeSolver.h"
#include "../Preconditioner.h"
#include "../../Discretizations/HHO/BiHarmonicMixedForm.h"
using namespace std;

class BiHarmonicFCG : public IterativeSolver
{
private:
	BiHarmonicMixedForm* _biHarPb;
	ToleranceStrategy _toleranceStgy = ToleranceStrategy::DynamicVariableStep;
	double _diffSolverStartingTol = 1e-3;
	double _diffSolverToleranceStep = 1e-3;

	struct Direction
	{
		Vector* d;
		Vector* Ad;
		double d_dot_Ad;
	};

	int _truncation = 1; // max number of previous search direction to use for A-orthogonalization
	int _restartPeriod = 0;

public:
	Preconditioner* Precond = nullptr;

	BiHarmonicFCG(BiHarmonicMixedForm* biHarPb, ToleranceStrategy toleranceStgy, double diffSolverToleranceStep = 1e-3, int truncation = 1, int restartPeriod = 0)
	{
		_biHarPb = biHarPb;
		_toleranceStgy = toleranceStgy;
		_diffSolverToleranceStep = diffSolverToleranceStep;
		_restartPeriod = restartPeriod;
	}

	void Serialize(ostream& os) const override
	{
		os << "Flexible Conjugate Gradient" << endl;
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

	void Solve(const Vector& b, Vector& x, bool xEquals0) override
	{
		this->SolvingComputationalWork = 0;

		//-----------------------------//
		// Flexible Conjugate Gradient //
		//-----------------------------//

		IterationResult result = CreateFirstIterationResult(b, x);

		// Init restart parameters
		list<Direction> previousDirections;

		Vector& r = this->Residual;
		if (xEquals0)
		{
			r = b;
			result.SetResidualAsB();
		}
		else
		{
			r = b - A(x);
			result.SetResidualNorm(r.norm());
		}

		this->IterationCount = 0;

		if (this->PrintIterationResults)
			cout << result << endl;

		while (!StoppingCriteriaReached(result))
		{
			result = IterationResult(result);

			// Apply the preconditioner
			Vector z, Az;
			z = Precond->Apply(r);

			// The final direction of research, d, is found by 
			// A-orthogonalizing z with the previous directions of research, 
			// according to the truncation / restart method chosen.

			Vector* d = new Vector(z);

			for (Direction const& directionk : previousDirections)
			{
				Vector& dk = *directionk.d;
				Vector& Adk = *directionk.Ad;
				double dk_dot_Adk = directionk.d_dot_Ad;

				// A-orthogonalization
				double z_dot_Adk = z.dot(Adk);
				*d -= (z_dot_Adk / dk_dot_Adk) * dk;
			}
			assert(d->norm() > 0);

			Vector* Ad = new Vector(A(*d));

			double d_dot_Ad = d->dot(*Ad);
			assert(d_dot_Ad != 0);

			// Step length in the direction of research
			double alpha = r.dot(*d) / d_dot_Ad;

			// Moving from the current solution to the next
			// by taking the step in the direction of research
			x += alpha * (*d);

			// New residual
			r -= alpha * (*Ad);

			// Restart?
			if (previousDirections.size() == this->_truncation)
				DeleteOldest(previousDirections);

			if (this->_restartPeriod > 0 && this->IterationCount % this->_restartPeriod == 0)
			{
				Clear(previousDirections);
				r = b - A(x);
				result.AddAtTheEndOfTheLine("R");
			}

			// The direction of research is saved 
			// (even after a restart, as advised by Notay)
			Direction direction = { d, Ad, d_dot_Ad };
			previousDirections.push_back(direction);

			//---------- End of iteration ----------//
			this->IterationCount++;

			result.SetX(x);
			result.SetResidualNorm(r.norm());                         result.AddWorkInFlops(Cost::Norm(r));

			if (this->PrintIterationResults)
				cout << result << endl;
		}

		Clear(previousDirections);

		if (this->PrintIterationResults)
			cout << endl;

		result.CopyInfoInto(this->LastIterationResult);
		this->SolvingComputationalWork = result.SolvingComputationalWork();
	}

private:
	Vector A(const Vector& x)
	{
		return _biHarPb->ProblemOperator(x);
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
			int cost1IterDiff = diffResult.IterationTimer().CPU().InMilliseconds();
			int cost1IterCGNoDiff = cgWithoutDiffTimer.CPU().InMilliseconds();
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

	void Clear(list<Direction>& previousDirections)
	{
		for (Direction& directionk : previousDirections)
		{
			delete directionk.d;
			delete directionk.Ad;
		}
		previousDirections.clear();
	}

	void DeleteOldest(list<Direction>& previousDirections)
	{
		Direction& oldestDirection = previousDirections.front();
		delete oldestDirection.d;
		delete oldestDirection.Ad;
		previousDirections.pop_front();
	}
};