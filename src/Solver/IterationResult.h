#pragma once
#include "../Utils/Timer.h"
using namespace std;

class IterationResult
{
private:
	double _bNorm;
	bool _computeError = false;
	Vector _exactSolution;

	double _oldResidualNorm = -1;
	double _iterationConvRate = 0;
	list<double> _previousItConvRates;
	double _asymptoticConvRate = 0;
	double _tolerance = 1e-8;

	BigNumber _iterationComputationalWork = 0;
	BigNumber _solvingComputationalWork = 0;
	Timer _solvingTimer;
	Vector x;
	Vector r;
	Vector e;
public:
	int IterationNumber = 0;
	double ResidualNorm = -1;
	double NormalizedResidualNorm = -1;
	double RelativeErrorNorm = -1;

	IterationResult()
	{
		this->_solvingTimer.Start();
	}

	IterationResult(const IterationResult& oldResult)
		: _solvingTimer(oldResult._solvingTimer)
	{
		this->IterationNumber = oldResult.IterationNumber + 1;
		this->_solvingComputationalWork = oldResult._solvingComputationalWork;
		this->_bNorm = oldResult._bNorm;
		this->_computeError = oldResult._computeError;
		this->_exactSolution = oldResult._exactSolution;
		this->_oldResidualNorm = oldResult.NormalizedResidualNorm;
		this->_previousItConvRates = oldResult._previousItConvRates;
		this->_tolerance = oldResult._tolerance;
	}

	BigNumber SolvingComputationalWork()
	{
		return _solvingComputationalWork;
	}

	void SetB(const Vector& b)
	{
		this->_bNorm = b.norm();
	}

	Vector X() const
	{
		return this->x;
	}

	void SetX(const Vector& x)
	{
		this->x = x;
		if (_computeError)
			ComputeError();
		this->_solvingTimer.Stop();
	}

	void SetExactSolution(const Vector exactSolution)
	{
		this->_exactSolution = exactSolution;
		this->_computeError = true;
	}

	void SetTolerance(double tol)
	{
		this->_tolerance = tol;
	}

	void SetResidual(const Vector& r)
	{
		this->r = r;
		this->ResidualNorm = r.norm();
		this->NormalizedResidualNorm = _bNorm > 0 ? ResidualNorm / _bNorm : ResidualNorm;

		if (_oldResidualNorm != -1)
		{
			this->_iterationConvRate = this->NormalizedResidualNorm / _oldResidualNorm;

			if (this->_previousItConvRates.size() == 5)
				this->_previousItConvRates.pop_front();
			this->_previousItConvRates.push_back(this->_iterationConvRate);
			this->_asymptoticConvRate = 1;
			for (double r : _previousItConvRates)
				this->_asymptoticConvRate *= r;
			this->_asymptoticConvRate = pow(this->_asymptoticConvRate, 1.0 / _previousItConvRates.size());
		}
	}

	bool IsResidualSet()
	{
		return ResidualNorm != -1;
	}

	void ComputeError()
	{
		this->e = _exactSolution - x;
		double exactSolNorm = _exactSolution.norm();
		this->RelativeErrorNorm = exactSolNorm > 0 ? e.norm() / _exactSolution.norm() : e.norm();
	}

	void AddCost(BigNumber cost)
	{
		_iterationComputationalWork += cost;
		_solvingComputationalWork += cost;
	}

	friend ostream& operator<<(ostream& os, const IterationResult& result)
	{
		int IterWidth = 3;
		int normalizedResWidth = 17;
		int relativeErrorWidth = 17;
		int convRateWidth = 16;
		int predictedIterationsWidth = 14;
		int computWorkWidth = 15;
		int cpuTimeWidth = 10;
		int remainingTimeWidth = 16;
		if (result.IterationNumber == 0)
		{
			os << setw(IterWidth);
			os << "It";

			os << setw(normalizedResWidth);
			os << "Norm. residual";

			if (result.RelativeErrorNorm != -1)
			{
				os << setw(relativeErrorWidth);
				os << "Relative err.";
			}

			//os << setw(convRateWidth);
			//os << "It. conv. rate";

			os << setw(convRateWidth);
			os << "Asymp. cv rate";

			os << setw(predictedIterationsWidth);
			os << "Predic. it.";

			os << setw(computWorkWidth);
			os << "Comput. work";

			os << setw(cpuTimeWidth);
			os << "CPU time";

			os << setw(remainingTimeWidth);
			os << "Remaining time";
		}
		else
		{
			os << setw(IterWidth);
			os << result.IterationNumber;

			os << setw(normalizedResWidth);
			os << std::scientific << result.NormalizedResidualNorm;
			
			if (result.RelativeErrorNorm != -1)
			{
				os << setw(relativeErrorWidth);
				os << result.RelativeErrorNorm;
			}

			//os << setw(convRateWidth);
			//os << std::defaultfloat << result._iterationConvRate;

			os << setw(convRateWidth);
			os << std::defaultfloat << result._asymptoticConvRate;

			os << setw(predictedIterationsWidth);
			int remainingIterations = abs(ceil(log(result._tolerance / result.NormalizedResidualNorm) / log(result._asymptoticConvRate)));
			os << result.IterationNumber + remainingIterations;

			os << setw(computWorkWidth);
			os << result._solvingComputationalWork;

			os << setw(cpuTimeWidth);
			os << result._solvingTimer.CPU().InMilliseconds;

			os << setw(remainingTimeWidth);
			Duration d(result._solvingTimer.CPU().InMilliseconds / result.IterationNumber * remainingIterations);
			stringstream ss;
			ss << d;
			os << ss.str().substr(0, 8);
		}
		return os;
	}
};