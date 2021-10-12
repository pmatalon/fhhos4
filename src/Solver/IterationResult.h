#pragma once
#include "../Utils/Timer.h"
#include "../Utils/Cost.h"
using namespace std;

class IterationResult
{
private:
	double _bNorm = -1;
	bool _computeError = false;
	Vector _exactSolution;
	MFlops _oneFineMatVecWork = 0;

	double _oldResidualNorm = -1;
	double _iterationConvRate = 0;
	list<double> _previousItConvRates;
	double _asymptoticConvRate = 0;
	double _tolerance = 1e-8;
	MFlops _iterationComputationalWork = 0;
	MFlops _solvingComputationalWork = 0; // total of all iterations
	Timer _solvingTimer;
	Vector e;
public:
	int IterationNumber = 0;
	double ResidualNorm = -1;
	double NormalizedResidualNorm = -1;
	double RelativeErrorNorm = -1;

	Vector Residual;
	Vector Ax;

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
		this->_oneFineMatVecWork = oldResult._oneFineMatVecWork;
		this->_oldResidualNorm = oldResult.NormalizedResidualNorm;
		this->_previousItConvRates = oldResult._previousItConvRates;
		this->_tolerance = oldResult._tolerance;
	}

	// Move assignment operator
	IterationResult& operator=(IterationResult&& result) = default;

	MFlops SolvingComputationalWork()
	{
		return _solvingComputationalWork;
	}

	void SetA(const SparseMatrix& A)
	{
		this->_oneFineMatVecWork = Cost::MatVec(A) * 1e-6;
	}

	void SetB(const Vector& b)
	{
		this->_bNorm = b.norm();
	}

	void SetX(const Vector& x)
	{
		if (_computeError)
			ComputeError(x);
		this->_solvingTimer.Stop();
	}

	void SetExactSolution(const Vector& exactSolution)
	{
		this->_exactSolution = exactSolution;
		this->_computeError = true;
	}

	void SetTolerance(double tol)
	{
		this->_tolerance = tol;
	}

	void SetResidualNorm(double rNorm)
	{
		this->ResidualNorm = rNorm;
		assert(_bNorm >= 0);
		this->NormalizedResidualNorm = _bNorm > 0 ? ResidualNorm / _bNorm : ResidualNorm;

		if (IterationNumber > 1 && _oldResidualNorm != -1)
		{
			this->_iterationConvRate = this->NormalizedResidualNorm / _oldResidualNorm;

			if (this->_previousItConvRates.size() == 5)
				this->_previousItConvRates.pop_front();
			this->_previousItConvRates.push_back(this->_iterationConvRate);
			this->_asymptoticConvRate = 1;
			for (double cr : _previousItConvRates)
				this->_asymptoticConvRate *= cr;
			this->_asymptoticConvRate = pow(this->_asymptoticConvRate, 1.0 / _previousItConvRates.size());
		}
	}

	void SetResidualAsB()
	{
		this->ResidualNorm = _bNorm;
		this->NormalizedResidualNorm = 1;
	}

	bool IsResidualSet()
	{
		return ResidualNorm != -1;
	}

	void ComputeError(const Vector& x)
	{
		this->e = _exactSolution - x;
		double exactSolNorm = _exactSolution.norm();
		this->RelativeErrorNorm = exactSolNorm > 0 ? e.norm() / _exactSolution.norm() : e.norm();
	}

	void AddWorkInFlops(Flops flops)
	{
		AddWorkInMFlops(flops * 1e-6); // conversion to Mflops
	}

	void AddWorkInMFlops(MFlops megaFlops)
	{
		_iterationComputationalWork += megaFlops;
		_solvingComputationalWork += megaFlops;
	}

	MFlops IterationComputationalWork()
	{
		return _iterationComputationalWork;
	}

	friend ostream& operator<<(ostream& os, const IterationResult& result)
	{
		int IterWidth = 3;
		int predictedIterationsWidth = 6;
		int normalizedResWidth = 12;
		int relativeErrorWidth = 12;
		int iterationConvRateWidth = 11;
		int asymptoticConvRateWidth = 12;
		int computWorkWidth = 15;
		int nFineMatVecWidth = 8;
		int cpuTimeWidth = 10;
		int remainingTimeWidth = 11;

		cout.precision(2);
		if (result.IterationNumber == 0)
		{
			// Label row 1
			os << setw(IterWidth + predictedIterationsWidth);
			os << "";
			os << setw(normalizedResWidth);
			os << "Relative";
			if (result.RelativeErrorNorm != -1)
			{
				os << setw(relativeErrorWidth);
				os << "Relative";
			}
			os << setw(iterationConvRateWidth);
			os << "Iteration";
			os << setw(asymptoticConvRateWidth);
			os << "Asymptotic";
			//os << setw(computWorkWidth);
			//os << "Computational";
			//os << setw(cpuTimeWidth);
			//os << "";
			os << setw(nFineMatVecWidth);
			os << "Fine";
			os << setw(remainingTimeWidth);
			os << "Remaining";
			os << endl;

			// Label row 2
			os << setw(IterWidth + predictedIterationsWidth);
			os << "Iteration";
			os << setw(normalizedResWidth);
			os << "residual";
			if (result.RelativeErrorNorm != -1)
			{
				os << setw(relativeErrorWidth);
				os << "error";
			}
			os << setw(iterationConvRateWidth);
			os << "cv rate";
			os << setw(asymptoticConvRateWidth);
			os << "cv rate";
			//os << setw(computWorkWidth);
			//os << "work";
			//os << setw(cpuTimeWidth);
			//os << "CPU time";
			os << setw(nFineMatVecWidth);
			os << "MatVec";
			os << setw(remainingTimeWidth);
			os << "time";
		}
		else
		{
			os << setw(IterWidth);
			os << result.IterationNumber;

			os << setw(predictedIterationsWidth);
			int remainingIterations = -1;
			if (result.IterationNumber == 1)
				os << "";
			else
			{
				if (result._asymptoticConvRate < 1)
				{
					remainingIterations = abs(ceil(log(result._tolerance / result.NormalizedResidualNorm) / log(result._asymptoticConvRate)));
					os << "/ " + to_string(result.IterationNumber + remainingIterations);
				}
				else
					os << "/ -";
			}

			os << setw(normalizedResWidth);
			os << std::scientific << result.NormalizedResidualNorm;
			
			if (result.RelativeErrorNorm != -1)
			{
				os << setw(relativeErrorWidth);
				os << result.RelativeErrorNorm;
			}

			os << setw(iterationConvRateWidth);
			if (result.IterationNumber == 1)
				os << " ";
			else
				os << std::defaultfloat << result._iterationConvRate;

			os << setw(asymptoticConvRateWidth);
			if (result.IterationNumber == 1)
				os << " ";
			else
				os << std::defaultfloat << result._asymptoticConvRate;

			//os << setw(computWorkWidth);
			//os << round(result._solvingComputationalWork);

			//os << setw(cpuTimeWidth);
			//os << result._solvingTimer.CPU().InMilliseconds;

			os << setw(nFineMatVecWidth);
			os << (int)round(result._solvingComputationalWork / result._oneFineMatVecWork);

			os << setw(remainingTimeWidth);
			if (result.IterationNumber == 1)
				os << " ";
			else if (remainingIterations == -1)
				os << "-";
			else
			{
				Duration d(result._solvingTimer.CPU().InMilliseconds / result.IterationNumber * remainingIterations);
				stringstream ss;
				ss << d;
				os << ss.str().substr(0, 8);
			}
		}
		return os;
	}
};