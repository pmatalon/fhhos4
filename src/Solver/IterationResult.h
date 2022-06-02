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
	MFlops _iterationComputationalWork = 0;
	MFlops _solvingComputationalWork = 0; // total of all iterations
	list<double> _previousItConvRates;
	double _tolerance = 0;
	Timer _solvingTimer;
	Vector e;
	string _textAtTheEnd;
public:
	function<void(IterationResult&, const Vector&)> OnNewSolution = nullptr;
	int IterationNumber = 0;
	double ResidualNorm = -1;
	double NormalizedResidualNorm = -1;
	double PreviousNormalizedResidualNorm = -1;
	double RelativeErrorNorm = -1;
	double L2Error = -1;
	double IterationConvRate = 0;
	double AsymptoticConvRate = 0;
	double BoundaryL2Norm = -1;

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
		this->PreviousNormalizedResidualNorm = oldResult.NormalizedResidualNorm;
		this->_previousItConvRates = oldResult._previousItConvRates;
		this->_tolerance = oldResult._tolerance;
		this->OnNewSolution = oldResult.OnNewSolution;
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
		if (OnNewSolution)
			OnNewSolution(*this, x);
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

		if (PreviousNormalizedResidualNorm != -1)
		{
			this->IterationConvRate = this->NormalizedResidualNorm / PreviousNormalizedResidualNorm;

			if (this->_previousItConvRates.size() == 5)
				this->_previousItConvRates.pop_front();
			this->_previousItConvRates.push_back(this->IterationConvRate);
			this->AsymptoticConvRate = 1;
			for (double cr : _previousItConvRates)
				this->AsymptoticConvRate *= cr;
			this->AsymptoticConvRate = pow(this->AsymptoticConvRate, 1.0 / _previousItConvRates.size());
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

	bool IsFineMatVecSet() const
	{
		return _oneFineMatVecWork > 0;
	}

	int NumberOfFineMatVec() const
	{
		if (IsFineMatVecSet())
			return (int)round(_solvingComputationalWork / _oneFineMatVecWork);
		return 0;
	}

	void AddAtTheEndOfTheLine(const string& text)
	{
		if (_textAtTheEnd.length() == 0)
			_textAtTheEnd = text;
		else
			_textAtTheEnd += " " + text;
	}

	friend ostream& operator<<(ostream& os, const IterationResult& result)
	{
		int IterWidth = 3;
		int predictedIterationsWidth = 6;
		int normalizedResWidth = 10;
		int relativeErrorWidth = 10;
		int l2ErrorWidth = 12;
		int boundaryL2NormWidth = 12;
		int iterationConvRateWidth = 9;
		int asymptoticConvRateWidth = 9;
		int computWorkWidth = 15;
		int nWorkUnitsWidth = 6;
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
			if (result.L2Error != -1)
			{
				os << setw(l2ErrorWidth);
				os << "L2-error";
			}
			if (result.BoundaryL2Norm != -1)
			{
				os << setw(boundaryL2NormWidth);
				os << "Boundary";
			}
			os << setw(iterationConvRateWidth);
			os << "Iter.";
			os << setw(asymptoticConvRateWidth);
			os << "Asymp.";
			//os << setw(computWorkWidth);
			//os << "Computational";
			//os << setw(cpuTimeWidth);
			//os << "";
			os << setw(nWorkUnitsWidth);
			os << "WUs";
			if (result._tolerance > 0)
			{
				os << setw(remainingTimeWidth);
				os << "Remaining";
			}
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
			if (result.L2Error != -1)
			{
				os << setw(l2ErrorWidth);
				os << "";
			}
			if (result.BoundaryL2Norm != -1)
			{
				os << setw(boundaryL2NormWidth);
				os << "L2-norm";
			}
			os << setw(iterationConvRateWidth);
			os << "cv rate";
			os << setw(asymptoticConvRateWidth);
			os << "cv rate";
			//os << setw(computWorkWidth);
			//os << "work";
			//os << setw(cpuTimeWidth);
			//os << "CPU time";
			os << setw(nWorkUnitsWidth);
			os << "";
			if (result._tolerance > 0)
			{
				os << setw(remainingTimeWidth);
				os << "time";
			}

			os << endl;
		}

		os << setw(IterWidth);
		os << result.IterationNumber;

		int remainingIterations = -1;
		os << setw(predictedIterationsWidth);
		if (result._tolerance > 0)
		{
			if (result.IterationNumber == 0)
				os << "";
			else
			{
				if (result.AsymptoticConvRate < 1)
				{
					remainingIterations = abs(ceil(log(result._tolerance / result.NormalizedResidualNorm) / log(result.AsymptoticConvRate)));
					os << "/ " + to_string(result.IterationNumber + remainingIterations);
				}
				else
					os << "/ -";
			}
		}
		else
			os << "";

		os << setw(normalizedResWidth);
		os << std::scientific << result.NormalizedResidualNorm;

		if (result.RelativeErrorNorm != -1)
		{
			os << setw(relativeErrorWidth);
			os << result.RelativeErrorNorm;
		}

		if (result.L2Error != -1)
		{
			os << setw(l2ErrorWidth);
			os << std::scientific << result.L2Error;
		}

		if (result.BoundaryL2Norm != -1)
		{
			os << setw(boundaryL2NormWidth);
			os << std::scientific << result.BoundaryL2Norm;
		}

		os << setw(iterationConvRateWidth);
		if (result.IterationNumber == 0)
			os << " ";
		else
			os << std::fixed << result.IterationConvRate;

		os << setw(asymptoticConvRateWidth);
		if (result.IterationNumber == 0)
			os << " ";
		else
			os << std::fixed << result.AsymptoticConvRate;

		//os << setw(computWorkWidth);
		//os << round(result._solvingComputationalWork);

		//os << setw(cpuTimeWidth);
		//os << result._solvingTimer.CPU().InMilliseconds;

		os << setw(nWorkUnitsWidth);
		os << result.NumberOfFineMatVec();

		if (result._tolerance > 0)
		{
			os << setw(remainingTimeWidth);
			if (result.IterationNumber == 0)
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

		if (result._textAtTheEnd.length() > 0)
			os << "  " << result._textAtTheEnd;

		return os;
	}
};