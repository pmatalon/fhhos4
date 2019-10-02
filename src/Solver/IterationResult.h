#pragma once
using namespace std;

class IterationResult
{
private:
	double _bNorm;
	bool _computeError = false;
	Vector _exactSolution;
	BigNumber _iterationComputationalWork = 0;
	BigNumber _solvingComputationalWork = 0;
	Vector x;
	Vector r;
	Vector e;
public:
	int IterationNumber = 0;
	double ResidualNorm = -1;
	double NormalizedResidualNorm = -1;
	double RelativeErrorNorm = -1;

	IterationResult()
	{}

	IterationResult(const IterationResult& oldResult)
	{
		this->IterationNumber = oldResult.IterationNumber + 1;
		this->_solvingComputationalWork = oldResult._solvingComputationalWork;
		this->_bNorm = oldResult._bNorm;
		this->_computeError = oldResult._computeError;
		this->_exactSolution = oldResult._exactSolution;
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
	}

	void SetExactSolution(const Vector exactSolution)
	{
		this->_exactSolution = exactSolution;
		this->_computeError = true;
	}

	void SetResidual(const Vector& r)
	{
		this->r = r;
		this->ResidualNorm = r.norm();
		this->NormalizedResidualNorm = _bNorm > 0 ? ResidualNorm / _bNorm : ResidualNorm;
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
		if (result.IterationNumber == 0)
		{
			os << "It.\tNormalized res";
			if (result.RelativeErrorNorm != -1)
				os << "\tRelative err";
			os << "\tComput. work";
		}
		else
		{
			os << result.IterationNumber << "\t" << std::scientific << result.NormalizedResidualNorm;
			if (result.RelativeErrorNorm != -1)
				os << "\t" << result.RelativeErrorNorm;
			os << "\t" << result._solvingComputationalWork;
		}
		return os;
	}
};