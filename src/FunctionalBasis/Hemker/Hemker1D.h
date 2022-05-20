#pragma once
#include "../BasisFunction.h"
#include "../../Utils/Utils.h"
#include <math.h>
#include <assert.h>
using namespace std;

class Hemker1D : public IBasisFunction1D
{
private:
	int _maxDegree;
	int _degree;
	int _i;
	double _coeff;
public:
	Hemker1D(int maxDegree, int i)
	{
		assert(maxDegree >= 0 && i <= maxDegree);
		this->LocalNumber = i;
		this->_maxDegree = maxDegree;
		this->_degree = ComputeDegree(maxDegree, i);
		this->_i = i;
		this->_coeff = ComputeCoeff();
	}

	int ComputeDegree(int maxDegree, int i)
	{
		if (maxDegree <= 1)
			return maxDegree;
		// maxDegree > 1
		if (i <= 1)
			return 1;
		// i > 1
		if (maxDegree == 2 && i == 2)
			return 2;
		// maxDegree > 2 or i != 2
		if (i <= 3)
			return 3;
		// i > 3
		return i;
	}

	int GetDegree() const
	{
		return this->_degree;
	}

	double Eval(double x) const
	{
		if (_maxDegree == 0)
			return 1;
		if (_i == 0)
			return x - 1;
		if (_i == 1)
			return x + 1;
		if (_i == 2)
		{
			if(_maxDegree == 2)
				return (x-1) * (x+1);
			return pow(x-1, 2) * (x+1);
		}
		if(_i == 3)
			return (x-1) * pow(x+1, 2);
		return Hemker(x);
	}

	double EvalDerivative(double x) const
	{
		if (_maxDegree == 0)
			return 0;
		if (_i <= 1)
			return 1;
		if (_i == 2)
		{
			if(_maxDegree == 2)
				return 2*x;
			return (x-1) * (3*x+1);
		}
		if(_i == 3)
			return (x+1) * (3*x-1);
		return DHemker(x);
	}

	string ToString()
	{
		return "";
	}

	string ToString(string var)
	{
		if (_maxDegree == 0)
			return "1";
		if (_i == 0)
			return "(" + var + "-1)";
		if (_i == 1)
			return "(" + var + "+1)";
		if (_i == 2)
		{
			if(_maxDegree == 2)
				return "(" + var + "-1) * (" + var + "+1)";
			return "(" + var + "-1)^2 * (" + var + "+1)";
		}
		if(_i == 3)
			return "(" + var + "-1) * (" + var + "+1)^2";
		return "Hemker(" + to_string(this->_i) + ", " + var + ")";
	}

private:
	double ComputeCoeff()
	{
		if (_degree > 3)
		{
			assert (_i == _degree);
			double coeff = sqrt(2*_degree+1) * sqrt(Utils::Factorial(_degree-4)) * sqrt(Utils::Factorial(_degree+4) ) / ( pow(2, _degree) * sqrt(2) * Utils::Factorial(_degree) );
			return coeff;
		}
		return 1;
	}

	double Hemker(double x) const
	{
		const int n = _degree - 4;
		assert(n >= 0);
		double jac = 0;
		for (int l = 0; l <= n; l++)
		{
			jac += pow(x-1, n-l) * pow(x+1, l) * Utils::Binomial(_degree, l) * Utils::Binomial(_degree, n-l);
		}
		jac *= _coeff * pow(x-1, 2) * pow(x+1, 2);
		return jac;
	}

	double DHemker(double x) const
	{
		const int n = _degree - 4;
		assert(n >= 0);
		double djac = 0;
		for (int l = 0; l <= n; l++)
		{
			const double d1 = (n-l+2) * pow(x-1, n-l+1) * pow(x+1, l+2);
			const double d2 = (l+2) * pow(x-1, n-l+2) * pow(x+1, l+1);
			djac += (d1 + d2) * Utils::Binomial(_degree, l) * Utils::Binomial(_degree, n-l);
		}
		djac *= _coeff;
		return djac;
	}

};
