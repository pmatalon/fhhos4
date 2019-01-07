#pragma once
#include "IBasisFunction.h"
#include <math.h>
using namespace std;

class TensorPolynomial2D : public IBasisFunction2D
{
private:
	IBasisFunction1D* _funcX;
	IBasisFunction1D* _funcY;
public:

	TensorPolynomial2D(IBasisFunction1D* funcX, IBasisFunction1D* funcY)
	{
		this->_funcX = funcX;
		this->_funcY = funcY;
	}

	RefInterval ReferenceInterval() { return this->_funcX->ReferenceInterval(); }

	int GetDegree()
	{
		return this->_funcX->GetDegree() + this->_funcY->GetDegree();
	}

	double Eval(double x, double y)
	{
		return this->_funcX->Eval(x) * this->_funcY->Eval(y);
	}

	double EvalGradX(double x, double y)
	{
		return this->_funcX->EvalDerivative(x) * this->_funcY->Eval(y);
	}

	double EvalGradY(double x, double y)
	{
		return this->_funcX->Eval(x) * this->_funcY->EvalDerivative(y);
	}

	string ToString()
	{
		string polyX = this->_funcX->ToString("X");
		string polyY = this->_funcY->ToString("Y");
		if (polyX.compare("1") == 0)
			return polyY;
		if (polyY.compare("1") == 0)
			return polyX;
		return polyX + " * " + polyY;
	}
};

