#pragma once
class BasisFunction1D
{
public:
	virtual double Eval(double x) = 0;
	virtual double EvalGrad(double x) = 0;
};