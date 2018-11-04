#pragma once
class IBasisFunction1D
{
public:
	virtual double Eval(double x) = 0;
	virtual double EvalDerivative(double x) = 0;
};