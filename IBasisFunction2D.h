#pragma once
class IBasisFunction2D
{
public:
	virtual double Eval(double x, double y) = 0;
	virtual double EvalGradX(double x, double y) = 0;
	virtual double EvalGradY(double x, double y) = 0;

	virtual int GetDegree() = 0;
	virtual RefInterval ReferenceInterval() = 0;
	virtual string ToString() = 0;
};