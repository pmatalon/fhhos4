#pragma once
#include <string>
using namespace std;

class RefInterval
{
public:
	static RefInterval &MinusOne_One()	{ static RefInterval interval(-1, 1); return interval; }
	static RefInterval &Zero_One()		{ static RefInterval interval( 0, 1); return interval; }
	double Left;
	double Right;
private:
	RefInterval(double left, double right)
	{
		this->Left = left;
		this->Right = right;
	}
};

class IBasisFunction1D
{
public:
	virtual double Eval(double x) = 0;
	virtual double EvalDerivative(double x) = 0;

	virtual int GetDegree() = 0;
	virtual RefInterval ReferenceInterval() = 0;
	virtual string ToString() = 0;
	virtual string ToString(string var) = 0;
};

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

class IBasisFunction3D
{
public:
	virtual double Eval(double x, double y, double z) = 0;
	virtual double EvalGradX(double x, double y, double z) = 0;
	virtual double EvalGradY(double x, double y, double z) = 0;
	virtual double EvalGradZ(double x, double y, double z) = 0;
	virtual double* Grad(double x, double y, double z) = 0;

	virtual int GetDegree() = 0;
	virtual RefInterval ReferenceInterval() = 0;
	virtual string ToString() = 0;
};