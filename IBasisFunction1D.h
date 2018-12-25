#include <cstdio>
#pragma once
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

/*class Transformation
{
public:
	virtual double F() = 0;
	virtual double JacobianF() = 0;
};

class Transform_X_to_T_in_MinusOne_One : public Transformation
{
private:
	double _a;
	double _b;
public:
	Transform_X_to_T_in_MinusOne_One(double a, double b)
	{
		this->_a = a;
		this->_b = b;
	}

	double T2X(double t)
	{
		return (this->_b - this->_a) / 2 * t + (this->_a + this->_b) / 2;
	}
	double T2XPrime()
	{
		return (this->_b - this->_a) / 2;
	}
	double X2TPrime()
	{
		return 2 / (this->_b - this->_a);
	}
};*/



class IBasisFunction1D
{
public:
	virtual double Eval(double x) = 0;
	virtual double EvalDerivative(double x) = 0;

	virtual int GetDegree() = 0;
	virtual RefInterval ReferenceInterval() = 0;
	virtual std::string ToString() = 0;
	virtual std::string ToString(std::string var) = 0;
};