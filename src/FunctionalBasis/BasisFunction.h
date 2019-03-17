#pragma once
#include <string>
#include "../Mesh/Point.h"
using namespace std;

template <short Dim>
class BasisFunction
{
public:
	int LocalNumber = -1;
	virtual double Eval(Point p) = 0;
	virtual double* Grad(Point p) = 0;
	virtual int GetDegree() = 0;
	virtual string ToString() = 0;
};

class IBasisFunction1D : public BasisFunction<1>
{
public:
	virtual double Eval(double x) = 0;
	virtual double EvalDerivative(double x) = 0;

	double* Grad(Point p) override
	{
		return this->Grad(p.X);
	}
	double* Grad(double x)
	{
		return new double[1]{ EvalDerivative(x) };
	}
	double Eval(Point p) override
	{
		return Eval(p.X);
	}

	virtual string ToString() = 0;
	virtual string ToString(string var) = 0;
};

class IBasisFunction2D : public BasisFunction<2>
{
public:
	virtual double Eval(double x, double y) = 0;
	virtual double EvalGradX(double x, double y) = 0;
	virtual double EvalGradY(double x, double y) = 0;
	double Eval(Point p) override
	{
		return Eval(p.X, p.Y);
	}
	double* Grad(double x, double y)
	{
		return new double[2]{ EvalGradX(x, y), EvalGradY(x, y) };
	}
	double* Grad(Point p) override
	{
		return this->Grad(p.X, p.Y);
	}
};

class IBasisFunction3D : public BasisFunction<3>
{
public:
	virtual double Eval(double x, double y, double z) = 0;
	virtual double EvalGradX(double x, double y, double z) = 0;
	virtual double EvalGradY(double x, double y, double z) = 0;
	virtual double EvalGradZ(double x, double y, double z) = 0;
	double Eval(Point p) override
	{
		return Eval(p.X, p.Y, p.Z);
	}
	double* Grad(double x, double y, double z)
	{
		return new double[3]{ EvalGradX(x, y, z), EvalGradY(x, y, z), EvalGradZ(x, y, z) };
	}
	double* Grad(Point p) override
	{
		return this->Grad(p.X, p.Y, p.Z);
	}
};