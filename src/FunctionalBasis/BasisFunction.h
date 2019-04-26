#pragma once
#include <string>
#include "../Mesh/Point.h"
using namespace std;

template <int Dim>
class BasisFunction
{
public:
	int LocalNumber = -1;
	virtual double Eval(Point p) = 0;
	virtual vector<double> Grad(Point p) = 0;
	virtual int GetDegree() = 0;
	virtual string ToString() = 0;
};

class IBasisFunction1D : public BasisFunction<1>
{
public:
	virtual double Eval(double x) = 0;
	virtual double EvalDerivative(double x) = 0;

	vector<double> Grad(Point p) override
	{
		return this->Grad(p.X);
	}
	vector<double> Grad(double x)
	{
		return vector<double>{ EvalDerivative(x) };
	}
	double Eval(Point p) override
	{
		return Eval(p.X);
	}

	void TestIsInReferenceInterval(double x)
	{
		assert(x >= -1.0000000000001 && x <= 1.0000000000001);
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
	vector<double> Grad(double x, double y)
	{
		return vector<double>{ EvalGradX(x, y), EvalGradY(x, y) };
	}
	vector<double> Grad(Point p) override
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
	vector<double> Grad(double x, double y, double z)
	{
		return vector<double>{ EvalGradX(x, y, z), EvalGradY(x, y, z), EvalGradZ(x, y, z) };
	}
	vector<double> Grad(Point p) override
	{
		return this->Grad(p.X, p.Y, p.Z);
	}
};