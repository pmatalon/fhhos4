#pragma once
#include <string>
#include "../Utils/Types.h"
using namespace std;

template <int Dim>
class BasisFunction
{
public:
	int LocalNumber = -1;
	virtual double Eval(RefPoint p) = 0;
	virtual DimVector<Dim> Grad(RefPoint p) = 0;
	virtual int GetDegree() = 0;
	virtual string ToString() = 0;
	virtual ~BasisFunction() {}
};

class BasisFunction0D : public BasisFunction<0>
{
public:
	BasisFunction0D()
	{
		this->LocalNumber = 0;
	}
	double Eval(RefPoint p) override
	{
		return 1;
	}
	virtual DimVector<0> Grad(RefPoint p) override
	{
		return DimVector<0>();
	}
	virtual int GetDegree() override
	{
		return 0;
	}
	virtual string ToString() override
	{
		return "1";
	}
};

class IBasisFunction1D : public BasisFunction<1>
{
public:
	virtual double Eval(double x) = 0;
	virtual double EvalDerivative(double x) = 0;

	DimVector<1> Grad(RefPoint p) override
	{
		return this->Grad(p.X);
	}
	DimVector<1> Grad(double x)
	{
		DimVector<1> g;
		g << EvalDerivative(x);
		return g;
	}
	double Eval(RefPoint p) override
	{
		return Eval(p.X);
	}

	void TestIsInReferenceInterval(double x)
	{
		assert(abs(x) < 1.1);
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
	double Eval(RefPoint p) override
	{
		return Eval(p.X, p.Y);
	}
	DimVector<2> Grad(double x, double y)
	{
		DimVector<2> g;
		g << EvalGradX(x, y), EvalGradY(x, y);
		return g;
	}
	DimVector<2> Grad(RefPoint p) override
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
	double Eval(RefPoint p) override
	{
		return Eval(p.X, p.Y, p.Z);
	}
	DimVector<3> Grad(double x, double y, double z)
	{
		DimVector<3> g;
		g << EvalGradX(x, y, z), EvalGradY(x, y, z), EvalGradZ(x, y, z);
		return g;
	}
	DimVector<3> Grad(RefPoint p) override
	{
		return this->Grad(p.X, p.Y, p.Z);
	}
};