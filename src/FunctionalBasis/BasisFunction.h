#pragma once
#include <string>
#include "../Utils/Types.h"
#include "../Geometry/Point.h"
using namespace std;

template <int Dim>
class BasisFunction
{
public:
	int LocalNumber = -1;
	virtual double         Eval(const RefPoint& p) const = 0;
	virtual DimVector<Dim> Grad(const RefPoint& p) const = 0;
	virtual int GetDegree() const = 0;
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
	double Eval(const RefPoint& p) const override
	{
		return 1;
	}
	DimVector<0> Grad(const RefPoint& p) const override
	{
		return DimVector<0>();
	}
	int GetDegree() const override
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
	virtual double Eval(double x) const = 0;
	virtual double EvalDerivative(double x) const = 0;

	DimVector<1> Grad(const RefPoint& p) const override
	{
		DimVector<1> g;
		g << EvalDerivative(p.X);
		return g;
	}
	double Eval(const RefPoint& p) const override
	{
		return Eval(p.X);
	}
protected:
	void TestIsInReferenceInterval(double x) const
	{
		//assert(abs(x) < 1.5); // x should be in [-1, 1], but apparently we need a big margin...
	}
public:
	//virtual string ToString() = 0;
	virtual string ToString(string var) = 0;
};

class IBasisFunction2D : public BasisFunction<2>
{
public:
	virtual double Eval(double x, double y) const = 0;
	virtual double EvalGradX(double x, double y) const = 0;
	virtual double EvalGradY(double x, double y) const = 0;
	double Eval(const RefPoint& p) const override
	{
		return Eval(p.X, p.Y);
	}
	DimVector<2> Grad(const RefPoint& p) const override
	{
		DimVector<2> g;
		g << EvalGradX(p.X, p.Y), EvalGradY(p.X, p.Y);
		return g;
	}
};

class IBasisFunction3D : public BasisFunction<3>
{
public:
	virtual double Eval(double x, double y, double z) const = 0;
	virtual double EvalGradX(double x, double y, double z) const = 0;
	virtual double EvalGradY(double x, double y, double z) const = 0;
	virtual double EvalGradZ(double x, double y, double z) const = 0;
	double Eval(const RefPoint& p) const override
	{
		return Eval(p.X, p.Y, p.Z);
	}
	DimVector<3> Grad(const RefPoint& p) const override
	{
		DimVector<3> g;
		g << EvalGradX(p.X, p.Y, p.Z), EvalGradY(p.X, p.Y, p.Z), EvalGradZ(p.X, p.Y, p.Z);
		return g;
	}
};