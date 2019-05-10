#pragma once
#include <functional>
#include "../Mesh/Point.h"

class SourceFunction
{
public:
	virtual double Eval(DomPoint p) = 0;
	virtual ~SourceFunction() {}
};

class SourceFunction1D : public SourceFunction
{
private:
	function<double(double)> _f;
public:
	SourceFunction1D(function<double(double)> f)
	{
		this->_f = f;
	}
	double Eval(DomPoint p)
	{
		return this->_f(p.X);
	}
};

class SourceFunction2D : public SourceFunction
{
private:
	function<double(double, double)> _f;
public:
	SourceFunction2D(function<double(double, double)> f)
	{
		this->_f = f;
	}
	double Eval(DomPoint p)
	{
		return this->_f(p.X, p.Y);
	}
};

class SourceFunction3D : public SourceFunction
{
private:
	function<double(double, double, double)> _f;
public:
	SourceFunction3D(function<double(double, double, double)> f)
	{
		this->_f = f;
	}
	double Eval(DomPoint p)
	{
		return this->_f(p.X, p.Y, p.Z);
	}
};