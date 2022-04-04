#pragma once
#include "BasisFunction.h"
#include <string>
#include <math.h>
using namespace std;

class LagrangeP1 : public IBasisFunction2D
{
public:
	LagrangeP1(int node)
	{
		this->LocalNumber = node;
	}

	static string Code() { return "lagrange"; };

	int GetDegree()
	{
		return 1;
	}

	string ToString()
	{
		return this->ToString("Lagrange_1");
	}

	string ToString(string var)
	{
		return "";
	}
};

class LagrangeP1_Node1 : public LagrangeP1
{
public:

	LagrangeP1_Node1() : LagrangeP1(0) {}

	double Eval(double x, double y) override
	{
		return 1 - x - y;
	}

	double EvalGradX(double x, double y) override
	{
		return -1;
	}

	double EvalGradY(double x, double y) override
	{
		return -1;
	}
};

class LagrangeP1_Node2 : public LagrangeP1
{
public:

	LagrangeP1_Node2() : LagrangeP1(1) {}

	double Eval(double x, double y) override
	{
		return x;
	}

	double EvalGradX(double x, double y) override
	{
		return 1;
	}

	double EvalGradY(double x, double y) override
	{
		return 0;
	}
};

class LagrangeP1_Node3 : public LagrangeP1
{
public:

	LagrangeP1_Node3() : LagrangeP1(2) {}

	double Eval(double x, double y) override
	{
		return y;
	}

	double EvalGradX(double x, double y) override
	{
		return 0;
	}

	double EvalGradY(double x, double y) override
	{
		return 1;
	}
};