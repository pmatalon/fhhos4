#pragma once
#include "TestCase.h"
#include "../Problem/DiffusionField.h"
using namespace std;

template <int Dim>
class BiHarmonicTestCase : public TestCase<Dim>
{
public:
	DomFunction SourceFunction = nullptr;
	//DiffusionField<Dim> DiffField;

	BiHarmonicTestCase()
	{}

protected:
	static double SineSource2D(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		return 4 * pow(4 * M_PI, 4) * sin(4 * M_PI * x)*sin(4 * M_PI * y);
	}
	static double SineSolution2D(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		return sin(4 * M_PI * x)*sin(4 * M_PI * y);
	}

public:
	virtual ~BiHarmonicTestCase()
	{}
};
