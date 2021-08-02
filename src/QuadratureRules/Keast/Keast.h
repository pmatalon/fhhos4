#pragma once
#include <functional>
#include <vector>
#include "../../Utils/Types.h"
#include "tetrahedron_keast_rule.hpp"
#include "../../Geometry/Point.h"
using namespace tetrahedron_keast_rule;
using namespace std;

class Keast
{
private:
	int _nPoints;
	vector<RefPoint> _points;
	vector<double> _weights;

public:
	Keast() : Keast(8)
	{}

	Keast(int degree)
	{
		if (degree == 0)
		{
			_nPoints = 1;
			_points.push_back(RefPoint());
			_weights.push_back(1);
			return;
		}
		else if (degree > MaxDegree())
		{
			Utils::Warning("The Keast quadradure rules only compute exact integrals of polynomials up to degree " + to_string(MaxDegree()) + ". Degree " + to_string(degree) + " has been requested.");
			degree = MaxDegree();
		}

		int rule, degreeRule;
		for (rule = 1; rule <= MaxRule(); rule++) {
			if (keast_degree(rule) >= degree)
				break;
		}
		rule++;
		rule = min(rule, MaxRule());

		_nPoints = keast_order_num(rule);
		vector<double> xyztab(3 * _nPoints), wtab(_nPoints);
		keast_rule(rule, _nPoints, &xyztab[0], &wtab[0]);

		_points.resize(_nPoints);
		_weights.resize(_nPoints);
		for (int i = 0; i < _nPoints; i++) {
			RefPoint p(xyztab[0 + i*3], xyztab[1 + i*3], xyztab[2 + i*3]);
			_points[i] = p;
			_weights[i] = wtab[i];
		}
	}
	
	double Quadrature(std::function<double(RefPoint)> func)
	{
		double sum = 0;
		for (int i = 0; i < _nPoints; i++)
			sum += func(_points[i]) * _weights[i];
		return sum;
	}

	inline vector<RefPoint> Points()
	{
		return _points;
	}

private:
	static constexpr int MaxRule()
	{
		return 10; // keast_rule_num();
	}

	static constexpr int MaxDegree()
	{
		return 8; // keast_degree(MaxRule());
	}
};