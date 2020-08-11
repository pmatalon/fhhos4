#pragma once
#include <functional>
#include <vector>
#include "../Types.h"
#include "triangle_dunavant_rule.hpp"
#include "../../Geometry/Point.h"
using namespace triangle_dunavant_rule;
using namespace std;

class Dunavant
{
private:
	int _nPoints;
	vector<RefPoint> _points;
	vector<double> _weights;

public:
	Dunavant() : Dunavant(10)
	{}

	Dunavant(int degree)
	{
		if (degree == 0)
		{
			_nPoints = 1;
			_points.push_back(RefPoint());
			_weights.push_back(1);
			return;
		}

		int maxRule = dunavant_rule_num();
		int rule, degreeRule;
		for (rule = 1; rule <= maxRule; rule++) {
			if (dunavant_degree(rule) >= degree)
				break;
		}
		rule++;
		rule = min(rule, maxRule);

		_nPoints = dunavant_order_num(rule);
		vector<double> xytab(2 * _nPoints), wtab(_nPoints);
		dunavant_rule(rule, _nPoints, &xytab[0], &wtab[0]);

		_points.resize(_nPoints);
		_weights.resize(_nPoints);
		for (int i = 0; i < _nPoints; i++) {
			RefPoint p(xytab[0 + i*2], xytab[1 + i*2]);
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
};