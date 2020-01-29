#pragma once
#include <functional>
#include <vector>
#include "../Types.h"
#include "tetrahedron_keast_rule.hpp"
#include "../../Mesh/Point.h"
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
		else if (degree > 8)
		{
			cout << "Warning: the Keast quadradure rules only compute exact integrals of polynomials up to degree 8." << endl;
			degree = 8;
		}

		int maxRule = keast_rule_num();
		int rule, degreeRule;
		for (rule = 1; rule <= maxRule; rule++) {
			if (keast_degree(rule) >= degree)
				break;
		}
		rule++;
		rule = min(rule, maxRule);

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
};