#pragma once
#include <functional>
#include <vector>
#include "../Types.h"
#include "triangle_dunavant_rule.hpp"
#include "../../Mesh/Point.h"
using namespace std;

class Dunavant
{
public:
	//static int MAX_POINTS;
private:
	int _nPoints;
	vector<RefPoint> _points;
	vector<double> _weights;

public:
	//Dunavant() : Dunavant()
	//{ }

	Dunavant(int degree)
	{
		if (degree == 0)
		{
			/*_nPoints = 1;
			_points.push_back(RefPoint());
			_weights.push_back(0);
			return;*/
			degree = 1;
		}
		//else 

		/*int rule_num = dunavant_rule_num();
		int rule, order_num, deg_rule;
		for (rule = 1; rule < rule_num; rule++) {
			deg_rule = dunavant_degree(rule);
			if (deg_rule >= degree) break;
		} // for rule
		if (degree > MAX_POINTS)
			degree = MAX_POINTS;
		assert(rule != rule_num or deg_rule >= degree);

		//rule = min(rule, MAX_POINTS);
		if (rule < 20)
			rule++;

		order_num = dunavant_order_num(rule);
		std::vector<double> xytab(2 * order_num), wtab(order_num);
		_nPoints = order_num;
		dunavant_rule(rule, order_num, &xytab[0], &wtab[0]);*/

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
};

//int Dunavant::MAX_POINTS = dunavant_rule_num() - 1;