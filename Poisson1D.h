#pragma once
class Poisson1D
{
private:
	int _n;
	int _polyDegree;
	double* _x;

	double* _rhs;
public:
	Poisson1D(int n, int polyDegree);
	~Poisson1D();
	void DiscretizeDG();
};

