#include <cstdio>
#include <iostream>
#include "Poisson1D.h"
#include <Eigen/Dense>
#include <functional>


int main(int argc, char* argv[])
{
	/*Eigen::VectorXd v(4);
	v << 1, 2, 3, 128;
	std::cout << v.transpose() << std::endl;
    printf("hello from DGHHO!\n");*/
	int n = 3;
	int polyDegree = 2;
	int penalizationCoefficient = 1000;
	std::function<double(double)> sourceFunction = [](double x) { return sin(4 * M_PI * x); };
	Poisson1D* problem = new Poisson1D(n, sourceFunction);
	problem->DiscretizeDG(polyDegree, penalizationCoefficient);
    return 0;
}