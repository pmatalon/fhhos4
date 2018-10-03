#include <cstdio>
#include <iostream>
#include "Poisson1D.h"
#include <Eigen/Dense>


int main(int argc, char* argv[])
{
	Eigen::VectorXd v(4);
	v << 1, 2, 3, 128;
	std::cout << v.transpose() << std::endl;
    printf("hello from DGHHO!\n");
	int n = 8;
	int polyDegree = 2;
	Poisson1D* problem = new Poisson1D(n, polyDegree);
	problem->DiscretizeDG();
    return 0;
}