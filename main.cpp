#include <cstdio>
#include <iostream>
#include "Poisson1D.h"
#include "Poisson2D.h"
//#include <Eigen/Dense>
#include <functional>
#include <getopt.h>
#include "MonomialBasis1DOLD.h"
#include "MonomialBasis1D.h"
#include "ReverseMonomialBasis1D.h"
#include "MonomialBasis2D.h"
using namespace std;


void print_usage(int d, int n, string b, int p, int z, string o) {
	cout << "--------------------------------------------------------" << endl;
	cout << "Arguments:" << endl;
	cout << "-d {1,2}:\t	space dimension (default: 1)\t--> " << d << endl;
	cout << "-n NUM:\t		number of subdivisions (default: 5)\t--> " << n << endl;
	cout << "-b {monomials,reversemonomials,oldmonomials}:	polynomial basis (default: monomials)\t--> " << b << endl;
	cout << "-p NUM:\t		max polynomial degree (default: 2)\t--> " << p << endl;
	cout << "-z NUM:\t		penalization coefficient (default: 100)\t--> " << z << endl;
	cout << "-o PATH:\t		output directory to export the system (default: ./)\t--> " << o << endl;
	cout << "--------------------------------------------------------" << endl;
}

int main(int argc, char* argv[])
{
	/*Eigen::VectorXd v(4);
	v << 1, 2, 3, 128;
	std::cout << v.transpose() << std::endl;
    printf("hello from DGHHO!\n");*/

	int dimension = 1;
	BigNumber n = 5;
	string basisCode = "monomials";
	int polyDegree = 2;
	int penalizationCoefficient = 100;
	string outputDirectory = "./";

	int option = 0;
	while ((option = getopt(argc, argv, "d:n:b:p:z:o:")) != -1) 
	{
		switch (option) 
		{
			case 'd': dimension = atoi(optarg);
				break;
			case 'n': n = stoul(optarg, nullptr, 0);
				break; 
			case 'b': basisCode = optarg;
				break;
			case 'p': polyDegree = atoi(optarg);
				break;
			case 'z': penalizationCoefficient = atoi(optarg);
				break;
			case 'o': outputDirectory = optarg;
				break;
			default: print_usage(dimension, n, basisCode, polyDegree, penalizationCoefficient, outputDirectory);
				exit(EXIT_FAILURE);
		}
	}
	print_usage(dimension, n, basisCode, polyDegree, penalizationCoefficient, outputDirectory);

	if (dimension == 1)
	{
		CartesianGrid1D* grid = new CartesianGrid1D(n);

		std::function<double(double)> sourceFunction = [](double x) { return sin(4 * M_PI * x); };
		Poisson1D* problem = new Poisson1D(sourceFunction);

		FunctionalBasisWithNumbers* basis;
		if (basisCode.compare("reversemonomials") == 0)
			basis = new ReverseMonomialBasis1D(polyDegree, grid, penalizationCoefficient, sourceFunction);
		else if (basisCode.compare("oldmonomials") == 0)
			basis = new MonomialBasis1DOLD(polyDegree, grid, penalizationCoefficient, sourceFunction);
		else
			basis = new MonomialBasis1D(polyDegree, grid, penalizationCoefficient, sourceFunction);

		problem->DiscretizeDG(grid, basis, penalizationCoefficient, outputDirectory);
		delete problem;
		delete basis;
		delete grid;
	}
	else if (dimension == 2)
	{
		CartesianGrid2D* grid = new CartesianGrid2D(n);

		//std::function<double(double, double)> sourceFunction = [](double x, double y) { return 32 * pow(M_PI,2) * sin(4 * M_PI * x)*sin(4 * M_PI * y); };
		std::function<double(double, double)> sourceFunction = [](double x, double y) { return 2 * y*(1 - y) + 2 * x*(1 - x); };
		Poisson2D* problem = new Poisson2D(sourceFunction);

		FunctionalBasisWithObjects* basis;
		if (basisCode.compare("oldmonomials") == 0)
			basis = new MonomialBasis2DOLD(polyDegree, grid, penalizationCoefficient, sourceFunction);
		else
			basis = new MonomialBasis2D(polyDegree, grid, penalizationCoefficient, sourceFunction);

		problem->DiscretizeDG(grid, basis, penalizationCoefficient, outputDirectory);
		delete problem;
		delete basis;
		delete grid;
	}
	cout << "-------------------------- DONE ------------------------" << endl;
    return 0;
}