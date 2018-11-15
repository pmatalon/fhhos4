#include <cstdio>
#include <iostream>
#include "Poisson1D.h"
#include "Poisson2D.h"
//#include <Eigen/Dense>
#include <functional>
#include <getopt.h>
#include "MonomialBasis1D.h"
#include "ReverseMonomialBasis1D.h"
#include "LegendreBasis1D.h"
#include "MonomialBasis2D.h"
using namespace std;


void print_usage(string s, int d, int n, string b, int p, int z, string o) {
	cout << "--------------------------------------------------------" << endl;
	cout << "Arguments:" << endl;
	cout << "-s {sine|poly}\t	solution (default: sine)\t--> " << s << endl;
	cout << "-d {1|2}:\t	space dimension (default: 1)\t--> " << d << endl;
	cout << "-n NUM:\t		number of subdivisions (default: 5)\t--> " << n << endl;
	cout << "-b {monomials|reversemonomials|legendre}:	polynomial basis (default: monomials)\t--> " << b << endl;
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

	string solution = "sine";
	int dimension = 1;
	BigNumber n = 5;
	string basisCode = "monomials";
	int polyDegree = 2;
	int penalizationCoefficient = 100;
	string outputDirectory = "./";

	int option = 0;
	while ((option = getopt(argc, argv, "s:d:n:b:p:z:o:")) != -1) 
	{
		switch (option) 
		{
			case 's': solution = optarg;
				break;
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
			default: print_usage(solution, dimension, n, basisCode, polyDegree, penalizationCoefficient, outputDirectory);
				exit(EXIT_FAILURE);
		}
	}
	print_usage(solution, dimension, n, basisCode, polyDegree, penalizationCoefficient, outputDirectory);

	if (dimension == 1)
	{
		CartesianGrid1D* grid = new CartesianGrid1D(n);

		std::function<double(double)> sourceFunction = [](double x) { return sin(4 * M_PI * x); };
		if (solution.compare("poly") == 0)
			sourceFunction = [](double x) { return 2; };

		Poisson1D* problem = new Poisson1D(solution, sourceFunction);

		FunctionalBasisWithNumbers* basis;
		if (basisCode.compare("monomials") == 0)
			basis = new MonomialBasis1D(polyDegree, grid, penalizationCoefficient, sourceFunction);
		else if (basisCode.compare("globalmonomials") == 0)
			basis = new MonomialGlobalBasis1D(polyDegree, grid, penalizationCoefficient, sourceFunction);
		else if (basisCode.compare("reversemonomials") == 0)
			basis = new ReverseMonomialBasis1D(polyDegree, grid, penalizationCoefficient, sourceFunction);
		else if (basisCode.compare("legendre") == 0)
			basis = new LegendreBasis1D(polyDegree, grid, penalizationCoefficient, sourceFunction);
		else if (basisCode.compare("globallegendre") == 0)
			basis = new GlobalLegendreBasis1D(polyDegree, grid, penalizationCoefficient, sourceFunction);
		else
		{
			cout << "Basis not managed!";
			exit(EXIT_FAILURE);
		}

		problem->DiscretizeDG(grid, basis, penalizationCoefficient, outputDirectory);
		delete problem;
		delete basis;
		delete grid;
	}
	else if (dimension == 2)
	{
		CartesianGrid2D* grid = new CartesianGrid2D(n);

		std::function<double(double, double)> sourceFunction = [](double x, double y) { return 2 * pow(4 * M_PI, 2) * sin(4 * M_PI * x)*sin(4 * M_PI * y); };
		if (solution.compare("poly") == 0)
			sourceFunction = [](double x, double y) { return 2 * y*(1 - y) + 2 * x*(1 - x); };
		Poisson2D<IBasisFunction2D>* problem = new Poisson2D<IBasisFunction2D>(solution, sourceFunction);

		FunctionalBasisWithObjects<IBasisFunction2D>* basis;
		if (basisCode.compare("monomials") == 0)
			basis = new MonomialBasis2D(polyDegree, penalizationCoefficient, sourceFunction);
		else if (basisCode.compare("globalmonomials") == 0)
			basis = new MonomialGlobalBasis2D(polyDegree, grid, penalizationCoefficient, sourceFunction);
		else
		{
			cout << "Basis not managed!";
			exit(EXIT_FAILURE);
		}

		problem->DiscretizeDG(grid, basis, penalizationCoefficient, outputDirectory);
		delete problem;
		delete basis;
		delete grid;
	}
	cout << "-------------------------- DONE ------------------------" << endl;
    return 0;
}