#include <cstdio>
#include <iostream>
#include "Poisson1D.h"
#include "Poisson2D.h"
//#include <Eigen/Dense>
#include <functional>
#include <getopt.h>
#include "FunctionalBasis.h"
#include "CartesianGrid2D.h"
#include "CartesianGrid3D.h"
#include "Poisson1D_DGTerms_LocalBasis.h"
#include "Poisson1D_DGTerms_GlobalBasis.h"
#include "Poisson2D_DGTerms_LocalBasis.h"
#include "Poisson2D_DGTerms_GlobalBasis.h"
#include "Poisson3D_DGTerms_LocalBasis.h"
using namespace std;


void print_usage(string s, int d, int n, string b, int p, int z, string o) {
	cout << "--------------------------------------------------------" << endl;
	cout << "Arguments:" << endl;
	cout << "-s {sine|poly}\t	solution (default: sine)\t--> " << s << endl;
	cout << "-d {1|2|3}:\t	space dimension (default: 1)\t--> " << d << endl;
	cout << "-n NUM:\t		number of subdivisions (default: 5)\t--> " << n << endl;
	cout << "-b {monomials|legendre|bernstein}:	polynomial basis (default: monomials)\t--> " << b << endl;
	cout << "-p NUM:\t		max polynomial degree (default: 2)\t--> " << p << endl;
	cout << "-z NUM:\t		penalization coefficient (default: 100)\t--> " << z << endl;
	cout << "-o PATH:\t		output directory to export the system (default: ./)\t--> " << o << endl;
	cout << "-e:\t\t		extract all components of the matrix in separate files" << endl;
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
	bool extractMatrixComponents = false;

	int option = 0;
	while ((option = getopt(argc, argv, "s:d:n:b:p:z:o:e")) != -1) 
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
			case 'e': extractMatrixComponents = true;
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
			//sourceFunction = [](double x) { return (-1)*(-6 * x*pow(x - 1, 3) - 3 * pow(x, 3) * (2 * x - 2) - 18 * pow(x, 2) * pow(x - 1, 2)); };

		Poisson1D* problem = new Poisson1D(solution, sourceFunction);

		IPoisson1D_DGTerms* dg = new Poisson1D_DGTerms_LocalBasis(grid, sourceFunction);

		FunctionalBasisWithNumbers* basis = new FunctionalBasis1D(basisCode, polyDegree);
		/*if (basisCode.compare("monomials") == 0)
			basis = new MonomialBasis1D(polyDegree);
		else if (basisCode.compare("globalmonomials") == 0)
		{
			basis = new MonomialBasis1D(polyDegree);
			dg = new Poisson1D_DGTerms_GlobalBasis(grid, sourceFunction);
		}
		else if (basisCode.compare("legendre") == 0)
			basis = new LegendreBasis1D(polyDegree);
		else if (basisCode.compare("bernstein") == 0)
			basis = new BernsteinBasis1D(polyDegree);
		else if (basisCode.compare("bernstein2") == 0)
			basis = new Bernstein2Basis1D(polyDegree);
		else
		{
			cout << "Basis not managed!";
			exit(EXIT_FAILURE);
		}*/

		problem->DiscretizeDG(grid, basis, dg, penalizationCoefficient, outputDirectory, extractMatrixComponents);
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
		Poisson2D<IBasisFunction2D>* problem = new Poisson2D<IBasisFunction2D>(solution);

		IPoisson_DGTerms<IBasisFunction2D>* dg = new Poisson2D_DGTerms_LocalBasis(sourceFunction);

		FunctionalBasisWithObjects<IBasisFunction2D>* basis = new FunctionalBasis2D(basisCode, polyDegree);
		/*if (basisCode.compare("monomials") == 0)
			basis = new MonomialBasis2D(polyDegree);
		else if (basisCode.compare("globalmonomials") == 0)
		{
			basis = new MonomialBasis2D(polyDegree);
			dg = new Poisson2D_DGTerms_GlobalBasis(sourceFunction);
		}
		else if (basisCode.compare("legendre") == 0)
			basis = new LegendreBasis2D(polyDegree);
		else if (basisCode.compare("bernstein") == 0)
			basis = new BernsteinBasis2D(polyDegree);
		else if (basisCode.compare("bernstein2") == 0)
			basis = new Bernstein2Basis2D(polyDegree);
		else
		{
			cout << "Basis not managed!";
			exit(EXIT_FAILURE);
		}*/

		problem->DiscretizeDG(grid, basis, dg, penalizationCoefficient, outputDirectory, extractMatrixComponents);
		delete problem;
		delete basis;
		delete grid;
	}
	else if (dimension == 3)
	{
		CartesianGrid3D* grid = new CartesianGrid3D(n);

		std::function<double(double, double, double)> sourceFunction = [](double x, double y, double z) { return 3 * pow(4 * M_PI, 2) * sin(4 * M_PI * x)*sin(4 * M_PI * y)*sin(4 * M_PI * z); };
		if (solution.compare("poly") == 0)
			sourceFunction = [](double x, double y, double z) { return 2 * (y*(1 - y)*z*(1 - z) + x * (1 - x)*z*(1 - z) + x * (1 - x)*y*(1 - y)); };

		Poisson2D<IBasisFunction3D>* problem = new Poisson2D<IBasisFunction3D>(solution);

		IPoisson_DGTerms<IBasisFunction3D>* dg = new Poisson3D_DGTerms_LocalBasis(sourceFunction);

		FunctionalBasisWithObjects<IBasisFunction3D>* basis = new FunctionalBasis3D(basisCode, polyDegree);

		problem->DiscretizeDG(grid, basis, dg, penalizationCoefficient, outputDirectory, extractMatrixComponents);
		delete problem;
		delete basis;
		delete grid;
	}
	cout << "-------------------------- DONE ------------------------" << endl;
    return 0;
}