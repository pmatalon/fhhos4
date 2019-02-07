#include <cstdio>
#include <iostream>
#include "Poisson_DG.h"
#include <functional>
#include <getopt.h>
#include <Eigen/Core>
#include "FunctionalBasis.h"
#include "CartesianGrid1D.h"
#include "CartesianGrid2D.h"
#include "CartesianGrid3D.h"
#include "Poisson1D_DGTerms_LocalBasis.h"
#include "Poisson2D_DGTerms_LocalBasis.h"
#include "Poisson2D_DGTerms_GlobalBasis.h"
#include "Poisson3D_DGTerms_LocalBasis.h"
using namespace std;


void print_usage(string s, int d, int n, string b, int p, bool f, int z, string o) {
	cout << "--------------------------------------------------------" << endl;
	cout << "Arguments:" << endl;
	cout << "-s {sine|poly}:\t\tsolution (default: sine)\t--> " << s << endl;
	cout << "-d {1|2|3}:\t\t	space dimension (default: 1)\t--> " << d << endl;
	cout << "-n NUM:\t\t		number of subdivisions (default: 5)\t--> " << n << endl;
	cout << "-b {monomials|legendre|bernstein}:	polynomial basis (default: monomials)\t--> " << b << endl;
	cout << "-p NUM:\t\t		max polynomial degree (default: 2)\t--> " << p << endl;
	cout << "-f:\t\t		full tensorization of the polynomials when d=2 or 3 (default: false)\t--> " << f << endl;
	cout << "-z NUM:\t\t\t	penalization coefficient (default: 100)\t--> " << z << endl;
	cout << "-e:\t\t		extract all components of the matrix in separate files" << endl;
	cout << "-m:\t\t		extract mass matrix" << endl;
	cout << "-o PATH:\t		output directory to export the system (default: ./)\t--> " << o << endl;
	cout << "--------------------------------------------------------" << endl;
}

int main(int argc, char* argv[])
{
	Eigen::initParallel();

	string solution = "sine";
	int dimension = 1;
	BigNumber n = 5;
	string basisCode = "monomials";
	int polyDegree = 2;
	bool fullTensorization = false;
	int penalizationCoefficient = 100;
	string outputDirectory = "./";
	bool extractMatrixComponents = false;
	bool extractMassMatrix = false;

	int option = 0;
	while ((option = getopt(argc, argv, "s:d:n:b:p:z:o:emf")) != -1) 
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
			case 'f': fullTensorization = true;
				break;
			case 'z': penalizationCoefficient = atoi(optarg);
				break;
			case 'o': outputDirectory = optarg;
				break;
			case 'e': extractMatrixComponents = true;
				break;
			case 'm': extractMassMatrix = true;
				break;
			default: print_usage(solution, dimension, n, basisCode, polyDegree, fullTensorization, penalizationCoefficient, outputDirectory);
				exit(EXIT_FAILURE);
		}
	}
	print_usage(solution, dimension, n, basisCode, polyDegree, fullTensorization, penalizationCoefficient, outputDirectory);

	IMesh* mesh;
	Poisson_DG* problem;
	Poisson_DGTerms* dg;
	SourceFunction* sourceFunction;

	if (dimension == 1)
	{
		mesh = new CartesianGrid1D(n);

		std::function<double(double)> exactSolution = [](double x) { return sin(4 * M_PI * x) / (16 * pow(M_PI, 2)); };
		//std::function<double(double)> sourceFunction = [](double x) { return sin(4 * M_PI * x); };
		sourceFunction = new SourceFunction1D([](double x) { return sin(4 * M_PI * x); });
		if (solution.compare("poly") == 0)
		{
			std::function<double(double)> exactSolution = [](double x) { return x * (1 - x); };
			//sourceFunction = [](double x) { return 2; };
			sourceFunction = new SourceFunction1D([](double x) { return 2; });
			//sourceFunction = [](double x) { return (-1)*(-6 * x*pow(x - 1, 3) - 3 * pow(x, 3) * (2 * x - 2) - 18 * pow(x, 2) * pow(x - 1, 2)); };
		}

		problem = new Poisson_DG(solution);
		FunctionalBasis1D* basis = new FunctionalBasis1D(basisCode, polyDegree);
		dg = new Poisson1D_DGTerms_LocalBasis(sourceFunction, basis);

		problem->Assemble(mesh, basis, dg, penalizationCoefficient, outputDirectory, extractMatrixComponents, extractMassMatrix);

		problem->Solve();
		double error = L2::Error(mesh, basis, problem->Solution, exactSolution);
		cout << "L2 Error = " << error << endl;

		delete basis;
	}
	else if (dimension == 2)
	{
		mesh = new CartesianGrid2D(n);

		std::function<double(double, double)> exactSolution = [](double x, double y) { return sin(4 * M_PI * x)*sin(4 * M_PI * y); };
		//std::function<double(double, double)> sourceFunction = [](double x, double y) { return 2 * pow(4 * M_PI, 2) * sin(4 * M_PI * x)*sin(4 * M_PI * y); };
		sourceFunction = new SourceFunction2D([](double x, double y) { return 2 * pow(4 * M_PI, 2) * sin(4 * M_PI * x)*sin(4 * M_PI * y); });
		if (solution.compare("poly") == 0)
		{
			exactSolution = [](double x, double y) { return x*(1 - x) * y*(1 - y); };
			//sourceFunction = [](double x, double y) { return 2 * y*(1 - y) + 2 * x*(1 - x); };
			sourceFunction = new SourceFunction2D([](double x, double y) { return 2 * y*(1 - y) + 2 * x*(1 - x); });
		}
		problem = new Poisson_DG(solution);

		FunctionalBasis2D* basis = new FunctionalBasis2D(basisCode, polyDegree, fullTensorization);

		dg = new Poisson2D_DGTerms_LocalBasis(sourceFunction, basis);

		problem->Assemble(mesh, basis, dg, penalizationCoefficient, outputDirectory, extractMatrixComponents, extractMassMatrix);

		problem->Solve();
		double error = L2::Error(mesh, basis, problem->Solution, exactSolution);
		cout << "L2 Error = " << error << endl;

		delete basis;
	}
	else if (dimension == 3)
	{
		mesh = new CartesianGrid3D(n);

		std::function<double(double, double, double)> exactSolution = [](double x, double y, double z) { return sin(4 * M_PI * x)*sin(4 * M_PI * y)*sin(4 * M_PI * z); };
		//std::function<double(double, double, double)> sourceFunction = [](double x, double y, double z) { return 3 * pow(4 * M_PI, 2) * sin(4 * M_PI * x)*sin(4 * M_PI * y)*sin(4 * M_PI * z); };
		sourceFunction = new SourceFunction3D([](double x, double y, double z) { return 3 * pow(4 * M_PI, 2) * sin(4 * M_PI * x)*sin(4 * M_PI * y)*sin(4 * M_PI * z); });
		if (solution.compare("poly") == 0)
		{
			exactSolution = [](double x, double y, double z) { return x * (1 - x)*y*(1 - y)*z*(1 - z); };
			//sourceFunction = [](double x, double y, double z) { return 2 * (y*(1 - y)*z*(1 - z) + x * (1 - x)*z*(1 - z) + x * (1 - x)*y*(1 - y)); };
			sourceFunction = new SourceFunction3D([](double x, double y, double z) { return 2 * (y*(1 - y)*z*(1 - z) + x * (1 - x)*z*(1 - z) + x * (1 - x)*y*(1 - y)); });
		}

		problem = new Poisson_DG(solution);

		FunctionalBasis3D* basis = new FunctionalBasis3D(basisCode, polyDegree, fullTensorization);

		dg = new Poisson3D_DGTerms_LocalBasis(sourceFunction, basis);

		problem->Assemble(mesh, basis, dg, penalizationCoefficient, outputDirectory, extractMatrixComponents, extractMassMatrix);

		problem->Solve();
		double error = L2::Error(mesh, basis, problem->Solution, exactSolution);
		cout << "L2 Error = " << error << endl;

		delete basis;
	}
	else
	{
		cout << "Dimension " << dimension << ", are you kidding?!";
		exit(EXIT_FAILURE);
	}

	delete sourceFunction;
	delete dg;
	delete problem;
	delete mesh;

	cout << "-------------------------- DONE ------------------------" << endl;
    return EXIT_SUCCESS;
}