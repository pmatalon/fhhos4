#include <cstdio>
#include <iostream>
#include <functional>
#include <getopt.h>
#include <Eigen/Core>
#include "DG/Poisson_DG.h"
#include "FunctionalBasis/FunctionalBasis.h"
#include "Mesh/CartesianGrid1D.h"
#include "Mesh/CartesianGrid2D.h"
#include "Mesh/CartesianGrid3D.h"
#include "Utils/Action.h"
#include "Utils/DiffusionPartition.h"
using namespace std;


void print_usage(string s, int d, double k, int n, string b, int p, bool f, int z, string a, string o) {
	cout << "--------------------------------------------------------" << endl;
	cout << "Arguments:" << endl;
	cout << "-d {1|2|3}:\t\t	space dimension (default: 1)\t--> " << d << endl;
	cout << "-k NUM:\t\t	diffusion coefficient k1 in the first part of the domain partition, while k2=1 in the second part (default: 1)\t--> " << k << endl;
	cout << "-s {sine|poly}:\t\t	analytical solution (default: sine)\t--> " << s << endl;
	cout << "\t\t\t 'sine' = sine solution" << endl;
	cout << "\t\t\t 'poly' = polynomial solution of global degree 2*d" << endl;
	cout << "-n NUM:\t\t		number of subdivisions in each spatial dimension (default: 5)\t--> " << n << endl;
	cout << "-b {monomials|legendre|bernstein}:	polynomial basis (default: monomials)\t--> " << b << endl;
	cout << "-p NUM:\t\t		max polynomial degree (default: 2)\t--> " << p << endl;
	cout << "-f:\t\t		full tensorization of the polynomials when d=2 or 3 (space Q) (default: false)\t--> " << f << endl;
	cout << "-z NUM:\t\t\t	penalization coefficient (default: -1 = automatic)\t--> " << z << endl;
	cout << "-a {e|c|m|s}+\t\t	action (default: es): " << endl;
	cout <<	"\t\t\t 'e' = extract system" << endl;
	cout << "\t\t\t 'c' = extract all components of the matrix in separate files" << endl;
	cout << "\t\t\t 'm' = extract mass matrix" << endl;
	cout << "\t\t\t 's' = solve system" << endl;
	cout << "-o PATH:\t		output directory to export files (default: ./)\t--> " << o << endl;
	cout << "--------------------------------------------------------" << endl;
}

int main(int argc, char* argv[])
{
	Eigen::initParallel();

	int dimension = 1;
	string solution = "sine";
	double kappa1 = 1;
	double kappa2 = 1;
	BigNumber n = 5;
	string basisCode = "monomials";
	int polyDegree = 2;
	bool fullTensorization = false;
	int penalizationCoefficient = -1;
	string a = "es";
	string outputDirectory = "./";

	int option = 0;
	while ((option = getopt(argc, argv, "d:k:s:n:b:p:z:a:o:f")) != -1) 
	{
		switch (option) 
		{
			case 'd': dimension = atoi(optarg);
				break;
			case 's': solution = optarg;
				break;
			case 'k': kappa1 = atof(optarg);
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
			case 'a': a = optarg;
				break;
			case 'o': outputDirectory = optarg;
				break;
			default: print_usage(solution, dimension, kappa1, n, basisCode, polyDegree, fullTensorization, penalizationCoefficient, a, outputDirectory);
				exit(EXIT_FAILURE);
		}
	}
	print_usage(solution, dimension, kappa1, n, basisCode, polyDegree, fullTensorization, penalizationCoefficient, a, outputDirectory);

	Action action = Action::None;
	for (int i = 0; i < a.length(); i++)
	{
		if (a[i] == 'e')
			action |= Action::ExtractSystem;
		else if (a[i] == 'c')
			action |= Action::ExtractComponentMatrices;
		else if (a[i] == 'm')
			action |= Action::ExtractMassMatrix;
		else if (a[i] == 's')
			action |= Action::SolveSystem;
	}

	Mesh* mesh;
	SourceFunction* sourceFunction;

	function<bool(Point)> isInPart1 = [](Point p) { return p.X < 0.5; };
	DiffusionPartition diffusionPartition(isInPart1, kappa1, kappa2);

	//------------//
	//     1D     //
	//------------//

	if (dimension == 1)
	{
		mesh = new CartesianGrid1D(n);

		function<double(Point)> exactSolution = [](Point p)
		{
			double x = p.X;
			return sin(4 * M_PI * x) / (16 * pow(M_PI, 2)); 
		};
		sourceFunction = new SourceFunction1D([&diffusionPartition](double x) { return diffusionPartition.Coefficient(x) * sin(4 * M_PI * x); });
		if (solution.compare("poly") == 0)
		{
			exactSolution = [](Point p) 
			{ 
				double x = p.X;
				return x * (1 - x); 
			};
			sourceFunction = new SourceFunction1D([&diffusionPartition](double x) { return diffusionPartition.Coefficient(x) * 2; });
			//sourceFunction = [](double x) { return (-1)*(-6 * x*pow(x - 1, 3) - 3 * pow(x, 3) * (2 * x - 2) - 18 * pow(x, 2) * pow(x - 1, 2)); };
		}

		//if (diffusionPartition.Kappa1 != diffusionPartition.Kappa2)
		//{
			exactSolution = [&diffusionPartition](Point p)
			{
				double x = p.X;
				double alpha = diffusionPartition.Kappa1;
				double a1 = -1 / (2 * alpha);
				double a2 = -0.5;
				double b1 = (1 + 3 * alpha)/(2 * alpha*(1 + alpha));
				double b2 = -(alpha + 3) / (2 * (1 + alpha));
				if (diffusionPartition.IsInPart1(p))
					return 4 * a1 *pow(x, 2) + 2 * b1 * x;
				else
					return 4 * a2 * pow(x - 1, 2) + 2 * b2 * (x - 1);
			};
			sourceFunction = new SourceFunction1D([&diffusionPartition](double x) {
				/*double alpha = diffusionPartition.Kappa1;
				if (diffusionPartition.IsInPart1(x))
					return 4.0 / alpha;
				else*/
					return 4.0;
			});
		//}

		Poisson_DG<1>* problem = new Poisson_DG<1>(solution);
		FunctionalBasis<1>* basis = new FunctionalBasis<1>(basisCode, polyDegree);
		Poisson_DGTerms<1>* dg = new Poisson_DGTerms<1>(sourceFunction, basis, diffusionPartition);

		problem->Assemble(mesh, basis, dg, penalizationCoefficient, outputDirectory, action);

		if ((action & Action::SolveSystem) == Action::SolveSystem)
		{
			problem->Solve();
			double error = L2::Error<1>(mesh, basis, problem->Solution, exactSolution);
			cout << "L2 Error = " << error << endl;
		}

		delete dg;
		delete problem;
		delete basis;
	}

	//------------//
	//     2D     //
	//------------//

	else if (dimension == 2)
	{
		mesh = new CartesianGrid2D(n);

		function<double(Point)> exactSolution = [](Point p) 
		{
			double x = p.X;
			double y = p.Y;
			return sin(4 * M_PI * x)*sin(4 * M_PI * y); 
		};
		sourceFunction = new SourceFunction2D([&diffusionPartition](double x, double y) { return diffusionPartition.Coefficient(x) * 2 * pow(4 * M_PI, 2) * sin(4 * M_PI * x)*sin(4 * M_PI * y); });
		if (solution.compare("poly") == 0)
		{
			exactSolution = [](Point p)
			{
				double x = p.X;
				double y = p.Y; 
				return x*(1 - x) * y*(1 - y); 
			};
			sourceFunction = new SourceFunction2D([&diffusionPartition](double x, double y) { return diffusionPartition.Coefficient(x) * 2 * y*(1 - y) + 2 * x*(1 - x); });
		}

		if (diffusionPartition.Kappa1 != diffusionPartition.Kappa2)
		{
			exactSolution = [&diffusionPartition](Point p)
			{
				double x = p.X;
				double y = p.Y;
				if (diffusionPartition.IsInPart1(p))
					return diffusionPartition.Kappa2 * (x - 0.5) * x * y * (1 - y);
				else
					return diffusionPartition.Kappa1 * (x - 0.5) * (x - 1) * y * (1 - y);
			};
			sourceFunction = new SourceFunction2D([&diffusionPartition](double x, double y) {
				if (diffusionPartition.IsInPart1(Point(x, y)))
					return -diffusionPartition.Coefficient(x)*2 * diffusionPartition.Kappa2 * (y * (1 - y) - (x - 0.5)*x);
				else
					return -diffusionPartition.Coefficient(x)*2 * diffusionPartition.Kappa1 * (y * (1 - y) - (x - 0.5)*(x - 1));
			});
		}

		Poisson_DG<2>* problem = new Poisson_DG<2>(solution);
		FunctionalBasis<2>* basis = new FunctionalBasis<2>(basisCode, polyDegree, fullTensorization);
		Poisson_DGTerms<2>* dg = new Poisson_DGTerms<2>(sourceFunction, basis, diffusionPartition);

		problem->Assemble(mesh, basis, dg, penalizationCoefficient, outputDirectory, action);

		if ((action & Action::SolveSystem) == Action::SolveSystem)
		{
			problem->Solve();
			double error = L2::Error<2>(mesh, basis, problem->Solution, exactSolution);
			cout << "L2 Error = " << error << endl;
		}

		delete dg;
		delete problem;
		delete basis;
	}

	//------------//
	//     3D     //
	//------------//

	else if (dimension == 3)
	{
		mesh = new CartesianGrid3D(n);

		function<double(Point)> exactSolution = [](Point p) 
		{ 
			double x = p.X;
			double y = p.Y;
			double z = p.Z;
			return sin(4 * M_PI * x)*sin(4 * M_PI * y)*sin(4 * M_PI * z);
		};
		sourceFunction = new SourceFunction3D([&diffusionPartition](double x, double y, double z) {  return diffusionPartition.Coefficient(x) * 3 * pow(4 * M_PI, 2) * sin(4 * M_PI * x)*sin(4 * M_PI * y)*sin(4 * M_PI * z); });
		if (solution.compare("poly") == 0)
		{
			exactSolution = [](Point p)
			{
				double x = p.X;
				double y = p.Y;
				double z = p.Z; 
				return x * (1 - x)*y*(1 - y)*z*(1 - z);
			};
			sourceFunction = new SourceFunction3D([&diffusionPartition](double x, double y, double z) { return diffusionPartition.Coefficient(x) * 2 * (y*(1 - y)*z*(1 - z) + x * (1 - x)*z*(1 - z) + x * (1 - x)*y*(1 - y)); });
		}

		Poisson_DG<3>* problem = new Poisson_DG<3>(solution);
		FunctionalBasis<3>* basis = new FunctionalBasis<3>(basisCode, polyDegree, fullTensorization);
		Poisson_DGTerms<3>* dg = new Poisson_DGTerms<3>(sourceFunction, basis, diffusionPartition);

		problem->Assemble(mesh, basis, dg, penalizationCoefficient, outputDirectory, action);

		if ((action & Action::SolveSystem) == Action::SolveSystem)
		{
			problem->Solve();
			double error = L2::Error<3>(mesh, basis, problem->Solution, exactSolution);
			cout << "L2 Error = " << error << endl;
		}

		delete dg;
		delete problem;
		delete basis;
	}
	else
	{
		cout << "Dimension " << dimension << ", are you kidding?!";
		exit(EXIT_FAILURE);
	}

	delete sourceFunction;
	delete mesh;

	cout << "-------------------------- DONE ------------------------" << endl;
    return EXIT_SUCCESS;
}