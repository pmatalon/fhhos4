#include <cstdio>
#include <iostream>
#include "Poisson1D.h"
#include <Eigen/Dense>
#include <functional>
#include <getopt.h>
using namespace std;


void print_usage(int n, int p, int z, string o) {
	cout << "--------------------------------------------------------" << endl;
	cout << "Arguments:" << endl;
	cout << "-n: number of subdivisions (default: 5)\t--> " << n << endl;
	cout << "-p: max polynomial degree (default: 2)\t--> " << p << endl;
	cout << "-z: penalization coefficient (default: 100)\t--> " << z << endl;
	cout << "-o: output directory to export the system (default: ./)\t--> " << o << endl;
	cout << "--------------------------------------------------------" << endl;
}

int main(int argc, char* argv[])
{
	/*Eigen::VectorXd v(4);
	v << 1, 2, 3, 128;
	std::cout << v.transpose() << std::endl;
    printf("hello from DGHHO!\n");*/

	int n = 5;
	int polyDegree = 2;
	int penalizationCoefficient = 100;
	string outputDirectory = "./";

	int option = 0;
	while ((option = getopt(argc, argv, "n:p:z:o:")) != -1) 
	{
		switch (option) 
		{
			case 'n': n = atoi(optarg);
				break;
			case 'p': polyDegree = atoi(optarg);
				break;
			case 'z': penalizationCoefficient = atoi(optarg);
				break;
			case 'o': outputDirectory = optarg;
				break;
			default: print_usage(n, polyDegree, penalizationCoefficient, outputDirectory);
				exit(EXIT_FAILURE);
		}
	}
	print_usage(n, polyDegree, penalizationCoefficient, outputDirectory);

	std::function<double(double)> sourceFunction = [](double x) { return sin(4 * M_PI * x); };
	Poisson1D* problem = new Poisson1D(n, sourceFunction);
	problem->DiscretizeDG(polyDegree, penalizationCoefficient, outputDirectory);
	cout << "-------------------------- DONE ------------------------" << endl;
    return 0;
}