#pragma once
#include <cstdio>
#include "Utils.h"
#include "../Mesh/Mesh.h"
#include "../FunctionalBasis/FunctionalBasis.h"
using namespace std;

class L2
{
public:
	template <int Dim>
	static double Error(Mesh<Dim>* mesh, FunctionalBasis<Dim>* basis, Eigen::VectorXd solution, function<double(DomPoint)> exactSolution)
	{
		double absoluteError = 0;
		double normExactSolution = 0;
		for (Element<Dim>* element : mesh->Elements)
		{
			auto approximate = basis->GetApproximateFunction(solution, element->Number * basis->NumberOfLocalFunctionsInElement(element));
			absoluteError += element->L2ErrorPow2(approximate, exactSolution);
			normExactSolution += element->IntegralGlobalFunction([exactSolution](DomPoint p) { return pow(exactSolution(p), 2); });
		}
		absoluteError = sqrt(absoluteError);
		normExactSolution = sqrt(normExactSolution);
		return absoluteError / normExactSolution;
	}
};