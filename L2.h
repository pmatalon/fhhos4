#pragma once
#include "Utils.h"
#include "IBasisFunction.h"
#include <Eigen/Sparse>
#include "IMesh.h"
#include "FunctionalBasis.h"
#include "Square.h"
#include "Cube.h"
using namespace std;

class L2
{
public:

	//--------//
	//   2D   //
	//--------//

	static double Error(IMesh* mesh, FunctionalBasis2D* basis, Eigen::VectorXd solution, function<double(double, double)> exactSolution)
	{
		double absoluteError = 0;
		double normExactSolution = 0;
		for (Element* element : mesh->Elements)
		{
			Square* square = static_cast<Square*>(element);

			double x1 = square->X;
			double x2 = square->X + square->Width;
			double y1 = square->Y;
			double y2 = square->Y + square->Width;

			RefInterval refInterval = basis->GetLocalBasisFunction(square, 0)->ReferenceInterval();

			auto approximate = basis->GetApproximateFunction(solution, element->Number * basis->NumberOfLocalFunctionsInElement(element));

			function<double(double, double)> errorFunction = NULL;
			if (refInterval.Left == -1 && refInterval.Right == 1)
			{
				errorFunction = [exactSolution, approximate, x1, x2, y1, y2](double t, double u) {
					return pow(exactSolution((x2 - x1) / 2 * t + (x2 + x1) / 2, (y2 - y1) / 2 * u + (y2 + y1) / 2) - approximate(t, u), 2);
				};
			}
			else
			{
				errorFunction = [exactSolution, approximate, x1, x2, y1, y2](double t, double u) {
					return pow(exactSolution((x2 - x1) * t + x1, (y2 - y1) * u + y1) - approximate(t, u), 2);
				};
			}

			double jacobian = (x2 - x1) * (y2 - y1) / pow(refInterval.Length, 2);
			absoluteError += jacobian * Utils::Integral(errorFunction, refInterval, refInterval);

			normExactSolution += Utils::Integral([exactSolution](double x, double y) { return pow(exactSolution(x, y), 2); }, x1, x2, y1, y2);
		}
		absoluteError = sqrt(absoluteError);
		normExactSolution = sqrt(normExactSolution);
		return absoluteError / normExactSolution;
	}

	//--------//
	//   3D   //
	//--------//

	static double Error(IMesh* mesh, FunctionalBasis3D* basis, Eigen::VectorXd solution, function<double(double, double, double)> exactSolution)
	{
		double absoluteError = 0;
		double normExactSolution = 0;
		for (Element* element : mesh->Elements)
		{
			Cube* cube = static_cast<Cube*>(element);

			double x1 = cube->X;
			double x2 = cube->X + cube->Width;
			double y1 = cube->Y;
			double y2 = cube->Y + cube->Width;
			double z1 = cube->Z;
			double z2 = cube->Z + cube->Width;

			RefInterval refInterval = basis->GetLocalBasisFunction(cube, 0)->ReferenceInterval();

			auto approximate = basis->GetApproximateFunction(solution, element->Number * basis->NumberOfLocalFunctionsInElement(element));

			function<double(double, double, double)> errorFunction = NULL;
			if (refInterval.Left == -1 && refInterval.Right == 1)
			{
				errorFunction = [exactSolution, approximate, x1, x2, y1, y2, z1, z2](double t, double u, double v) {
					return pow(exactSolution((x2 - x1) / 2 * t + (x2 + x1) / 2, (y2 - y1) / 2 * u + (y2 + y1) / 2, (z2 - z1) / 2 * v + (z2 + z1) / 2) - approximate(t, u, v), 2);
				};
			}
			else
			{
				errorFunction = [exactSolution, approximate, x1, x2, y1, y2, z1, z2](double t, double u, double v) {
					return pow(exactSolution((x2 - x1) * t + x1, (y2 - y1) * u + y1, (z2 - z1) * v + z1) - approximate(t, u, v), 2);
				};
			}

			double jacobian = (x2 - x1) * (y2 - y1) * (z2 - z1) / pow(refInterval.Length, 3);
			absoluteError += jacobian * Utils::Integral(errorFunction, refInterval, refInterval, refInterval);

			normExactSolution += Utils::Integral([exactSolution](double x, double y, double z) { return pow(exactSolution(x, y, z), 2); }, x1, x2, y1, y2, z1, z2);
		}
		absoluteError = sqrt(absoluteError);
		normExactSolution = sqrt(normExactSolution);
		return absoluteError / normExactSolution;
	}
};