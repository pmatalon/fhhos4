#pragma once
#include "../DiffusionTestCase.h"
using namespace std;

class KelloggTestCase : public DiffusionTestCase<2>
{
public:
	KelloggTestCase(ProblemArguments pb) :
		DiffusionTestCase()
	{
		if (!Utils::Contains(pb.GeoCode, "square4quadrants"))
			Utils::FatalError("The Kellogg test case is only compatible with the 'square4quadrants' geometry.");
		
		// Diffusion field
		if (pb.HeterogeneityRatio != 1)
			Utils::Warning("The heterogeneity ratio of the Kellogg problem is fixed. Ignoring -heterog argument.");
		if (pb.AnisotropyRatio != 1)
			Utils::Warning("The Kellogg problem is isotropic. Ignoring -aniso argument.");
		if (pb.AnisotropyAngle != 0)
			Utils::Warning("The Kellogg problem is isotropic. Ignoring -aniso-angle argument.");

		double anisotropyRatio = 1;
		double anisotropyAngle = 0;
		double kappa1 = 1;
		double kappa2 = 161.4476387975881;
		double heterogeneityRatio = kappa1 / kappa2;

		this->DiffField = Square4quadrantsGeometry::DiffField(heterogeneityRatio, anisotropyRatio, anisotropyAngle);

		// Source function
		this->SourceFunction = this->Zero;

		// Boundary conditions
		if (pb.BCCode.compare("d") != 0)
			Utils::FatalError("The requested boundary conditions are not defined in this test case.");

		this->BC.DirichletFunction = this->Solution;
		if (this->DiffField.IsIsotropic)
			this->ExactSolution = this->Solution;
	}

	string Code() override
	{
		return "kellogg";
	}
	string Description() override
	{
		return "Kellogg";
	}

private:
	/*static double Source(const DomPoint& p)
	{
		return 0;
	}*/

	static double Solution(const DomPoint& p)
	{
		// Conversion from [0,1]x[0,1] to [-1,1]x[-1,1]
		double x = 2 * p.X - 1;
		double y = 2 * p.Y - 1;

		// Conversion to polar coordinates (r, t)
		double r = sqrt(x*x + y * y);
		double t = 0;
		if (x > 0 && y >= 0)
			t = atan(y / x);
		else if (x > 0 && y < 0)
			t = atan(y / x) + 2 * M_PI;
		else if (x < 0)
			t = atan(y / x) + M_PI;
		else if (x == 0 && y > 0)
			t = M_PI / 2;
		else if (x == 0 && y < 0)
			t = 3 * M_PI / 2;
		else
			assert(false);

		// Function in polar coordinates
		double eps = 0.1;
		double nu = M_PI / 4;
		double ksi = -14.9225565104455152;

		if (t >= 0 && t <= M_PI / 2)
			return pow(r, eps) * cos((M_PI / 2 - ksi)*eps) * cos((t - M_PI / 2 + nu)*eps);
		if (t >= M_PI / 2 && t <= M_PI)
			return pow(r, eps) * cos(nu*eps) * cos((t - M_PI + ksi)*eps);
		if (t >= M_PI && t <= 3 * M_PI / 2)
			return pow(r, eps) * cos(ksi*eps) * cos((t - M_PI - nu)*eps);
		if (t >= 3 * M_PI / 2 && t <= 2 * M_PI)
			return pow(r, eps) * cos((M_PI / 2 - nu)*eps) * cos((t - 3 * M_PI / 2 - ksi)*eps);

		assert(false);
		return 0.0;
	}
};
