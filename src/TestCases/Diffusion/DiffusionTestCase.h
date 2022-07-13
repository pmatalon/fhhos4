#pragma once
#include "../TestCase.h"
#include "DiffusionField.h"
using namespace std;

template <int Dim>
class DiffusionTestCase : public TestCase<Dim>
{
public:
	DomFunction SourceFunction = nullptr;
	DiffusionField<Dim> DiffField;
	BoundaryConditions BC;

	// Normal derivative on the boundary
	DomFunction ExactSolution_Neumann = nullptr;

	DiffusionTestCase()
	{}

	void PrintPhysicalProblem()
	{
		cout << "Problem: Diffusion " << Dim << "D";
		if (DiffField.IsHomogeneous && DiffField.IsIsotropic)
			cout << " (homogeneous and isotropic)" << endl;
		else
		{
			cout << endl;
			if (DiffField.IsHomogeneous)
			{
				cout << "    Homogeneous coefficient" << endl;
				cout << "    Anisotropic: ratio = " << DiffField.LargestAnisotropyRatio() << endl;
			}
			else
			{
				cout << "    Heterogeneous coefficient: ratio = " << scientific << DiffField.LargestHeterogeneityRatio() << fixed << endl;
				if (DiffField.IsIsotropic)
					cout << "    Isotropic" << endl;
				else
					cout << "    Anisotropic: ratio = " << DiffField.LargestAnisotropyRatio() << endl;
			}
		}

		//cout << "    Geometry           : " << this->_mesh->GeometryDescription() << endl;

		cout << "    Test case          : " << this->Description() << endl;
		cout << "    Boundary conditions: " << this->BC.Description << endl;
	}

	string FilePrefix()
	{
		string heterogeneityString = "";
		if (this->Code().compare("kellogg") != 0 && !DiffField.IsHomogeneous)
		{
			char res[32];
			sprintf(res, "_heterog%g", DiffField.LargestHeterogeneityRatio());
			heterogeneityString = res;
		}
		return "Diff" + to_string(Dim) + "D_" + this->Code() + heterogeneityString;
	}

protected:
	static double SineSource2D(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		return 2 * pow(4 * M_PI, 2) * sin(4 * M_PI * x)*sin(4 * M_PI * y);
	}
	static double SineSolution2D(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		double a = 1;//anisotropyCoefficients1[0];
		double b = 1;//anisotropyCoefficients1[1];
		return 2 / (a + b) * sin(4 * M_PI * x)*sin(4 * M_PI * y);
	}
	static double SineSolution2D_Neumann(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		if (abs(x) < Utils::NumericalZero) // x = 0
			return -4 * M_PI * sin(4 * M_PI * y);
		else if (abs(x - 1) < Utils::NumericalZero) // x = 1
			return 4 * M_PI * sin(4 * M_PI * y);
		else if (abs(y) < Utils::NumericalZero) // y = 0
			return -4 * M_PI * sin(4 * M_PI * x);
		else if (abs(y - 1) < Utils::NumericalZero) // y = 1
			return 4 * M_PI * sin(4 * M_PI * x);
		Utils::FatalError("Should never happen");
		return 0;
	}

	static double PolySource2D(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		return 2 * (y*(1 - y) + x * (1 - x));
	}
	static double PolySolution2D(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		return x * (1 - x) * y*(1 - y);
	}

	static double ExpSource2D(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		return (-pow(y, 4) - 2 * x*(1 + 2 * x*y*y))*exp(x*y*y);
	}
	static double ExpSolution2D(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		return exp(x*y*y);
	}

	static double SineSource3D(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		double z = p.Z;
		return 3 * pow(4 * M_PI, 2) * sin(4 * M_PI * x)*sin(4 * M_PI * y)*sin(4 * M_PI * z);
	}
	static double SineSolution3D(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		double z = p.Z;
		return sin(4 * M_PI * x)*sin(4 * M_PI * y)*sin(4 * M_PI * z);
	}

	static double PolySource3D(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		double z = p.Z;
		return 2 * ((y*(1 - y)*z*(1 - z) + x * (1 - x)*z*(1 - z) + x * (1 - x)*y*(1 - y)));
	}
	static double PolySolution3D(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		double z = p.Z;
		return x * (1 - x)*y*(1 - y)*z*(1 - z);
	}

	static double ExpSource3D(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		double z = p.Z;
		return -(pow(y, 4)*pow(z, 6) + 2 * x*pow(z, 3) + 4 * x*x*y*y*pow(z, 6) + 6 * x*y*y*z + 9 * x*x*pow(y, 4)*pow(z, 4))*exp(x*y*y*z*z*z);
	}
	static double ExpSolution3D(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		double z = p.Z;
		return exp(x*y*y*z*z*z);
	}

	static double DiscontinuousSource(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		double z = p.Z;
		return x * x + y * y + z * z <= 0.5 ? 1.0 : 0.0;
	};

public:
	virtual ~DiffusionTestCase()
	{}
};
