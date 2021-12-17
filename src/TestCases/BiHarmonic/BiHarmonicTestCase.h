#pragma once
#include "../TestCase.h"
#include "../Diffusion/DiffusionField.h"
using namespace std;

template <int Dim>
class BiHarmonicTestCase : public TestCase<Dim>
{
public:
	DomFunction SourceFunction = nullptr;
	BoundaryConditions DirichletBC;
	BoundaryConditions* NeumannBC = nullptr;
	BoundaryConditions* LaplacianDirichletBC = nullptr;

	DomFunction MinusLaplacianOfSolution = nullptr;

	BiHarmonicTestCase()
	{}

	void PrintPhysicalProblem() override
	{
		cout << "Problem: BiHarmonic " << Dim << "D" << endl;
		//cout << "    Geometry           : " << this->_mesh->GeometryDescription() << endl;
		cout << "    Test case          : " << this->Description() << endl;
		cout << "    Boundary conditions: Dirichlet + ";
		if (NeumannBC)
			cout << "Neumann";
		else
			cout << "Dirichlet on laplacian";
		cout << endl;
	}

	string FilePrefix()
	{
		return "BiHar" + to_string(Dim) + "D_" + this->Code();
	}

protected:
	// With homogeneous Dirichlet B.C. for both diffusion problems
	static double SineSource2D(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		return 4 * pow(4 * M_PI, 4) * sin(4 * M_PI * x) * sin(4 * M_PI * y);
	}
	static double SineMinusLaplacianOfSolution2D(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		return 2 * pow(4 * M_PI, 2) * sin(4 * M_PI * x) * sin(4 * M_PI * y);
	}
	static double SineSolution2D(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		return sin(4 * M_PI * x) * sin(4 * M_PI * y);
	}

	// With mixed Dirichlet/Neumann B.C.
	static double PolySource2D(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		return 56400 * (1 - 10*x + 15*x*x) * pow(1 - y, 2) * pow(y, 4) + 18800 * x*x * (6 - 20 * x + 15 * x*x) * y*y* (6 - 20*y + 15*y*y) + 56400*pow(1-x, 2) * pow(x, 4) * (1 - 10*y + 15*y*y);
	}
	static double PolySolution2D(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		return 2350*pow(x, 4)*pow(x-1, 2)*pow(y, 4)*pow(y-1, 2);
	}

public:
	virtual ~BiHarmonicTestCase()
	{}
};
