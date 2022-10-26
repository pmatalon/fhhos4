#pragma once
#include "../TestCase.h"
#include "../Diffusion/DiffusionField.h"
using namespace std;

template <int Dim>
class BiHarmonicDDTestCase : public TestCase<Dim>
{
public:
	DomFunction SourceFunction = nullptr;
	BoundaryConditions DirichletBC;
	BoundaryConditions* LaplacianDirichletBC = nullptr;

	DomFunction MinusLaplacianOfSolution = nullptr;

	BiHarmonicDDTestCase()
	{}

	void PrintPhysicalProblem() override
	{
		cout << "Problem: BiHarmonicDD " << Dim << "D" << endl;
		//cout << "    Geometry           : " << this->_mesh->GeometryDescription() << endl;
		cout << "    Test case          : " << this->Description() << endl;
		cout << "    Boundary conditions: Dirichlet + Dirichlet on laplacian";
		cout << endl;
	}

	string FilePrefix() override
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

public:
	virtual ~BiHarmonicDDTestCase()
	{}
};
