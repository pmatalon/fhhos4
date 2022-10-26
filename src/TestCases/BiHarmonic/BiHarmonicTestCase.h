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
	BoundaryConditions NeumannBC;

	DomFunction MinusLaplacianOfSolution = nullptr;
	DomFunction MinusLaplacianOfSolution_Dirichlet = nullptr;

	BiHarmonicTestCase()
	{}

	void PrintPhysicalProblem() override
	{
		cout << "Problem: BiHarmonic " << Dim << "D" << endl;
		//cout << "    Geometry           : " << this->_mesh->GeometryDescription() << endl;
		cout << "    Test case          : " << this->Description() << endl;
		cout << "    Boundary conditions: homogeneous Dirichlet and Neumann";
		cout << endl;
	}

	string FilePrefix() override
	{
		return "BiHar" + to_string(Dim) + "D_" + this->Code();
	}

protected:

public:
	virtual ~BiHarmonicTestCase()
	{}
};
