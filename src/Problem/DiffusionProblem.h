#pragma once
#include "Problem.h"
#include "../TestCases/TestCase.h"
using namespace std;

template <int Dim>
class DiffusionProblem : public Problem<Dim>
{
protected:
	TestCase<Dim>* _testCase;
	DiffusionField<Dim>* _diffusionField;
	DomFunction _sourceFunction;
	BoundaryConditions* _boundaryConditions;
public:

	DiffusionProblem(Mesh<Dim>* mesh, TestCase<Dim>* testCase, string outputDirectory)
		: Problem<Dim>(mesh, outputDirectory)
	{
		this->_testCase = testCase;
		this->_diffusionField = &testCase->DiffField;
		this->_sourceFunction = testCase->SourceFunction;
		this->_boundaryConditions = &testCase->BC;

		string heterogeneityString = "";
		if (testCase->Code().compare("kellogg") != 0 && !this->_diffusionField->IsHomogeneous)
		{
			char res[32];
			sprintf(res, "_heterog%g", this->_diffusionField->HeterogeneityRatio);
			heterogeneityString = res;
		}
		this->_fileName = "Diff" + to_string(Dim) + "D_" + testCase->Code() + heterogeneityString;
		if (!this->_mesh->FileNamePart().empty())
			this->_fileName += "_" + this->_mesh->FileNamePart();
	}

	void PrintPhysicalProblem() override
	{
		cout << "Problem: Diffusion " << Dim << "D";
		if (this->_diffusionField->IsHomogeneous && this->_diffusionField->IsIsotropic)
			cout << " (homogeneous and isotropic)" << endl;
		else
		{
			cout << endl;
			if (this->_diffusionField->IsHomogeneous)
			{
				cout << "    Homogeneous coefficient" << endl;
				cout << "    Anisotropic: ratio = " << this->_diffusionField->K1->AnisotropyRatio << endl;
			}
			else
			{
				cout << "    Heterogeneous coefficient: ratio = " << scientific << this->_diffusionField->HeterogeneityRatio << fixed << endl;
				if (this->_diffusionField->IsIsotropic)
					cout << "    Isotropic" << endl;
				else
					cout << "    Anisotropic: ratio = " << this->_diffusionField->K1->AnisotropyRatio << endl;
			}
		}

		cout << "    Geometry           : " << this->_mesh->GeometryDescription() << endl;

		cout << "    Test case          : " << _testCase->Description() << endl;
		cout << "    Boundary conditions: " << _testCase->BC.Description << endl;
	}

	virtual void AssertSchemeConvergence(double l2Error) {}
};