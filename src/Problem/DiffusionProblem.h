#pragma once
#include "Problem.h"
#include "../TestCases/TestCase.h"
using namespace std;

template <int Dim>
class DiffusionProblem : public Problem<Dim>
{
protected:
	TestCase<Dim>* _testCase;
	DiffusionPartition<Dim>* _diffusionPartition;
	DomFunction _sourceFunction;
	BoundaryConditions* _boundaryConditions;
public:

	DiffusionProblem(Mesh<Dim>* mesh, TestCase<Dim>* testCase, string outputDirectory)
		: Problem<Dim>(mesh, outputDirectory)
	{
		this->_testCase = testCase;
		this->_diffusionPartition = testCase->DiffPartition;
		this->_sourceFunction = testCase->SourceFunction;
		this->_boundaryConditions = &testCase->BC;

		string heterogeneityString = "";
		if (testCase->Code().compare("kellogg") != 0 && !this->_diffusionPartition->IsHomogeneous)
		{
			char res[32];
			sprintf(res, "_heterog%g", this->_diffusionPartition->HeterogeneityRatio);
			heterogeneityString = res;
		}
		this->_fileName = "Diffusion" + to_string(Dim) + "D" + testCase->Code() + heterogeneityString + "_" + this->_mesh->FileNamePart();
	}

	void PrintPhysicalProblem() override
	{
		cout << "Problem: Diffusion " << Dim << "D";
		if (this->_diffusionPartition->IsHomogeneous && this->_diffusionPartition->IsIsotropic)
			cout << " (homogeneous and isotropic)" << endl;
		else
		{
			cout << endl;
			if (this->_diffusionPartition->IsHomogeneous)
			{
				cout << "    Homogeneous coefficient" << endl;
				cout << "    Anisotropic: ratio = " << this->_diffusionPartition->K1->AnisotropyRatio << endl;
			}
			else
			{
				cout << "    Heterogeneous coefficient: partition = " << this->_diffusionPartition->Partition << endl;
				cout << "                               ratio     = " << scientific << this->_diffusionPartition->HeterogeneityRatio << fixed << endl;
				if (this->_diffusionPartition->IsIsotropic)
					cout << "    Isotropic" << endl;
				else
					cout << "    Anisotropic: ratio = " << this->_diffusionPartition->K1->AnisotropyRatio << endl;
			}
		}

		cout << "    Geometry           : " << this->_mesh->GeometryDescription() << endl;

		cout << "    Test case          : " << _testCase->Description() << endl;
		cout << "    Boundary conditions: " << _testCase->BC.Description << endl;
	}

	virtual void AssertSchemeConvergence(double l2Error) {}
};