#pragma once
#include "Problem.h"
#include "DiffusionPartition.h"
#include "SourceFunction.h"
using namespace std;

template <int Dim>
class DiffusionProblem : public Problem<Dim>
{
protected:
	DiffusionPartition<Dim>* _diffusionPartition;
	string _rhsCode;
	DomFunction _sourceFunction;
	BoundaryConditions* _boundaryConditions;
public:

	DiffusionProblem(Mesh<Dim>* mesh, DiffusionPartition<Dim>* diffusionPartition, string rhsCode, DomFunction sourceFunction, BoundaryConditions* bc, string outputDirectory)
		: Problem<Dim>(mesh, outputDirectory)
	{
		this->_diffusionPartition = diffusionPartition;
		this->_rhsCode = rhsCode; 
		this->_sourceFunction = sourceFunction;
		this->_boundaryConditions = bc;

		string heterogeneityString = "";
		if (rhsCode.compare("kellogg") != 0 && !this->_diffusionPartition->IsHomogeneous)
		{
			char res[32];
			sprintf(res, "_heterog%g", this->_diffusionPartition->HeterogeneityRatio);
			heterogeneityString = res;
		}
		this->_fileName = "Diffusion" + to_string(Dim) + "D" + this->_rhsCode + heterogeneityString + "_" + this->_mesh->FileNamePart();
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

		cout << "    Analytical solution: ";
		if (this->_rhsCode.compare("sine") == 0)
			cout << "sine function";
		else if (this->_rhsCode.compare("poly") == 0)
			cout << "polynomial function";
		else if (this->_rhsCode.compare("one") == 0)
			cout << "constant 1";
		else if (this->_rhsCode.compare("x") == 0)
			cout << "x";
		else if (this->_rhsCode.compare("hetero") == 0)
			cout << "heterogeneous-specific piecewise polynomial function";
		else if (this->_rhsCode.compare("kellogg") == 0)
			cout << "Kellogg";
		else
			cout << "unknown";
		cout << endl;
	}

	virtual void AssertSchemeConvergence(double l2Error) {}
};