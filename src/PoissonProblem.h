#pragma once
#include "Problem.h"
#include "Utils/DiffusionPartition.h"
#include "Utils/SourceFunction.h"
using namespace std;

template <int Dim>
class PoissonProblem : public Problem<Dim>
{
protected:
	DiffusionPartition<Dim>* _diffusionPartition;
	string _rhsCode;
	SourceFunction* _sourceFunction;
public:

	PoissonProblem(Mesh<Dim>* mesh, DiffusionPartition<Dim>* diffusionPartition, string rhsCode, SourceFunction* sourceFunction, string outputDirectory)
		: Problem<Dim>(mesh, outputDirectory)
	{
		this->_diffusionPartition = diffusionPartition;
		this->_rhsCode = rhsCode; 
		this->_sourceFunction = sourceFunction;

		string heterogeneityString = "";
		if (!this->_diffusionPartition->IsHomogeneous)
		{
			char res[16];
			sprintf(res, "_heterog%g", this->_diffusionPartition->HeterogeneityRatio);
			heterogeneityString = res;
		}
		this->_fileName = "Poisson" + to_string(Dim) + "D" + this->_rhsCode + heterogeneityString + "_" + this->_mesh->FileNamePart();
	}

	void PrintPhysicalProblem() override
	{
		cout << "Problem: Poisson " << Dim << "D";
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

		cout << "Analytical solution: ";
		if (this->_rhsCode.compare("sine") == 0)
			cout << "sine function";
		else if (this->_rhsCode.compare("poly") == 0)
			cout << "polynomial function";
		else if (this->_rhsCode.compare("hetero") == 0)
			cout << "heterogeneous-specific piecewise polynomial function";
		else
			cout << "unknown";
		cout << endl;
	}
};