#pragma once
#include "../TestCases/BiHarmonicDD/BiHarmonicDDTestCase.h"
#include "Diffusion_HHO.h"
#include "../TestCases/Diffusion/VirtualDiffusionTestCase.h"
using namespace std;

template<int Dim>
class BiHarmonicDDMixedForm_HHO
{
private:
	Mesh<Dim>* _mesh;
	bool _saveMatrixBlocks = true;
	DiffusionField<Dim> _diffField;
	VirtualDiffusionTestCase<Dim> _diffPbTestCase;
	Diffusion_HHO<Dim> _diffPb1;
public:
	HHOParameters<Dim>* HHO;

	BiHarmonicDDMixedForm_HHO(Mesh<Dim>* mesh, BiHarmonicDDTestCase<Dim>* testCase, HHOParameters<Dim>* hho, bool saveMatrixBlocks)
	{
		_mesh = mesh;
		HHO = hho;
		_diffField = DiffusionField<Dim>(new Tensor<Dim>());
		mesh->SetDiffusionField(&_diffField);
		_diffPbTestCase = VirtualDiffusionTestCase<Dim>(testCase->SourceFunction, _diffField);
		_diffPb1 = Diffusion_HHO<Dim>(mesh, &_diffPbTestCase, HHO, true, saveMatrixBlocks);
		_saveMatrixBlocks = saveMatrixBlocks;
	}

	Diffusion_HHO<Dim>& DiffPb1()
	{
		return _diffPb1;
	}

	Diffusion_HHO<Dim>& DiffPb2()
	{
		return _diffPb1;
	}

	void AssembleDiffPb1()
	{
		ActionsArguments diffActions;
		//diffActions.LogAssembly = true;
		_diffPb1.Assemble(diffActions);
	}

	void AssembleDiffPb2(const Vector& solutionDiffPb1)
	{
		_diffPb1.ChangeSourceFunction(solutionDiffPb1);
		_diffPb1.SetCondensedRHS();
	}

	~BiHarmonicDDMixedForm_HHO()
	{
	}
};
