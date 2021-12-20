#pragma once
#include "../TestCases/BiHarmonic/BiHarmonicTestCase.h"
#include "Diffusion_HHO.h"
#include "../TestCases/Diffusion/VirtualDiffusionTestCase.h"
#include "../Solver/Solver.h"
using namespace std;

template<int Dim>
class BiHarmonicMixedForm_HHO
{
private:
	Mesh<Dim>* _mesh;
	bool _saveMatrixBlocks = true;
	DiffusionField<Dim> _diffField;
	VirtualDiffusionTestCase<Dim> _diffPbTestCase;
	Diffusion_HHO<Dim> _diffPb;
	Solver* _lapSolver = nullptr;
public:
	HHOParameters<Dim>* HHO;

	BiHarmonicMixedForm_HHO(Mesh<Dim>* mesh, BiHarmonicTestCase<Dim>* testCase, HHOParameters<Dim>* hho, bool saveMatrixBlocks)
	{
		_mesh = mesh;
		HHO = hho;
		_diffField = DiffusionField<Dim>(new Tensor<Dim>());
		mesh->SetDiffusionField(&_diffField);
		_diffPbTestCase = VirtualDiffusionTestCase<Dim>(testCase->SourceFunction, _diffField);
		_diffPbTestCase.BC = BoundaryConditions::HomogeneousNeumannEverywhere();
		_diffPb = Diffusion_HHO<Dim>(mesh, &_diffPbTestCase, HHO, true, saveMatrixBlocks);
		_saveMatrixBlocks = saveMatrixBlocks;
	}

	Diffusion_HHO<Dim>& DiffPb()
	{
		return _diffPb;
	}

	void AssembleDiffPb()
	{
		ActionsArguments diffActions;
		//diffActions.LogAssembly = true;
		_diffPb.Assemble(diffActions);
	}

	void SetLaplacianSolver(Solver* solver)
	{
		_lapSolver = solver;
	}

	void Solve()
	{

	}

	~BiHarmonicMixedForm_HHO()
	{
	}
};
