#pragma once
#include "../TestCases/BiHarmonic/BiHarmonicTestCase.h"
#include "Diffusion_HHO.h"
#include "../TestCases/Diffusion/VirtualDiffusionTestCase.h"
using namespace std;

template<int Dim>
class BiHarmonicMixedForm_HHO
{
private:
	Mesh<Dim>* _mesh;
	bool _saveMatrixBlocks = true;
	DiffusionField<Dim> _diffField;
	VirtualDiffusionTestCase<Dim> _diffPb1TestCase;
	Diffusion_HHO<Dim> _diffPb1;
public:
	HHOParameters<Dim>* HHO;

	BiHarmonicMixedForm_HHO(Mesh<Dim>* mesh, BiHarmonicTestCase<Dim>* testCase, HHOParameters<Dim>* hho, bool saveMatrixBlocks)
	{
		_mesh = mesh;
		HHO = hho;
		_diffField = DiffusionField<Dim>(new Tensor<Dim>());
		mesh->SetDiffusionField(&_diffField);
		_diffPb1TestCase = VirtualDiffusionTestCase<Dim>(testCase->SourceFunction, _diffField);
		_diffPb1 = Diffusion_HHO<Dim>(mesh, &_diffPb1TestCase, HHO, true, saveMatrixBlocks);
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
		DiffPb1().Assemble(diffActions);
	}

	void AssembleDiffPb2(const Vector& solutionDiffPb1, bool isOfHigherDegree)
	{
		if (isOfHigherDegree)
		{
			ElementParallelLoop<Dim> parallelLoop(this->_mesh->Elements);
			DiffPb2().B_T = Vector(HHO->nTotalCellUnknowns);

			parallelLoop.Execute([this, &solutionDiffPb1](Element<Dim>* e)
				{
					Diff_HHOElement<Dim>* element = this->DiffPb2().HHOElement(e);
					DiffPb2().B_T.segment(e->Number * HHO->nCellUnknowns, HHO->nCellUnknowns) = element->ApplyCellReconstructMassMatrix(solutionDiffPb1.segment(e->Number * HHO->nReconstructUnknowns, HHO->nReconstructUnknowns));
				}
			);

			// Static condensation
			DiffPb2().b = /*B_ndF*/ -DiffPb2().A_T_ndF.transpose() * DiffPb2().Solve_A_T_T(DiffPb2().B_T);
		}
		else
			Utils::FatalError("RHS without the reconstruction non implemented");
	}

	~BiHarmonicMixedForm_HHO()
	{
	}
};
