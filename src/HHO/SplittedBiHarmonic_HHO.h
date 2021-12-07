#pragma once
#include "../Problem/BiHarmonicProblem.h"
#include "Diffusion_HHO.h"
#include "../TestCases/VirtualDiffusionTestCase.h"
using namespace std;

template<int Dim>
class SplittedBiHarmonic_HHO : public BiHarmonicProblem<Dim>
{
private:
	bool _saveMatrixBlocks = true;
	VirtualDiffusionTestCase<Dim> _diffPb1TestCase;
	VirtualDiffusionTestCase<Dim> _diffPb2TestCase;
public:
	HHOParameters<Dim>* HHO;
	Diffusion_HHO<Dim> DiffPb1;
	Diffusion_HHO<Dim> DiffPb2;

	SplittedBiHarmonic_HHO(Mesh<Dim>* mesh, BiHarmonicTestCase<Dim>* testCase, HHOParameters<Dim>* hho, bool saveMatrixBlocks, string outputDirectory)
		: BiHarmonicProblem<Dim>(mesh, testCase, outputDirectory)
	{
		HHO = hho;
		this->_diffPb1TestCase = VirtualDiffusionTestCase<Dim>(this->_sourceFunction, DiffusionField<Dim>());
		this->DiffPb1 = Diffusion_HHO<Dim>(this->_mesh, &_diffPb1TestCase, HHO, true, saveMatrixBlocks, "");
		_saveMatrixBlocks = saveMatrixBlocks;
	}

	void Assemble(ActionsArguments actions) override
	{
	}

	void AssembleDiffPb1()
	{
		ActionsArguments diffActions;
		//diffActions.LogAssembly = true;
		this->DiffPb1.Assemble(diffActions);
	}

	void AssembleDiffPb2()
	{
		DiffPb1.ReconstructHigherOrderApproximation();

		AssembleRHSForDiffPb2(DiffPb1.ReconstructedSolution);

		this->_diffPb2TestCase = VirtualDiffusionTestCase<Dim>(this->_sourceFunction, DiffusionField<Dim>());
		this->DiffPb2 = Diffusion_HHO<Dim>(this->_mesh, &_diffPb1TestCase, HHO, true, _saveMatrixBlocks, "");
	}

private:
	void AssembleRHSForDiffPb2(const Vector& reconstructedSolutionDiff1)
	{
		ElementParallelLoop<Dim> parallelLoop(this->_mesh.Elements);
		Vector b_T = Vector(HHO->nTotalCellUnknowns);

		/*if (HHO->OrthogonalizeElemBases())
		{
			DenseMatrix Nt = this->CellReconstructMassMatrix(this->CellBasis, this->ReconstructionBasis);
		}
		else
		{

		}

		parallelLoop.Execute([this](Element<Dim>* e)
			{
				Diff_HHOElement<Dim>* element = this->DiffPb1.HHOElement(e);
				element->
			}
		):*/
	}

public:
	double L2Error(DomFunction exactSolution) override
	{

	}

	void ExportSolutionToGMSH() override
	{

	}

	void ExportErrorToGMSH(const Vector& coeffs) override
	{

	}

	~SplittedBiHarmonic_HHO()
	{
	}

};
