#pragma once
#include "BiHarmonicMixedForm_HHO.h"
#include "../../TestCases/BiHarmonic/BiHarmonicTestCase.h"
#include "../../TestCases/Diffusion/VirtualDiffusionTestCase.h"
#include "HigherOrderBoundary.h"
using namespace std;

template<int Dim>
class BiHarmonicMixedFormGlowinski_HHO : public BiHarmonicMixedForm_HHO<Dim>
{
private:
	Mesh<Dim>* _mesh;
	BiHarmonicTestCase<Dim>* _testCase;

	//bool _saveMatrixBlocks = true;

	DiffusionField<Dim> _diffField;
	VirtualDiffusionTestCase<Dim> _diffPbTestCase;
	Diffusion_HHO<Dim> _diffPb;

	SparseMatrix _PTranspose_Stiff;
	SparseMatrix _PTranspose_Mass;
	SparseMatrix _normalDerivativeMatrix;
	HHOParameters<Dim>* HHO;

	Vector _b_fSource;
	Vector _zeroDirichlet;
public:
	bool UseIntegrationByParts = true;

	BiHarmonicMixedFormGlowinski_HHO(Mesh<Dim>* mesh, BiHarmonicTestCase<Dim>* testCase, HHOParameters<Dim>* hho, bool useIntegrationByParts, bool saveMatrixBlocks)
	{
		_mesh = mesh;
		_testCase = testCase;
		HHO = hho;
		_diffField = DiffusionField<Dim>();
		mesh->SetDiffusionField(&_diffField);
		_diffPbTestCase = VirtualDiffusionTestCase<Dim>(testCase->SourceFunction, _diffField);
		_diffPbTestCase.BC = BoundaryConditions::HomogeneousDirichletEverywhere();
		_diffPb = Diffusion_HHO<Dim>(mesh, &_diffPbTestCase, HHO, true, true);
		//_saveMatrixBlocks = saveMatrixBlocks;
		UseIntegrationByParts = useIntegrationByParts;
	}

	Diffusion_HHO<Dim>& DiffPb() override
	{
		return _diffPb;
	}

	void Setup() override
	{
		ActionsArguments diffActions;
		diffActions.AssembleRightHandSide = false;
		diffActions.LogAssembly = true;
		_diffPb.Assemble(diffActions);

		if (UseIntegrationByParts)
		{
			_PTranspose_Stiff = _diffPb.PTranspose_Stiff();
			_PTranspose_Mass = _diffPb.PTranspose_Mass();
		}
		else
			_normalDerivativeMatrix = _diffPb.NormalDerivativeMatrix();

		_b_fSource = _diffPb.AssembleSourceTerm(_testCase->SourceFunction);
		_zeroDirichlet = Vector::Zero(HHO->nDirichletCoeffs);
	}

	Vector FindCompatibleTheta() override
	{
		return _zeroDirichlet;
	}

	// Solve problem 1 (f=source, Dirich=<dirichlet>)
	Vector Solve1stDiffProblemWithFSource(const Vector& dirichlet) override
	{
		// Define problem
		Vector b_T   = _diffPb.ComputeB_T(_b_fSource, dirichlet);
		Vector b_ndF = _diffPb.ComputeB_ndF_noNeumann(dirichlet);
		Vector rhs = _diffPb.CondensedRHS(b_T, b_ndF);

		// Solve
		Vector faceSolution = this->_diffSolver->Solve(rhs);
		this->CheckDiffSolverConvergence();

		// Reconstruct the higher-order polynomial
		return _diffPb.ReconstructHigherOrderApproximationFromFaceCoeffs(faceSolution, dirichlet, b_T);
	}

	// Solve problem 1 (f=0, Dirich=<dirichlet>)
	Vector Solve1stDiffProblemWithZeroSource(const Vector& dirichlet) override
	{
		// Define problem
		Vector b_T = _diffPb.ComputeB_T_zeroSource(dirichlet);
		Vector b_ndF = _diffPb.ComputeB_ndF_noNeumann(dirichlet);
		Vector rhs = _diffPb.CondensedRHS(b_T, b_ndF);

		// Solve
		Vector faceSolution = this->_diffSolver->Solve(rhs);
		this->CheckDiffSolverConvergence();

		// Reconstruct the higher-order polynomial
		return _diffPb.ReconstructHigherOrderApproximationFromFaceCoeffs(faceSolution, dirichlet, b_T);
	}

	// Solve problem 2 (f=<source>, Dirich=0)
	Vector Solve2ndDiffProblem(const Vector& source, bool returnBoundaryNormalDerivative = false) override
	{
		// Define problem
		Vector b_source = _diffPb.AssembleSourceTerm(source);
		Vector rhs = _diffPb.CondensedRHS_noNeumannZeroDirichlet(b_source);

		// Solve
		Vector faceSolution = this->_diffSolver->Solve(rhs);
		this->CheckDiffSolverConvergence();

		if (returnBoundaryNormalDerivative)
		{
			if (UseIntegrationByParts)
			{
				/*Vector cellSolution = _diffPb.SolveCellUnknowns(faceSolution, b_source);
				Vector normalDerivative = _diffPb.A_T_dF.transpose() * cellSolution + _diffPb.A_ndF_dF.transpose() * faceSolution;
				normalDerivative -= _PTranspose_Mass * source;
				return normalDerivative;*/
				Vector reconstructedElemBoundary = _diffPb.ReconstructHigherOrderOnBoundaryOnly(faceSolution, _zeroDirichlet, b_source);
				return _PTranspose_Stiff * reconstructedElemBoundary - _PTranspose_Mass * source;
			}
			else
			{
				Vector reconstructedElemBoundary = _diffPb.ReconstructHigherOrderOnBoundaryOnly(faceSolution, _zeroDirichlet, b_source);
				return _normalDerivativeMatrix * reconstructedElemBoundary;
			}
		}
		else
			return _diffPb.ReconstructHigherOrderApproximationFromFaceCoeffs(faceSolution, _zeroDirichlet, b_source);
	}

	double L2InnerProdOnBoundary(const Vector& v1, const Vector& v2) override
	{
		return _diffPb.BoundarySpace.L2InnerProd(v1, v2);
	}

	Vector ComputeSolution(const Vector& theta) override
	{
		// Solve problem 1 (f=source, Dirich=<theta>)
		Vector lambda = Solve1stDiffProblemWithFSource(theta);

		// Solve problem 2 (f=<lambda>, Dirich=0)
		return Solve2ndDiffProblem(lambda);
	}
};
