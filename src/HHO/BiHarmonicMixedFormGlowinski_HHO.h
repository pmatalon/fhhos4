#pragma once
#include "BiHarmonicMixedForm_HHO.h"
#include "../TestCases/BiHarmonic/BiHarmonicTestCase.h"
#include "Diffusion_HHO.h"
#include "../TestCases/Diffusion/VirtualDiffusionTestCase.h"
#include "../Solver/Solver.h"
#include "HigherOrderBoundary.h"
using namespace std;

template<int Dim>
class BiHarmonicMixedFormGlowinski_HHO : public BiHarmonicMixedForm_HHO<Dim>
{
private:
	Mesh<Dim>* _mesh;
	BiHarmonicTestCase<Dim>* _testCase;

	bool _saveMatrixBlocks = true;
	bool _reconstructHigherOrderBoundary = false;

	DiffusionField<Dim> _diffField;
	VirtualDiffusionTestCase<Dim> _diffPbTestCase;
	Diffusion_HHO<Dim> _diffPb;

	HigherOrderBoundary<Dim> _higherOrderBoundary;
	HHOBoundarySpace<Dim>* _boundarySpace = nullptr;
	HHOParameters<Dim>* HHO;
public:

	BiHarmonicMixedFormGlowinski_HHO(Mesh<Dim>* mesh, BiHarmonicTestCase<Dim>* testCase, HHOParameters<Dim>* hho, bool reconstructHigherOrderBoundary, bool saveMatrixBlocks)
	{
		_mesh = mesh;
		_testCase = testCase;
		HHO = hho;
		_diffField = DiffusionField<Dim>();
		mesh->SetDiffusionField(&_diffField);
		_diffPbTestCase = VirtualDiffusionTestCase<Dim>(testCase->SourceFunction, _diffField);
		_diffPbTestCase.BC = BoundaryConditions::HomogeneousDirichletEverywhere();
		_diffPb = Diffusion_HHO<Dim>(mesh, &_diffPbTestCase, HHO, true, saveMatrixBlocks);
		_saveMatrixBlocks = saveMatrixBlocks;
		_reconstructHigherOrderBoundary = reconstructHigherOrderBoundary;
	}

	Diffusion_HHO<Dim>& DiffPb() override
	{
		return _diffPb;
	}

	void Setup() override
	{
		ActionsArguments diffActions;
		diffActions.AssembleRightHandSide = false;
		diffActions.LogAssembly = false;
		_diffPb.Assemble(diffActions);


		_higherOrderBoundary = HigherOrderBoundary<Dim>(&_diffPb);
		_higherOrderBoundary.Setup(false, true);
		if (_reconstructHigherOrderBoundary)
		{
			//_higherOrderBoundary = HigherOrderBoundary<Dim>(&_diffPb);
			//_higherOrderBoundary.Setup(false, true);
			_boundarySpace = &_higherOrderBoundary.BoundarySpace;
		}
		else
			_boundarySpace = &_diffPb.BoundarySpace;
	}

	Vector FindCompatibleTheta() override
	{
		if (_reconstructHigherOrderBoundary)
			return Vector::Zero(_higherOrderBoundary.HHO->nDirichletCoeffs);
		else
			return Vector::Zero(HHO->nDirichletCoeffs);
	}

	// Solve problem 1 (f=source, Dirich=<dirichlet>)
	Vector Solve1stDiffProblem(const Vector& dirichlet) override
	{
		// Define problem
		Vector b_source = _diffPb.AssembleSourceTerm(_testCase->SourceFunction);
		Vector dirichletCoeffs = _reconstructHigherOrderBoundary ? _higherOrderBoundary.AssembleDirichletTerm(dirichlet) : dirichlet;
		Vector b_noNeumann = Vector::Zero(HHO->nTotalFaceUnknowns);
		Vector b_T   = _diffPb.ComputeB_T(b_source, dirichletCoeffs);
		Vector b_ndF = _diffPb.ComputeB_ndF(b_noNeumann, dirichletCoeffs);
		Vector rhs = _diffPb.CondensedRHS(b_T, b_ndF);

		// Solve
		Vector faceSolution = this->_diffSolver->Solve(rhs);
		this->CheckDiffSolverConvergence();

		// Reconstruct the higher-order polynomial
		return _diffPb.ReconstructHigherOrderApproximationFromFaceCoeffs(faceSolution, dirichletCoeffs, b_T);
	}

	// Solve problem 1 (f=0, Dirich=<dirichlet>)
	Vector Solve1stDiffProblemWithZeroSource(const Vector& dirichlet) override
	{
		// Define problem
		Vector b_source = Vector::Zero(HHO->nTotalCellUnknowns);
		Vector dirichletCoeffs = _reconstructHigherOrderBoundary ? _higherOrderBoundary.AssembleDirichletTerm(dirichlet) : dirichlet;
		Vector b_noNeumann = Vector::Zero(HHO->nTotalFaceUnknowns);
		Vector b_T = _diffPb.ComputeB_T(b_source, dirichletCoeffs);
		Vector b_ndF = _diffPb.ComputeB_ndF(b_noNeumann, dirichletCoeffs);
		Vector rhs = _diffPb.CondensedRHS(b_T, b_ndF);

		// Solve
		Vector faceSolution = this->_diffSolver->Solve(rhs);
		this->CheckDiffSolverConvergence();

		// Reconstruct the higher-order polynomial
		return _diffPb.ReconstructHigherOrderApproximationFromFaceCoeffs(faceSolution, dirichletCoeffs, b_T);
	}

	// Solve problem 2 (f=<source>, Dirich=0)
	Vector Solve2ndDiffProblem(const Vector& source, bool returnBoundaryNormalDerivative = false) override
	{
		// Define problem
		Vector dirichletCoeffs = Vector::Zero(HHO->nDirichletCoeffs);
		Vector b_source = _diffPb.AssembleSourceTerm(source);
		Vector b_noNeumann = Vector::Zero(HHO->nTotalFaceUnknowns);
		Vector rhs = _diffPb.CondensedRHS(b_source, b_noNeumann);

		// Solve
		Vector faceSolution = this->_diffSolver->Solve(rhs);
		this->CheckDiffSolverConvergence();

		if (returnBoundaryNormalDerivative)
		{
			Vector normalDerivative;
			if (_reconstructHigherOrderBoundary)
			{
				Vector reconstructedElemBoundary = _diffPb.ReconstructHigherOrderOnBoundaryOnly(faceSolution, dirichletCoeffs, b_source);
				normalDerivative = _higherOrderBoundary.NormalDerivative(reconstructedElemBoundary);
			}
			else
			{
				//Utils::FatalError("Non implemented");
				Vector reconstructedElemBoundary = _diffPb.ReconstructHigherOrderOnBoundaryOnly(faceSolution, dirichletCoeffs, b_source);
				normalDerivative = _higherOrderBoundary.NormalDerivative(reconstructedElemBoundary);
				normalDerivative = _higherOrderBoundary.AssembleDirichletTerm(normalDerivative);
			}
			return normalDerivative;
		}
		else
			return _diffPb.ReconstructHigherOrderApproximationFromFaceCoeffs(faceSolution, dirichletCoeffs, b_source);
	}

	double L2InnerProdOnBoundary(const Vector& v1, const Vector& v2) override
	{
		return _boundarySpace->L2InnerProd(v1, v2);
		//return v1.dot(v2);
	}

	Vector ComputeSolution(const Vector& theta) override
	{
		// Solve problem 1 (f=source, Dirich=<theta>)
		Vector lambda = Solve1stDiffProblem(theta);

		// Solve problem 2 (f=<lambda>, Dirich=0)
		return Solve2ndDiffProblem(lambda);
	}
};
