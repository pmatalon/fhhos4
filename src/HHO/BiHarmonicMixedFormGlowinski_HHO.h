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
		_diffPb.Assemble(diffActions);

		if (_reconstructHigherOrderBoundary)
		{
			_higherOrderBoundary = HigherOrderBoundary<Dim>(&_diffPb);
			_higherOrderBoundary.Setup();
			_boundarySpace = &_higherOrderBoundary.BoundarySpace;
		}
		else
			_boundarySpace = &_diffPb.BoundarySpace;
	}

	Vector FindCompatibleTheta() override
	{
		
	}

	// Solve problem 1 (f=source, Dirich=<dirichlet>)
	Vector Solve1stDiffProblem(const Vector& dirichlet) override
	{
		// Define problem
		_diffPb.ChangeSourceFunction(_testCase->SourceFunction);
		/*if (_reconstructHigherOrderBoundary)
		{
			Vector b_ndF = _higherOrderBoundary.AssembleNeumannTerm(neumann);
			_diffPb.SetNeumannTerm(b_ndF);
		}
		else
			_diffPb.ChangeNeumannFunction(neumann);*/
		Vector& rhs = _diffPb.SetCondensedRHS();

		// Solve
		Vector faceSolution = this->_diffSolver->Solve(rhs);
		this->CheckDiffSolverConvergence();

		// Reconstruct the higher-order polynomial
		Vector lambda = _diffPb.ReconstructHigherOrderApproximationFromFaceCoeffs(faceSolution);

		return lambda;
	}

	// Solve problem 1 (f=0, Dirich=<dirichlet>)
	Vector Solve1stDiffProblemWithZeroSource(const Vector& dirichlet) override
	{
		// Define problem
		_diffPb.ChangeSourceFunctionToZero();
		/*if (_reconstructHigherOrderBoundary)
		{
			Vector b_ndF = _higherOrderBoundary.AssembleNeumannTerm(neumann);
			_diffPb.SetNeumannTerm(b_ndF);
		}
		else
			_diffPb.ChangeNeumannFunction(neumann);*/
		Vector& rhs = _diffPb.SetCondensedRHS();

		// Solve
		Vector faceSolution = this->_diffSolver->Solve(rhs);
		this->CheckDiffSolverConvergence();

		// Reconstruct the higher-order polynomial
		return _diffPb.ReconstructHigherOrderApproximationFromFaceCoeffs(faceSolution);
	}

	// Solve problem 2 (f=<source>, Dirich=0)
	Vector Solve2ndDiffProblem(const Vector& source, bool returnBoundaryOnly = false) override
	{
		// Define problem
		_diffPb.ChangeSourceFunction(source);
		_diffPb.ChangeNeumannFunctionToZero();
		Vector& rhs = _diffPb.SetCondensedRHS();

		// Solve
		Vector faceSolution = this->_diffSolver->Solve(rhs);
		this->CheckDiffSolverConvergence();

		if (returnBoundaryOnly)
		{
			Vector boundary;
			if (_reconstructHigherOrderBoundary)
			{
				Vector reconstructedElemBoundary = _diffPb.ReconstructHigherOrderOnBoundaryOnly(faceSolution);
				boundary = _higherOrderBoundary.Trace(reconstructedElemBoundary);
			}
			else
				boundary = faceSolution.tail(HHO->nBoundaryFaces * HHO->nFaceUnknowns); // keep only the boundary unknowns
			return boundary;
		}
		else
		{
			Vector reconstructedSolution = _diffPb.ReconstructHigherOrderApproximationFromFaceCoeffs(faceSolution);
			return reconstructedSolution;
		}
	}

	double L2InnerProdOnBoundary(const Vector& v1, const Vector& v2) override
	{
		return _boundarySpace->L2InnerProd(v1, v2);
		//return v1.dot(v2);
	}

	Vector ComputeSolution(const Vector& theta) override
	{
		// Solve problem 1 (f=source, Neum=<theta>)
		Vector lambda = Solve1stDiffProblem(theta);

		// Solve problem 2 (f=<lambda>, Neum=0)
		return Solve2ndDiffProblem(lambda);
	}
};
