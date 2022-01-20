#pragma once
#include "../TestCases/BiHarmonic/BiHarmonicTestCase.h"
#include "Diffusion_HHO.h"
#include "../TestCases/Diffusion/VirtualDiffusionTestCase.h"
#include "../Solver/Solver.h"
#include "HigherOrderBoundary.h"
#include "ZeroMeanEnforcer.h"
using namespace std;

template<int Dim>
class BiHarmonicMixedForm_HHO
{
private:
	Mesh<Dim>* _mesh;
	BiHarmonicTestCase<Dim>* _testCase;
	bool _saveMatrixBlocks = true;
	DiffusionField<Dim> _diffField;
	VirtualDiffusionTestCase<Dim> _diffPbTestCase;
	Diffusion_HHO<Dim> _diffPb;
	Solver* _diffSolver;
	ZeroMeanEnforcerFromReconstructCoeffs<Dim> _integralZeroOnDomain;
	ZeroMeanEnforcerFromBoundaryFaceCoeffs<Dim> _integralZeroOnBoundary;
	ZeroMeanEnforcerFromHigherOrderBoundary<Dim> _integralZeroOnHigherOrderBoundary;
	HigherOrderBoundary<Dim> _higherOrderBoundary;
	double _integralSource = 0;
	bool _reconstructHigherOrderBoundary = true;
public:
	HHOParameters<Dim>* HHO;

	BiHarmonicMixedForm_HHO(Mesh<Dim>* mesh, BiHarmonicTestCase<Dim>* testCase, HHOParameters<Dim>* hho, bool saveMatrixBlocks)
	{
		_mesh = mesh;
		_testCase = testCase;
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
		diffActions.AssembleRightHandSide = false;
		_diffPb.Assemble(diffActions);

		_integralSource = _diffPb.IntegralOverDomain(_testCase->SourceFunction);

		_higherOrderBoundary = HigherOrderBoundary<Dim>(&_diffPb);
		_higherOrderBoundary.Setup();

		_integralZeroOnDomain = ZeroMeanEnforcerFromReconstructCoeffs<Dim>(&_diffPb);
		_integralZeroOnBoundary = ZeroMeanEnforcerFromBoundaryFaceCoeffs<Dim>(&_diffPb);
		_integralZeroOnHigherOrderBoundary = ZeroMeanEnforcerFromHigherOrderBoundary<Dim>(&_higherOrderBoundary);

		_integralZeroOnDomain.Setup();
		if (_reconstructHigherOrderBoundary)
			_integralZeroOnHigherOrderBoundary.Setup();
		else
			_integralZeroOnBoundary.Setup();
	}

	void SetDiffSolver(Solver* solver)
	{
		_diffSolver = solver;
		IterativeSolver* iter = dynamic_cast<IterativeSolver*>(_diffSolver);
		if (iter)
		{
			iter->PrintIterationResults = false;
			iter->MaxIterations = 50;
		}
	}

	// Find theta verifying the compatibility condition
	//               (source|1) + <theta|1> = 0
	Vector FindCompatibleTheta()
	{
		// We define 
		//       theta := -(source|1) / |\partial \Omega| * 1
		// where |\partial \Omega| is the measure of the boundary.
		double boundaryMeasure = _diffPb._mesh->BoundaryMeasure();
		Vector theta;
		if (_reconstructHigherOrderBoundary)
			theta = -_integralSource / boundaryMeasure * _higherOrderBoundary.ProjectOnBoundaryDiscreteSpace(Utils::ConstantFunctionOne);
		else
			theta = -_integralSource / boundaryMeasure * _diffPb.ProjectOnBoundaryDiscreteSpace(Utils::ConstantFunctionOne);
		return theta;
	}

	Vector Solve1stDiffProblem(const Vector& neumann)
	{
#ifndef NDEBUG
		// Check compatibility condition
		double integralNeumann = 0;
		if (_reconstructHigherOrderBoundary)
			integralNeumann = _higherOrderBoundary.IntegralOverBoundaryFromFaceCoeffs(neumann);
		else
			integralNeumann = _diffPb.IntegralOverBoundaryFromFaceCoeffs(neumann);

		assert(abs(_integralSource + integralNeumann) < Utils::Eps);
#endif
		// Define problem
		_diffPb.ChangeSourceFunction(_testCase->SourceFunction);
		if (_reconstructHigherOrderBoundary)
		{
			Vector b_ndF = _higherOrderBoundary.AssembleNeumannTerm(neumann);
			_diffPb.SetNeumannTerm(b_ndF);
		}
		else
			_diffPb.ChangeNeumannFunction(neumann);
		Vector& rhs = _diffPb.SetCondensedRHS();

		// Solve
		Vector faceSolution = _diffSolver->Solve(rhs);
		CheckDiffSolverConvergence();

#ifndef NDEBUG
		// Check that mean value = 0 (on the faces)
		double integralSkeleton = _diffPb.IntegralOverSkeletonFromFaceCoeffs(faceSolution);
		assert(abs(integralSkeleton) < Utils::Eps);
#endif
		// Reconstruct the higher-order polynomial
		Vector lambda = _diffPb.ReconstructHigherOrderApproximationFromFaceCoeffs(faceSolution);

		// Even if the mean value is 0 on the faces, the mean value of the reconstructed polynomial is not necessarily 0,
		// so we enforce it:
		_integralZeroOnDomain.Enforce(lambda);
		return lambda;
	}

	Vector Solve1stDiffProblemWithZeroSource(const Vector& neumann)
	{
#ifndef NDEBUG
		/*
		// Check compatibility condition
		double integralNeumann = _diffPb.IntegralOverBoundaryFromFaceCoeffs(neumann);
		assert(abs(integralNeumann) < Utils::Eps);
		// Probably the same thing:
		*/
		if (_reconstructHigherOrderBoundary)
			assert(_integralZeroOnHigherOrderBoundary.Check(neumann));
		else
			assert(_integralZeroOnBoundary.Check(neumann));
#endif
		// Define problem
		_diffPb.ChangeSourceFunctionToZero();
		if (_reconstructHigherOrderBoundary)
		{
			Vector b_ndF = _higherOrderBoundary.AssembleNeumannTerm(neumann);
			_diffPb.SetNeumannTerm(b_ndF);
		}
		else
			_diffPb.ChangeNeumannFunction(neumann);
		Vector& rhs = _diffPb.SetCondensedRHS();

		// Solve
		Vector faceSolution = _diffSolver->Solve(rhs);
		CheckDiffSolverConvergence();

		// Reconstruct the higher-order polynomial
		Vector lambda = _diffPb.ReconstructHigherOrderApproximationFromFaceCoeffs(faceSolution);

		// Enforce (lambda|1) = 0
		_integralZeroOnDomain.Enforce(lambda);
		return lambda;
	}

	Vector Solve2ndDiffProblem(const Vector& source, bool boundaryOnly = false)
	{
#ifndef NDEBUG
		// Check compatibility condition: (source|1) = 0
		double integralSource = _diffPb.IntegralOverDomainFromReconstructedCoeffs(source);
		assert(abs(integralSource) < Utils::Eps);
		// Probably the same thing:
		assert(_integralZeroOnDomain.Check(source));
#endif
		// Define problem
		_diffPb.ChangeSourceFunction(source);
		_diffPb.ChangeNeumannFunctionToZero();
		Vector& rhs = _diffPb.SetCondensedRHS();

		// Solve
		Vector faceSolution = _diffSolver->Solve(rhs);
		CheckDiffSolverConvergence();

		if (boundaryOnly)
		{
			if (_reconstructHigherOrderBoundary)
			{
				Vector reconstructedElemBoundary = _diffPb.ReconstructHigherOrderOnBoundaryOnly(faceSolution);
				Vector boundary = _higherOrderBoundary.Trace(reconstructedElemBoundary);
				_integralZeroOnHigherOrderBoundary.Enforce(boundary);
				return boundary;
			}
			else
			{
				Vector boundary = faceSolution.tail(HHO->nBoundaryFaces * HHO->nFaceUnknowns); // keep only the boundary unknowns
				_integralZeroOnBoundary.Enforce(boundary);
				return boundary;
			}
		}
		else
		{
			/*_integralZeroOnBoundary.Enforce(faceSolution);
			return _diffPb.ReconstructHigherOrderApproximationFaceCoeffs(faceSolution);*/
			return _diffPb.ReconstructHigherOrderApproximationFromFaceCoeffs(faceSolution);
		}
	}

	double L2InnerProdOnBoundary(const Vector& v1, const Vector& v2)
	{
		if (_reconstructHigherOrderBoundary)
			return _higherOrderBoundary.L2InnerProdOnBoundary(v1, v2);
		else
			return _diffPb.L2InnerProdOnBoundary(v1, v2);
		//return v1.dot(v2);
	}

private:
	void CheckDiffSolverConvergence()
	{
		IterativeSolver* iterSolver = dynamic_cast<IterativeSolver*>(_diffSolver);
		if (iterSolver)
		{
			if (iterSolver->IterationCount == iterSolver->MaxIterations)
				Utils::Warning("The diffusion solver has reached the max number of iterations (" + to_string(dynamic_cast<IterativeSolver*>(_diffSolver)->MaxIterations) + ")");
			//assert(dynamic_cast<IterativeSolver*>(_diffSolver)->IterationCount < dynamic_cast<IterativeSolver*>(_diffSolver)->MaxIterations);
		}
	}

public:
	~BiHarmonicMixedForm_HHO()
	{
	}
};
