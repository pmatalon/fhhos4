#pragma once
#include "../TestCases/BiHarmonic/BiHarmonicTestCase.h"
#include "Diffusion_HHO.h"
#include "../TestCases/Diffusion/VirtualDiffusionTestCase.h"
#include "../Solver/Solver.h"
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
	ZeroMeanEnforcerFromReconstructCoeffs<Dim> _integralZeroOverDomain;
	ZeroMeanEnforcerFromBoundaryFaceCoeffs<Dim> _integralOnBoundary;
public:
	HHOParameters<Dim>* HHO;
	double IntegralSource = 0;

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
		_integralZeroOverDomain = ZeroMeanEnforcerFromReconstructCoeffs<Dim>(&_diffPb);
		_integralOnBoundary = ZeroMeanEnforcerFromBoundaryFaceCoeffs<Dim>(&_diffPb);
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

		IntegralSource = _diffPb.IntegralOverDomain(_testCase->SourceFunction);

		_integralZeroOverDomain.Setup();
		_integralOnBoundary.Setup();
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
		Vector theta = -IntegralSource / boundaryMeasure * _diffPb.ProjectOnBoundaryDiscreteSpace(Utils::ConstantFunctionOne);
		return theta;
	}

	Vector Solve1stDiffProblem(const Vector& neumann)
	{
#ifndef NDEBUG
		// Check compatibility condition
		double integralNeumann = _diffPb.IntegralOverBoundaryFromFaceCoeffs(neumann);
		assert(abs(IntegralSource + integralNeumann) < Utils::Eps);
#endif
		// Define problem
		_diffPb.ChangeSourceFunction(_testCase->SourceFunction);
		_diffPb.ChangeNeumannFunction(neumann);
		Vector& rhs = _diffPb.SetCondensedRHS();

		// Solve
		_diffPb.SystemSolution = _diffSolver->Solve(rhs);
		assert(dynamic_cast<IterativeSolver*>(_diffSolver)->IterationCount < dynamic_cast<IterativeSolver*>(_diffSolver)->MaxIterations);

#ifndef NDEBUG
		// Check that mean value = 0 (on the faces)
		double integralSkeleton = _diffPb.IntegralOverSkeletonFromFaceCoeffs(_diffPb.SystemSolution);
		assert(abs(integralSkeleton) < Utils::Eps);
#endif
		// Reconstruct the higher-order polynomial
		_diffPb.ReconstructHigherOrderApproximation(false);
		Vector lambda = std::move(_diffPb.ReconstructedSolution);

		// Even if the mean value is 0 on the faces, the mean value of the reconstructed polynomial is not necessarily 0,
		// so we enforce it:
		_integralZeroOverDomain.Enforce(lambda);
		return lambda;
	}

	Vector Solve1stDiffProblemWithZeroSource(const Vector& neumann)
	{
#ifndef NDEBUG
		// Check compatibility condition
		double integralNeumann = _diffPb.IntegralOverBoundaryFromFaceCoeffs(neumann);
		assert(abs(integralNeumann) < Utils::Eps);
#endif
		// Define problem
		_diffPb.ChangeSourceFunctionToZero();
		_diffPb.ChangeNeumannFunction(neumann);
		Vector& rhs = _diffPb.SetCondensedRHS();

		// Solve
		_diffPb.SystemSolution = _diffSolver->Solve(rhs);
		assert(dynamic_cast<IterativeSolver*>(_diffSolver)->IterationCount < dynamic_cast<IterativeSolver*>(_diffSolver)->MaxIterations);

		// Reconstruct the higher-order polynomial
		_diffPb.ReconstructHigherOrderApproximation(false);
		Vector lambda = std::move(_diffPb.ReconstructedSolution);

		// Enforce (lambda|1) = 0
		_integralZeroOverDomain.Enforce(lambda);
		return lambda;
	}

	Vector Solve2ndDiffProblem(const Vector& source, bool boundaryUnknownsOnly = false)
	{
#ifndef NDEBUG
		// Check compatibility condition: (source|1) = 0
		double integralSource = _diffPb.IntegralOverDomainFromReconstructedCoeffs(source);
		assert(abs(integralSource) < Utils::Eps);
#endif
		// Define problem
		_diffPb.ChangeSourceFunction(source);
		_diffPb.ChangeNeumannFunctionToZero();
		Vector rhs = _diffPb.SetCondensedRHS();

		// Solve
		_diffPb.SystemSolution = _diffSolver->Solve(rhs);
		assert(dynamic_cast<IterativeSolver*>(_diffSolver)->IterationCount < dynamic_cast<IterativeSolver*>(_diffSolver)->MaxIterations);

		if (boundaryUnknownsOnly)
		{
			Vector boundary = _diffPb.SystemSolution.tail(_diffPb.HHO->nBoundaryFaces * _diffPb.HHO->nFaceUnknowns); // keep only the boundary unknowns
			_integralOnBoundary.Enforce(boundary);
			return boundary;
		}
		else
		{
			/*_integralOnBoundary.Enforce(_diffPb.SystemSolution);
			_diffPb.ReconstructHigherOrderApproximation(false);
			return _diffPb.ReconstructedSolution;*/
			_diffPb.ReconstructHigherOrderApproximation();
			return std::move(_diffPb.ReconstructedSolution);
		}
	}

	~BiHarmonicMixedForm_HHO()
	{
	}
};
