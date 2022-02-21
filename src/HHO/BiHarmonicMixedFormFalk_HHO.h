#pragma once
#include "BiHarmonicMixedForm_HHO.h"
#include "../TestCases/BiHarmonic/BiHarmonicTestCase.h"
#include "Diffusion_HHO.h"
#include "../TestCases/Diffusion/VirtualDiffusionTestCase.h"
#include "../Solver/Solver.h"
#include "HigherOrderBoundary.h"
#include "ZeroMeanEnforcer.h"
#include "NumericImageEnforcer.h"
using namespace std;

template<int Dim>
class BiHarmonicMixedFormFalk_HHO : public BiHarmonicMixedForm_HHO<Dim>
{
private:
	Mesh<Dim>* _mesh;
	BiHarmonicTestCase<Dim>* _testCase;

	bool _saveMatrixBlocks = true;
	bool _reconstructHigherOrderBoundary = false;
	//bool _enforceDirichletBCInLastPb = true;

	DiffusionField<Dim> _diffField;
	VirtualDiffusionTestCase<Dim> _diffPbTestCase;
	Diffusion_HHO<Dim> _diffPb;

	HigherOrderBoundary<Dim> _higherOrderBoundary;
	HHOBoundarySpace<Dim>* _boundarySpace = nullptr;

	ZeroMeanEnforcer _integralZeroOnDomain;
	ZeroMeanEnforcer _integralZeroOnBoundary;
	NumericImageEnforcerFromFaceCoeffs<Dim> _imageEnforcer;

	double _integralSource = 0;
	HHOParameters<Dim>* HHO;
public:

	BiHarmonicMixedFormFalk_HHO(Mesh<Dim>* mesh, BiHarmonicTestCase<Dim>* testCase, HHOParameters<Dim>* hho, bool reconstructHigherOrderBoundary, /*bool enforceDirichletBCInLastPb,*/ bool saveMatrixBlocks)
	{
		_mesh = mesh;
		_testCase = testCase;
		HHO = hho;
		_diffField = DiffusionField<Dim>();
		mesh->SetDiffusionField(&_diffField);
		_diffPbTestCase = VirtualDiffusionTestCase<Dim>(testCase->SourceFunction, _diffField);
		_diffPbTestCase.BC = BoundaryConditions::HomogeneousNeumannEverywhere();
		_diffPb = Diffusion_HHO<Dim>(mesh, &_diffPbTestCase, HHO, true, saveMatrixBlocks);
		_saveMatrixBlocks = saveMatrixBlocks;
		_reconstructHigherOrderBoundary = reconstructHigherOrderBoundary;
		//_enforceDirichletBCInLastPb = enforceDirichletBCInLastPb;
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

		_integralSource = _diffPb.IntegralOverDomain(_testCase->SourceFunction);

		if (_reconstructHigherOrderBoundary)
		{
			_higherOrderBoundary = HigherOrderBoundary<Dim>(&_diffPb);
			_higherOrderBoundary.Setup(true, false);
			_boundarySpace = &_higherOrderBoundary.BoundarySpace;
		}
		else
			_boundarySpace = &_diffPb.BoundarySpace;

		_integralZeroOnBoundary = ZeroMeanEnforcer(_boundarySpace);
		_integralZeroOnBoundary.Setup();

		_integralZeroOnDomain = ZeroMeanEnforcer(&_diffPb.ReconstructSpace);
		_integralZeroOnDomain.Setup();

		_imageEnforcer = NumericImageEnforcerFromFaceCoeffs<Dim>(&_diffPb);
		_imageEnforcer.Setup();
	}

	// Find theta verifying the compatibility condition
	//               (source|1) + <theta|1> = 0
	Vector FindCompatibleTheta() override
	{
		// We define 
		//       theta := -(source|1) / |\partial \Omega| * 1
		// where |\partial \Omega| is the measure of the boundary.
		Vector theta = -_integralSource / _boundarySpace->Measure() * _boundarySpace->Project(Utils::ConstantFunctionOne);
		return theta;
	}

	// Solve problem 1 (f=source, Neum=<neumann>)
	Vector Solve1stDiffProblem(const Vector& neumann) override
	{
#ifndef NDEBUG
		// Check compatibility condition
		double integralNeumann = _boundarySpace->Integral(neumann);
		assert(abs(_integralSource + integralNeumann) < Utils::Eps);
#endif
		// Define problem
		Vector b_source = _diffPb.AssembleSourceTerm(_testCase->SourceFunction);
		Vector b_neumann = _reconstructHigherOrderBoundary ? _higherOrderBoundary.AssembleNeumannTerm(neumann) : _diffPb.AssembleNeumannTerm(neumann);
		Vector rhs = _diffPb.CondensedRHS(b_source, b_neumann);

		// Solve
		_imageEnforcer.ProjectOntoImage(rhs); // enforce numerical compatibility
		Vector faceSolution = this->_diffSolver->Solve(rhs);
		this->CheckDiffSolverConvergence();

		// Reconstruct the higher-order polynomial
		Vector lambda = _diffPb.ReconstructHigherOrderApproximationFromFaceCoeffs(faceSolution, Vector(), b_source);

		// Even if the mean value is 0 on the faces, the mean value of the reconstructed polynomial is not necessarily 0,
		// so we enforce it:
		_integralZeroOnDomain.Enforce(lambda);
		return lambda;
	}

	// Solve problem 1 (f=0, Neum=<neumann>)
	Vector Solve1stDiffProblemWithZeroSource(const Vector& neumann) override
	{
#ifndef NDEBUG
		// Check compatibility condition
		assert(_integralZeroOnBoundary.Check(neumann));
#endif
		// Define problem
		Vector b_source = Vector::Zero(HHO->nTotalCellUnknowns);
		Vector b_neumann = _reconstructHigherOrderBoundary ? _higherOrderBoundary.AssembleNeumannTerm(neumann) : _diffPb.AssembleNeumannTerm(neumann);
		Vector& rhs = b_neumann;

		// Solve
		_imageEnforcer.ProjectOntoImage(rhs); // enforce numerical compatibility
		Vector faceSolution = this->_diffSolver->Solve(rhs);
		this->CheckDiffSolverConvergence();

		// Reconstruct the higher-order polynomial
		Vector lambda = _diffPb.ReconstructHigherOrderApproximationFromFaceCoeffs(faceSolution, Vector(), b_source);

		// Enforce (lambda|1) = 0
		_integralZeroOnDomain.Enforce(lambda);
		return lambda;
	}

	// Solve problem 2 (f=<source>, Neum=0)
	Vector Solve2ndDiffProblem(const Vector& source, bool returnBoundaryOnly = false) override
	{
#ifndef NDEBUG
		// Check compatibility condition: (source|1) = 0
		double integralSource = _diffPb.ReconstructSpace.Integral(source);
		assert(abs(integralSource) < Utils::Eps);
		// Probably the same thing:
		assert(_integralZeroOnDomain.Check(source));
#endif
		// Define problem
		Vector b_source = _diffPb.AssembleSourceTerm(source);
		Vector b_neumann = Vector::Zero(HHO->nTotalFaceUnknowns);
		Vector rhs = _diffPb.CondensedRHS(b_source, b_neumann);

		// Solve
		_imageEnforcer.ProjectOntoImage(rhs); // enforce numerical compatibility
		Vector faceSolution = this->_diffSolver->Solve(rhs);
		this->CheckDiffSolverConvergence();

		Vector boundary;
		if (_reconstructHigherOrderBoundary)
		{
			Vector reconstructedElemBoundary = _diffPb.ReconstructHigherOrderOnBoundaryOnly(faceSolution, Vector(), b_source);
			boundary = _higherOrderBoundary.Trace(reconstructedElemBoundary);
		}
		else
			boundary = faceSolution.tail(HHO->nBoundaryFaces * HHO->nFaceUnknowns); // keep only the boundary unknowns

		if (returnBoundaryOnly)
		{
			_integralZeroOnBoundary.Enforce(boundary);
			return boundary;
		}
		else
		{
			double orthogonalityFactor = _integralZeroOnBoundary.OrthogonalityFactor(boundary);

			Vector one_skeleton = _diffPb.SkeletonSpace.Project(Utils::ConstantFunctionOne);
			faceSolution -= orthogonalityFactor * one_skeleton;

			Vector reconstructedSolution = _diffPb.ReconstructHigherOrderApproximationFromFaceCoeffs(faceSolution, Vector(), b_source);
			return reconstructedSolution;
		}
	}

	double L2InnerProdOnBoundary(const Vector& v1, const Vector& v2) override
	{
		return _boundarySpace->L2InnerProd(v1, v2);
		//return v1.dot(v2);
	}
	/*
	//------------------------//
	//      Last problem      //
	//------------------------//

private:
	VirtualDiffusionTestCase<Dim>* _diffLastPbTestCase = nullptr;
	Diffusion_HHO<Dim>* _lastPb = nullptr;
	Solver* _lastPbSolver = nullptr;
public:
	Diffusion_HHO<Dim>* LastPb()
	{
		return _lastPb;
	}

	void SetLastPbSolver(Solver* solver)
	{
		_lastPbSolver = solver;
		IterativeSolver* iter = dynamic_cast<IterativeSolver*>(_lastPbSolver);
		if (iter)
		{
			iter->PrintIterationResults = false;
			iter->MaxIterations = 50;
		}
	}

	void SetupLastPb(Mesh<Dim>* mesh)
	{
		if (_enforceDirichletBCInLastPb)
		{
			_diffLastPbTestCase = new VirtualDiffusionTestCase<Dim>(Utils::ConstantFunctionZero, _diffField);

			auto FullDirichlet = BoundaryConditions::HomogeneousDirichletEverywhere();
			mesh->SetBoundaryConditions(&FullDirichlet, true);
			mesh->SetDiffusionField(&_diffField);

			HHOParameters<Dim>* hhoLast = new HHOParameters<Dim>(mesh, HHO->Stabilization, HHO->ReconstructionBasis, HHO->CellBasis, HHO->FaceBasis, HHO->OrthogonalizeElemBasesCode, HHO->OrthogonalizeFaceBasesCode);

			_lastPb = new Diffusion_HHO<Dim>(mesh, _diffLastPbTestCase, hhoLast, true, false);
			ActionsArguments diffActions;
			diffActions.LogAssembly = false;
			diffActions.AssembleRightHandSide = true;
			_lastPb->Assemble(diffActions);
		}
	}*/

	Vector ComputeSolution(const Vector& theta)
	{
		// Solve problem 1 (f=source, Neum=<theta>)
		Vector lambda = Solve1stDiffProblem(theta);

		// Solve problem 2 (f=<lambda>, Neum=0)
		Vector reconstructedSolution;
		//if (!_enforceDirichletBCInLastPb)
			reconstructedSolution = Solve2ndDiffProblem(lambda);
		/*else
		{
			Vector dirichletCoeffs = Vector::Zero(_lastPb->HHO->nDirichletCoeffs);
			Vector b_lambdaSource = _lastPb->AssembleSourceTerm(lambda);
			Vector b_noNeumann = Vector::Zero(_lastPb->HHO->nTotalFaceUnknowns);
			Vector rhs = _lastPb->CondensedRHS(b_lambdaSource, b_noNeumann);

			Vector faceSolution = _lastPbSolver->Solve(rhs);
			reconstructedSolution = _lastPb->ReconstructHigherOrderApproximationFromFaceCoeffs(faceSolution, dirichletCoeffs, b_lambdaSource);
		}*/
		return reconstructedSolution;
	}
};
