#pragma once
#include "BiHarmonicMixedForm_HHO.h"
#include "../TestCases/BiHarmonic/BiHarmonicTestCase.h"
#include "Diffusion_HHO.h"
#include "../TestCases/Diffusion/VirtualDiffusionTestCase.h"
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

	Vector _b_fSource;
	Vector _b_zeroSource;
	Vector _b_zeroNeumann;
	Vector _noDirichlet;
	Vector _one_skeleton;
public:

	BiHarmonicMixedFormFalk_HHO(Mesh<Dim>* mesh, BiHarmonicTestCase<Dim>* testCase, HHOParameters<Dim>* hho, bool reconstructHigherOrderBoundary, bool saveMatrixBlocks)
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
			_higherOrderBoundary.Setup();
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

		_b_fSource = _diffPb.AssembleSourceTerm(_testCase->SourceFunction);
		_b_zeroSource = Vector::Zero(HHO->nTotalCellUnknowns);
		_b_zeroNeumann = Vector::Zero(HHO->nTotalFaceUnknowns);
		_noDirichlet = Vector();
		_one_skeleton = _diffPb.SkeletonSpace.Project(Utils::ConstantFunctionOne);
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
	Vector Solve1stDiffProblemWithFSource(const Vector& neumann) override
	{
#ifndef NDEBUG
		// Check compatibility condition
		double integralNeumann = _boundarySpace->Integral(neumann);
		assert(abs(_integralSource + integralNeumann) < Utils::Eps);
#endif
		// Define problem
		Vector b_neumann = _reconstructHigherOrderBoundary ? _higherOrderBoundary.AssembleNeumannTerm(neumann) : _diffPb.AssembleNeumannTerm(neumann);
		Vector rhs = _diffPb.CondensedRHS(_b_fSource, b_neumann);

		// Solve
		_imageEnforcer.ProjectOntoImage(rhs); // enforce numerical compatibility
		Vector faceSolution = this->_diffSolver->Solve(rhs);
		this->CheckDiffSolverConvergence();

		// Reconstruct the higher-order polynomial
		Vector lambda = _diffPb.ReconstructHigherOrderApproximationFromFaceCoeffs(faceSolution, _noDirichlet, _b_fSource);

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
		Vector b_neumann = _reconstructHigherOrderBoundary ? _higherOrderBoundary.AssembleNeumannTerm(neumann) : _diffPb.AssembleNeumannTerm(neumann);
		Vector& rhs = b_neumann;

		// Solve
		_imageEnforcer.ProjectOntoImage(rhs); // enforce numerical compatibility
		Vector faceSolution = this->_diffSolver->Solve(rhs);
		this->CheckDiffSolverConvergence();

		// Reconstruct the higher-order polynomial
		Vector lambda = _diffPb.ReconstructHigherOrderApproximationFromFaceCoeffs(faceSolution, _noDirichlet, _b_zeroSource);

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
		Vector rhs = _diffPb.CondensedRHS(b_source, _b_zeroNeumann);

		// Solve
		_imageEnforcer.ProjectOntoImage(rhs); // enforce numerical compatibility
		Vector faceSolution = this->_diffSolver->Solve(rhs);
		this->CheckDiffSolverConvergence();

		Vector boundary;
		if (_reconstructHigherOrderBoundary)
		{
			Vector reconstructedElemBoundary = _diffPb.ReconstructHigherOrderOnBoundaryOnly(faceSolution, _noDirichlet, b_source);
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
			faceSolution -= orthogonalityFactor * _one_skeleton;
			return _diffPb.ReconstructHigherOrderApproximationFromFaceCoeffs(faceSolution, _noDirichlet, b_source);
		}
	}

	double L2InnerProdOnBoundary(const Vector& v1, const Vector& v2) override
	{
		return _boundarySpace->L2InnerProd(v1, v2);
	}

	Vector ComputeSolution(const Vector& theta)
	{
		// Solve problem 1 (f=source, Neum=<theta>)
		Vector lambda = Solve1stDiffProblemWithFSource(theta);

		// Solve problem 2 (f=<lambda>, Neum=0)
		return Solve2ndDiffProblem(lambda);
	}
};
