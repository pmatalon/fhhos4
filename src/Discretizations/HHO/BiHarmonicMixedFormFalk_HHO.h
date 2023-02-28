#pragma once
#include "BiHarmonicMixedForm_HHO.h"
#include "../../TestCases/BiHarmonic/BiHarmonicTestCase.h"
#include "../../TestCases/Diffusion/VirtualDiffusionTestCase.h"
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

	//bool _saveMatrixBlocks = true;

	DiffusionField<Dim> _diffField;
	VirtualDiffusionTestCase<Dim> _diffPbTestCase;
	Diffusion_HHO<Dim> _diffPb;

	HigherOrderBoundary<Dim> _higherOrderBoundary;
	HHOBoundarySpace<Dim>* _boundarySpace = nullptr;

	ZeroMeanEnforcer _integralZeroOnDomain;
	ZeroMeanEnforcer _integralZeroOnBoundary;
	NumericImageEnforcer _imageEnforcer;

	double _integralSource = 0;
	HHOParameters<Dim>* HHO;

	Vector _b_fSource;
	Vector _b_zeroSource;
	Vector _noDirichlet;
	Vector _one_skeleton;
public:
	bool ReconstructHigherOrderBoundary = false;

	BiHarmonicMixedFormFalk_HHO(Mesh<Dim>* mesh, BiHarmonicTestCase<Dim>* testCase, HHOParameters<Dim>* hho, bool reconstructHigherOrderBoundary, bool saveMatrixBlocks)
	{
		_mesh = mesh;
		_testCase = testCase;
		HHO = hho;
		_diffField = DiffusionField<Dim>();
		mesh->SetDiffusionField(&_diffField);
		_diffPbTestCase = VirtualDiffusionTestCase<Dim>(testCase->SourceFunction, _diffField);
		_diffPbTestCase.BC = BoundaryConditions::HomogeneousNeumannEverywhere();
		_diffPb = Diffusion_HHO<Dim>(mesh, &_diffPbTestCase, HHO, true, true);
		//_saveMatrixBlocks = saveMatrixBlocks;
		ReconstructHigherOrderBoundary = reconstructHigherOrderBoundary;
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

		_integralSource = _diffPb.CellSpace.Integral(_testCase->SourceFunction);

		if (ReconstructHigherOrderBoundary)
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

		_imageEnforcer = NumericImageEnforcer(&_diffPb.SkeletonSpace);
		_imageEnforcer.Setup();

		_b_fSource = _diffPb.AssembleSourceTerm(_testCase->SourceFunction);
		_b_zeroSource = Vector::Zero(HHO->nTotalCellUnknowns);
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
	Vector Solve1stDiffProblem(const Vector& neumann) override
	{
#ifndef NDEBUG
		// Check compatibility condition
		double integralNeumann = _boundarySpace->Integral(neumann);
		assert(abs(_integralSource + integralNeumann) < Utils::Eps);
#endif
		// Define problem
		Vector b_neumann = ReconstructHigherOrderBoundary ? _higherOrderBoundary.AssembleNeumannTerm(neumann) : _diffPb.AssembleNeumannTerm(neumann);
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
	Vector Solve1stDiffProblem_Homogeneous(const Vector& neumann)
	{
#ifndef NDEBUG
		// Check compatibility condition
		assert(_integralZeroOnBoundary.Check(neumann));
#endif
		// Define problem
		Vector b_neumann = ReconstructHigherOrderBoundary ? _higherOrderBoundary.AssembleNeumannTerm(neumann) : _diffPb.AssembleNeumannTerm(neumann);
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
		Vector rhs = _diffPb.CondensedRHS_noDirichletZeroNeumann(b_source); // TODO: Update for non-homogeneous Dirichlet BC

		// Solve
		_imageEnforcer.ProjectOntoImage(rhs); // enforce numerical compatibility
		Vector faceSolution = this->_diffSolver->Solve(rhs);
		this->CheckDiffSolverConvergence();

		Vector boundary;
		if (ReconstructHigherOrderBoundary)
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

	Vector Solve2ndDiffProblem_Homogeneous(const Vector& source)
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
		Vector rhs = _diffPb.CondensedRHS_noDirichletZeroNeumann(b_source);

		// Solve
		_imageEnforcer.ProjectOntoImage(rhs); // enforce numerical compatibility
		Vector faceSolution = this->_diffSolver->Solve(rhs);
		this->CheckDiffSolverConvergence();

		Vector boundary;
		if (ReconstructHigherOrderBoundary)
		{
			Vector reconstructedElemBoundary = _diffPb.ReconstructHigherOrderOnBoundaryOnly(faceSolution, _noDirichlet, b_source);
			boundary = _higherOrderBoundary.Trace(reconstructedElemBoundary);
		}
		else
			boundary = faceSolution.tail(HHO->nBoundaryFaces * HHO->nFaceUnknowns); // keep only the boundary unknowns

		_integralZeroOnBoundary.Enforce(boundary);
		return boundary;
	}

	double L2InnerProdOnBoundary(const Vector& v1, const Vector& v2) override
	{
		return _boundarySpace->L2InnerProd(v1, v2);
	}

	pair<Vector, Vector> ComputeSolution(const Vector& theta) override
	{
		pair<Vector, Vector> p;
		auto& [lambda, solution] = p;

		// Solve problem 1 (f=source, Neum=<theta>)
		lambda = Solve1stDiffProblem(theta);

		// Solve problem 2 (f=<lambda>, Neum=0)
		solution = Solve2ndDiffProblem(lambda);
		
		return p;
	}

	/*DenseMatrix BasisChangeMatrix()
	{
		auto n = _boundarySpace->Dimension();
		DenseMatrix Change(n, n);
		NonZeroCoefficients coeffs;
		for (int i = 0; i < n; i++)
		{
			Vector e_i = Vector::Zero(n);
			e_i[i] = 1;
			_integralZeroOnBoundary.Enforce(e_i);
			Vector lambda = Solve1stDiffProblemWithZeroSource(e_i);
			Change.col(i) = -Solve2ndDiffProblem(lambda, true);
		}
		return Change;
	}*/

	Vector ProblemOperator(const Vector& x) override
	{
		auto delta = Solve1stDiffProblem_Homogeneous(x);
		return -Solve2ndDiffProblem_Homogeneous(delta);
	}

	DenseMatrix Matrix() override
	{
		int n = _boundarySpace->Dimension();
		DenseMatrix A(n, n);
		for (int i = 0; i < n; i++)
		{
			Vector e_i = Vector::Zero(n);
			e_i[i] = 1;
			_integralZeroOnBoundary.Enforce(e_i);
			Vector lambda = Solve1stDiffProblem_Homogeneous(e_i);
			A.col(i) = -Solve2ndDiffProblem_Homogeneous(lambda);
		}
		return A;
	}
};
