#pragma once
#include "BiHarmonicMixedForm_FEM.h"
#include "../../TestCases/BiHarmonic/BiHarmonicTestCase.h"
#include "../../TestCases/Diffusion/VirtualDiffusionTestCase.h"
using namespace std;

template<int Dim>
class BiHarmonicMixedFormGlowinski_FEM : public BiHarmonicMixedForm_FEM<Dim>
{
private:
	Mesh<Dim>* _mesh;
	BiHarmonicTestCase<Dim>* _testCase;

	DiffusionField<Dim> _diffField;
	VirtualDiffusionTestCase<Dim> _diffPbTestCase;
	Diffusion_FEM<Dim> _diffPb;

	SparseMatrix _massMatrix;
public:

	BiHarmonicMixedFormGlowinski_FEM(Mesh<Dim>* mesh, BiHarmonicTestCase<Dim>* testCase, FunctionalBasis<Dim>* basis)
	{
		_mesh = mesh;
		_testCase = testCase;
		_diffField = DiffusionField<Dim>();
		mesh->SetDiffusionField(&_diffField);
		_diffPbTestCase = VirtualDiffusionTestCase<Dim>(testCase->SourceFunction, _diffField);
		_diffPbTestCase.BC = BoundaryConditions::HomogeneousDirichletEverywhere();
		_diffPb = Diffusion_FEM<Dim>(mesh, &_diffPbTestCase, basis);
	}

	Diffusion_FEM<Dim>& DiffPb() override
	{
		return _diffPb;
	}

	void Setup() override
	{
		ActionsArguments diffActions;
		diffActions.AssembleRightHandSide = false;
		diffActions.LogAssembly = false;
		_diffPb.Assemble(diffActions);

		_massMatrix = _diffPb.MassMatrix();
	}

	Vector FindCompatibleTheta() override
	{
		return Vector::Zero(_mesh->BoundaryVertices.size());
	}

	// Solve problem 1 (f=source, Dirich=<dirichlet>)
	Vector Solve1stDiffProblemWithFSource(const Vector& dirichlet) override
	{
		// Define problem
		Vector b_source = _diffPb.AssembleSourceTerm(_testCase->SourceFunction);
		Vector rhs = _diffPb.RHS(b_source, dirichlet);

		// Solve
		Vector solution = this->_diffSolver->Solve(rhs);
		this->CheckDiffSolverConvergence();

		return _diffPb.BuildCompleteSolution(solution, dirichlet);
		//return solution;
	}

	// Solve problem 1 (f=0, Dirich=<dirichlet>)
	Vector Solve1stDiffProblemWithZeroSource(const Vector& dirichlet) override
	{
		// Define problem
		Vector b_source = Vector::Zero(_mesh->InteriorVertices.size());
		Vector rhs = _diffPb.RHS(b_source, dirichlet);

		// Solve
		Vector solution = this->_diffSolver->Solve(rhs);
		this->CheckDiffSolverConvergence();

		return _diffPb.BuildCompleteSolution(solution, dirichlet);
		//return solution;
	}

	// Solve problem 2 (f=<source>, Dirich=0)
	Vector Solve2ndDiffProblem(const Vector& sourceNodeValues, bool returnBoundaryNormalDerivative = false) override
	{
		// Define problem
		Vector dirichlet = Vector::Zero(_mesh->BoundaryVertices.size());
		Vector b_source = _diffPb.AssembleSourceTerm(sourceNodeValues, _massMatrix);
		Vector rhs = _diffPb.RHS(b_source, dirichlet);

		// Solve
		Vector solution = this->_diffSolver->Solve(rhs);
		this->CheckDiffSolverConvergence();

		if (returnBoundaryNormalDerivative)
		{
			return _diffPb.A_id.transpose() * solution - _massMatrix.bottomRightCorner(_mesh->BoundaryVertices.size(), _mesh->BoundaryVertices.size()) * sourceNodeValues.tail(_mesh->BoundaryVertices.size());
		}
		else
			return _diffPb.BuildCompleteSolution(solution, dirichlet);
	}

	double L2InnerProdOnBoundary(const Vector& v1, const Vector& v2) override
	{
		//assert(false);
		return v1.dot(v2);
	}

	pair<Vector, Vector> ComputeSolution(const Vector& theta) override
	{
		pair<Vector, Vector> p;
		auto& [lambda, solution] = p;

		// Solve problem 1 (f=source, Dirich=<theta>)
		lambda = Solve1stDiffProblemWithFSource(theta);

		// Solve problem 2 (f=<lambda>, Dirich=0)
		solution = Solve2ndDiffProblem(lambda);

		return p;
	}
};
