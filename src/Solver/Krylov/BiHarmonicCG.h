#pragma once
#include "../IterativeSolver.h"
#include "../../HHO/BiHarmonicMixedForm_HHO.h"
using namespace std;

template <int Dim>
class BiHarmonicCG : public IterativeSolver
{
public:
	BiHarmonicTestCase<Dim>* _testCase;
	Diffusion_HHO<Dim>& _diffPb;
	Solver* _diffSolver = nullptr;

	BiHarmonicCG(BiHarmonicTestCase<Dim>* testCase, Diffusion_HHO<Dim>& diffPb, Solver* diffSolver)
		: _diffPb(diffPb)
	{ 
		_testCase = testCase;
		_diffSolver = diffSolver;
	}

	virtual void Serialize(ostream& os) const override
	{
		os << "Conjugate Gradient for the bi-harmonic problem";
	}

	Vector Solve()
	{
		this->SolvingComputationalWork = 0;

		// 1. Find initial theta verifying the compatibility condition
		//               (source|1) + <theta|1> = 0
		//    We define theta := -(source|1) / |\partial \Omega| * 1
		//    where |\partial \Omega| is the measure of the boundary.

		double integralSource = _diffPb.IntegralOverDomain(_testCase->SourceFunction);
		double boundaryMeasure = _diffPb._mesh->BoundaryMeasure();
		Vector theta = -integralSource / boundaryMeasure * _diffPb.ProjectOnBoundaryDiscreteSpace(Utils::ConstantFunctionOne);

		IterationResult result = CreateFirstIterationResult(Vector::Zero(theta.rows()), theta);

		// Check compatibility condition
		double integralTheta = _diffPb.IntegralOverBoundaryFromFaceCoeffs(theta);
		assert(abs(integralSource + integralTheta) < Utils::Eps);

		//---------------------------------------//
		// 2. 1st problem (f=source, Neum=theta) //
		//---------------------------------------//

		_diffPb.ChangeSourceFunction(_testCase->SourceFunction);
		_diffPb.ChangeNeumannFunction(theta);
		Vector& rhs = _diffPb.SetCondensedRHS();

		// Solve
		_diffPb.SystemSolution = _diffSolver->Solve(rhs);

		// Check that mean value = 0 (on the faces)
		double integralSkeleton = _diffPb.IntegralOverSkeletonFromFaceCoeffs(_diffPb.SystemSolution);
		assert(abs(integralSkeleton) < Utils::Eps);

		// Reconstruct the higher-order polynomial
		_diffPb.ReconstructHigherOrderApproximation();
		Vector& lambda = _diffPb.ReconstructedSolution;

		// Even if the mean value is 0 on the faces, the mean value of the reconstructed polynomial is not necessarily 0,
		// so we enforce it: lambda <- lambda - (lambda|1)/(1|1)*1
		Vector innerProdsWithOne = _diffPb.InnerProdWithReconstructBasis(Utils::ConstantFunctionOne);
		Vector one = _diffPb.SolveReconstructMassMatrix(innerProdsWithOne);
		double normOneSquare = _diffPb._mesh->Measure();
		lambda -= lambda.dot(innerProdsWithOne) / normOneSquare * one;


		//----------------------------------//
		// 3. 2nd problem (f=lamda, Neum=0) //
		//----------------------------------//

		// Verify compatibility condition: (lambda|1) = 0
		double meanValue = _diffPb.MeanValueFromReconstructedCoeffs(lambda);
		assert(abs(meanValue) < Utils::Eps);

		_diffPb.ChangeSourceFunction(lambda);
		_diffPb.ChangeNeumannFunctionToZero();
		rhs = _diffPb.SetCondensedRHS();

		_diffPb.SystemSolution = _diffSolver->Solve(rhs);

		//-----------------------//
		// 4. Conjugate Gradient //
		//-----------------------//

		// Residual on the faces (we want homogeneous Dirichlet BC)
		Vector r = _diffPb.SystemSolution.tail(_diffPb.HHO->nInteriorFaces * _diffPb.HHO->nFaceUnknowns); // keep only the boundary unknowns

		result.SetResidualNorm(r.norm());               result.AddWorkInFlops(Cost::Norm(r));

		// Direction of research
		Vector p = r;

		double r_dot_r = r.dot(r);

		this->IterationCount = 0;

		if (this->PrintIterationResults)
			cout << result << endl;

		while (!StoppingCriteriaReached(result))
		{
			result = IterationResult(result);

			// Solve 1st diffusion problem (f=0, Neum=p)
			Vector& delta = SolveDiffPbSourceZero(p);

			// Solve 2nd diffusion problem (f=delta, Neum=0)
			Vector gamma = SolveDiffPbNeumannZero(delta);

			// Step for theta in the direction of research
			double rho = p.dot(r) / gamma.tail(_diffPb.HHO->nInteriorFaces * _diffPb.HHO->nFaceUnknowns).dot(p);
			// Move theta in the direction of research
			theta += rho * p;

			double r_dot_r_old = r_dot_r; // save the dot product before overwriting r and p

			// Update residual
			r -= rho * gamma.tail(_diffPb.HHO->nInteriorFaces * _diffPb.HHO->nFaceUnknowns);

			r_dot_r = r.dot(r);

			// Step for the direction of research
			double q = r_dot_r / r_dot_r_old;
			// Update the direction of research
			p = r + q * p;

			this->IterationCount++;

			result.SetX(theta);
			result.SetResidualNorm(r.norm());                   result.AddWorkInFlops(Cost::Norm(r));

			if (this->PrintIterationResults)
				cout << result << endl;
		}

		if (this->PrintIterationResults)
			cout << endl;

		this->SolvingComputationalWork = result.SolvingComputationalWork();

		return theta;
	}

private:
	// Solve 1st diffusion problem (f=0, Neum=<neumann>)
	Vector& SolveDiffPbSourceZero(const Vector& neumann)
	{
		_diffPb.ChangeSourceFunctionToZero();
		_diffPb.ChangeNeumannFunction(neumann);
		Vector& rhs = _diffPb.SetCondensedRHS();
		_diffPb.SystemSolution = _diffSolver->Solve(rhs);
		_diffPb.ReconstructHigherOrderApproximation();
		return _diffPb.ReconstructedSolution;
	}

	// Solve 2nd diffusion problem (f=<source>, Neum=0)
	Vector SolveDiffPbNeumannZero(const Vector& source)
	{
		_diffPb.ChangeSourceFunction(source);
		_diffPb.ChangeNeumannFunctionToZero();
		Vector& rhs = _diffPb.SetCondensedRHS();
		return _diffSolver->Solve(rhs);
	}
};