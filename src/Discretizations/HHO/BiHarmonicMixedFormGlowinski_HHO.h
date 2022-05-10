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


	SparseMatrix _CellMass;
	SparseMatrix _DomainSolveCellUknTranspose;
	SparseMatrix _SolveCellUknTranspose;
	SparseMatrix _A_T_T;

	HHOParameters<Dim>* HHO;

	Vector _b_fSource;
	Vector _zeroDirichlet;

	int _option = 0;
public:

	HHOBoundarySpace<Dim>* ThetaSpace = nullptr;

	BiHarmonicMixedFormGlowinski_HHO(Mesh<Dim>* mesh, BiHarmonicTestCase<Dim>* testCase, HHOParameters<Dim>* hho, bool saveMatrixBlocks)
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

		ThetaSpace = &_diffPb.BoundarySpace;

		_option = Utils::ProgramArgs.Actions.Option1;

		_DomainSolveCellUknTranspose = _diffPb.DomainSolveCellUknTransposeOnBoundary();
		_SolveCellUknTranspose = _diffPb.SolveCellUknTransposeOnBoundary();
		_A_T_T = _diffPb.A_T_T_Matrix();
		_CellMass = _diffPb.CellMassMatrixOnBoundaryElements();

		_b_fSource = _diffPb.AssembleSourceTerm(_testCase->SourceFunction);
		_zeroDirichlet = Vector::Zero(HHO->nDirichletCoeffs);
	}

	Vector FindCompatibleTheta() override
	{
		return ThetaSpace->ZeroVector();
	}

private:
	Vector BuildDirichlet(const Vector& dirichletArg)
	{
		return dirichletArg;
	}

public:
	// Solve problem 1 (source=f, Dirich=<dirichlet>)
	Vector Solve1stDiffProblemWithFSource(const Vector& dirichletArg) override
	{
		Vector dummy;
		return Solve1stDiffProblemWithFSource(dirichletArg, false, dummy);
	}

private:
	// Solve problem 1 (source=f, Dirich=<dirichlet>)
	Vector Solve1stDiffProblemWithFSource(const Vector& dirichletArg, bool reconstruct, Vector& reconstruction)
	{
		// Define problem
		Vector dirichlet = BuildDirichlet(dirichletArg);
		Vector b_T   = _diffPb.ComputeB_T(_b_fSource, dirichlet);
		Vector b_ndF = _diffPb.ComputeB_ndF_noNeumann(dirichlet);
		Vector rhs = _diffPb.CondensedRHS(b_T, b_ndF);

		// Solve
		Vector faceSolution = this->_diffSolver->Solve(rhs);
		this->CheckDiffSolverConvergence();

		Vector cellSolution = _diffPb.SolveCellUnknowns(faceSolution, b_T);

		if (reconstruct)
		{
			// Reconstruct the higher-order polynomial
			reconstruction = _diffPb.ReconstructHigherOrderApproximationFromFaceCoeffs(faceSolution, dirichlet, b_T);

			//ExportModule out(Utils::ProgramArgs.OutputDirectory, "", Utils::ProgramArgs.Actions.Export.ValueSeparator);
			//_diffPb.ExportReconstructedVectorToGMSH(reconstruction, out, "lambda_f");
		}

		return cellSolution;
	}

	// Solve problem 1 (source=0, Dirich=<dirichlet>)
	Vector Solve1stDiffProblemWithZeroSource(const Vector& dirichletArg) override
	{
		// Define problem
		Vector dirichlet = BuildDirichlet(dirichletArg);
		Vector b_T = _diffPb.ComputeB_T_zeroSource(dirichlet);
		Vector b_ndF = _diffPb.ComputeB_ndF_noNeumann(dirichlet);
		Vector rhs = _diffPb.CondensedRHS(b_T, b_ndF);

		// Solve
		Vector faceSolution = this->_diffSolver->Solve(rhs);
		this->CheckDiffSolverConvergence();

		Vector cellSolution = _diffPb.SolveCellUnknowns(faceSolution, b_T);

		/*Vector reconstruction = _diffPb.ReconstructHigherOrderApproximationFromFaceCoeffs(faceSolution, dirichlet, b_T);
		ExportModule out(Utils::ProgramArgs.OutputDirectory, "", Utils::ProgramArgs.Actions.Export.ValueSeparator);
		_diffPb.ExportReconstructedVectorToGMSH(reconstruction, out, "lambda_0");*/

		return cellSolution;
	}

	// Solve problem 2 (source=<source>, Dirich=0)
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
			Vector sourceElemBoundary = _diffPb.ExtractElemBoundary(source);
			if (_option == 0)
			{
				Vector cellSolution = _diffPb.SolveCellUnknowns(faceSolution, b_source);

				Vector normalDerivative = _diffPb.A_T_dF.transpose() * cellSolution + _diffPb.A_ndF_dF.transpose() * faceSolution;
				normalDerivative += _DomainSolveCellUknTranspose * (_A_T_T * cellSolution + _diffPb.A_T_ndF * faceSolution);

				normalDerivative -= _SolveCellUknTranspose * _CellMass * sourceElemBoundary;

				return ThetaSpace->SolveMassMatrix(normalDerivative);
			}
			else if (_option == 1) // close to 0
			{
				//Vector cellSolution = _diffPb.SolveCellUnknowns(faceSolution, b_source);

				//Vector normalDerivative = _diffPb.A_T_dF.transpose() * cellSolution + _diffPb.A_ndF_dF.transpose() * faceSolution;
				//normalDerivative += _DomainSolveCellUknTranspose * (_A_T_T * cellSolution + _diffPb.A_T_ndF * faceSolution);
				Vector normalDerivative = (_diffPb.A * faceSolution).tail(HHO->nBoundaryFaces * HHO->nFaceUnknowns);

				normalDerivative -= _SolveCellUknTranspose * _CellMass * sourceElemBoundary;
				return ThetaSpace->SolveMassMatrix(normalDerivative);
			}
			else
			{
				Utils::FatalError("Unmanaged -opt " + to_string(_option));
				return Vector();
			}
		}
		else
			return _diffPb.ReconstructHigherOrderApproximationFromFaceCoeffs(faceSolution, _zeroDirichlet, b_source);
	}

	double L2InnerProdOnBoundary(const Vector& v1, const Vector& v2) override
	{
		return _diffPb.BoundarySpace.L2InnerProd(v1, v2);
	}

	pair<Vector, Vector> ComputeSolution(const Vector& theta) override
	{
		pair<Vector, Vector> p;
		auto& [lambda, solution] = p;

		// Solve problem 1 (source=f, Dirich=<theta>)
		Vector lambdaCells = Solve1stDiffProblemWithFSource(theta, true, lambda);

		// Solve problem 2 (source=<lambda>, Dirich=0)
		solution = Solve2ndDiffProblem(lambdaCells);

		return p;
	}
};
