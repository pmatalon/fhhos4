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

	HigherOrderBoundary<Dim> _higherOrderBoundary;

	SparseMatrix _ReconstructStiff;
	SparseMatrix _ReconstructMass;
	SparseMatrix _PTranspose;

	SparseMatrix _CellStiff;
	SparseMatrix _CellMass;
	SparseMatrix _SolveCellUknTranspose;
	
	SparseMatrix _Trace;
	SparseMatrix _normalDerivativeMatrix;

	EigenSparseCholesky _PT_M_P_solver;

	HHOParameters<Dim>* HHO;

	Vector _b_fSource;
	Vector _zeroDirichlet;
public:
	bool UseIntegrationByParts = true;

	BiHarmonicMixedFormGlowinski_HHO(Mesh<Dim>* mesh, BiHarmonicTestCase<Dim>* testCase, HHOParameters<Dim>* hho, bool useIntegrationByParts, bool saveMatrixBlocks)
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
		UseIntegrationByParts = useIntegrationByParts;
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

		if (UseIntegrationByParts)
		{
			_ReconstructStiff = _diffPb.ReconstructStiffnessMatrixOnBoundaryElements();
			_ReconstructMass = _diffPb.ReconstructMassMatrixOnBoundaryElements();
			_PTranspose = _diffPb.PTransposeOnBoundary();
			
			if (Utils::ProgramArgs.Actions.Option == 3)
				_Trace = _diffPb.LowerOrderTraceMatrixOnBoundary();

			if (Utils::ProgramArgs.Actions.Option == 1)
			{
				_higherOrderBoundary = HigherOrderBoundary<Dim>(&_diffPb);
				_higherOrderBoundary.Setup();

				SparseMatrix& hoTrace = _higherOrderBoundary.TraceMatrix();
				SparseMatrix M = _PTranspose * hoTrace.transpose() * _higherOrderBoundary.BoundaryFaceMassMatrix() * hoTrace * _PTranspose.transpose();
				_PT_M_P_solver.Setup(M);
			}
			if (Utils::ProgramArgs.Actions.Option == 5)
			{
				_CellStiff = _diffPb.CellStiffnessMatrixOnBoundaryElements();
				_CellMass = _diffPb.CellMassMatrixOnBoundaryElements();
				_SolveCellUknTranspose = _diffPb.SolveCellUknTransposeOnBoundary();
			}
		}
		else
			_normalDerivativeMatrix = _diffPb.NormalDerivativeMatrix();

		_b_fSource = _diffPb.AssembleSourceTerm(_testCase->SourceFunction);
		_zeroDirichlet = Vector::Zero(HHO->nDirichletCoeffs);
	}

	Vector FindCompatibleTheta() override
	{
		if (Utils::ProgramArgs.Actions.Option == 2)
			return Vector::Zero(_mesh->NBoundaryElements() * HHO->nReconstructUnknowns);
		else
			return _zeroDirichlet;
	}

	// Solve problem 1 (f=source, Dirich=<dirichlet>)
	Vector Solve1stDiffProblemWithFSource(const Vector& dirichletArg) override
	{
		// Define problem
		Vector dirichlet;
		if (Utils::ProgramArgs.Actions.Option == 0)
			dirichlet = dirichletArg;
		else if (Utils::ProgramArgs.Actions.Option == 1)
			dirichlet = _diffPb.ReconstructAndAssembleDirichletTerm(dirichletArg);
		else if (Utils::ProgramArgs.Actions.Option == 2)
			dirichlet = _diffPb.AssembleDirichletTermFromReconstructedBoundaryElem(dirichletArg);
		else if (Utils::ProgramArgs.Actions.Option == 6)
			dirichlet = _diffPb.SolveFaceMassMatrixOnBoundary(dirichletArg);
		else
			dirichlet = dirichletArg;

		Vector b_T   = _diffPb.ComputeB_T(_b_fSource, dirichlet);
		Vector b_ndF = _diffPb.ComputeB_ndF_noNeumann(dirichlet);
		Vector rhs = _diffPb.CondensedRHS(b_T, b_ndF);

		// Solve
		Vector faceSolution = this->_diffSolver->Solve(rhs);
		this->CheckDiffSolverConvergence();

		if (Utils::ProgramArgs.Actions.Option == 5)
			return _diffPb.SolveCellUnknowns(faceSolution, b_T);
		else
		{
			// Reconstruct the higher-order polynomial
			Vector lambda = _diffPb.ReconstructHigherOrderApproximationFromFaceCoeffs(faceSolution, dirichlet, b_T);

			//ExportModule out(Utils::ProgramArgs.OutputDirectory, "", Utils::ProgramArgs.Actions.Export.ValueSeparator);
			//_diffPb.ExportReconstructedVectorToGMSH(lambda, out, "lambda_f");
			
			return lambda;
		}
	}

	// Solve problem 1 (f=0, Dirich=<dirichlet>)
	Vector Solve1stDiffProblemWithZeroSource(const Vector& dirichletArg) override
	{
		// Define problem
		Vector dirichlet;
		if (Utils::ProgramArgs.Actions.Option == 0)
			dirichlet = dirichletArg;
		else if (Utils::ProgramArgs.Actions.Option == 1)
			dirichlet = _diffPb.ReconstructAndAssembleDirichletTerm(dirichletArg);
		else if (Utils::ProgramArgs.Actions.Option == 2)
			dirichlet = _diffPb.AssembleDirichletTermFromReconstructedBoundaryElem(dirichletArg);
		else if (Utils::ProgramArgs.Actions.Option == 6 || Utils::ProgramArgs.Actions.Option == 4)
			dirichlet = _diffPb.SolveFaceMassMatrixOnBoundary(dirichletArg);
		else
			dirichlet = dirichletArg;

		Vector b_T = _diffPb.ComputeB_T_zeroSource(dirichlet);
		Vector b_ndF = _diffPb.ComputeB_ndF_noNeumann(dirichlet);
		Vector rhs = _diffPb.CondensedRHS(b_T, b_ndF);

		// Solve
		Vector faceSolution = this->_diffSolver->Solve(rhs);
		this->CheckDiffSolverConvergence();

		if (Utils::ProgramArgs.Actions.Option == 5)
			return _diffPb.SolveCellUnknowns(faceSolution, b_T);
		else
		{
			// Reconstruct the higher-order polynomial
			Vector lambda = _diffPb.ReconstructHigherOrderApproximationFromFaceCoeffs(faceSolution, dirichlet, b_T);

			//ExportModule out(Utils::ProgramArgs.OutputDirectory, "", Utils::ProgramArgs.Actions.Export.ValueSeparator);
			//_diffPb.ExportReconstructedVectorToGMSH(lambda, out, "lambda_0");

			return lambda;
		}
	}

	// Solve problem 2 (f=<source>, Dirich=0)
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
			if (UseIntegrationByParts)
			{
				//Vector reconstructedSolution = _diffPb.ReconstructHigherOrderApproximationFromFaceCoeffs(faceSolution, _zeroDirichlet, b_source);
				//return _diffPb.SolveFaceMassMatrixOnBoundary(_PTranspose * (_Stiff * reconstructedSolution - _Mass * source));

				Vector reconstructedElemBoundary = _diffPb.ReconstructHigherOrderOnBoundaryOnly(faceSolution, _zeroDirichlet, b_source);
				Vector sourceElemBoundary = _diffPb.ExtractElemBoundary(source);
				if (Utils::ProgramArgs.Actions.Option == 0)
					return _diffPb.SolveFaceMassMatrixOnBoundary(_PTranspose * (_ReconstructStiff * reconstructedElemBoundary - _ReconstructMass * sourceElemBoundary));
				else if (Utils::ProgramArgs.Actions.Option == 1)
					return _PT_M_P_solver.Solve(_PTranspose * (_ReconstructStiff * reconstructedElemBoundary - _ReconstructMass * sourceElemBoundary));
				else if (Utils::ProgramArgs.Actions.Option == 2)
					return _ReconstructStiff * reconstructedElemBoundary - _ReconstructMass * sourceElemBoundary;
				else if (Utils::ProgramArgs.Actions.Option == 3)
					return _diffPb.SolveFaceMassMatrixOnBoundary(_Trace * (_ReconstructStiff * reconstructedElemBoundary - _ReconstructMass * sourceElemBoundary));
				else if (Utils::ProgramArgs.Actions.Option == 4)
				{
					Vector cellSolution = _diffPb.SolveCellUnknowns(faceSolution, b_source);
					Vector normalDerivative = _diffPb.A_T_dF.transpose() * cellSolution + _diffPb.A_ndF_dF.transpose() * faceSolution;
					normalDerivative -= _PTranspose * _ReconstructMass * sourceElemBoundary;
					return normalDerivative;
				}
				else if (Utils::ProgramArgs.Actions.Option == 5)
				{
					Vector cellSolution = _diffPb.SolveCellUnknownsOnBoundaryOnly(faceSolution, b_source);
					return _diffPb.SolveFaceMassMatrixOnBoundary(_SolveCellUknTranspose * (_CellStiff * cellSolution - _CellMass * sourceElemBoundary));
				}
				else if (Utils::ProgramArgs.Actions.Option == 6)
					return _PTranspose * (_ReconstructStiff * reconstructedElemBoundary - _ReconstructMass * sourceElemBoundary);
				else
				{
					Utils::FatalError("Unmanaged -opt " + to_string(Utils::ProgramArgs.Actions.Option));
					return Vector();
				}
			}
			else
			{
				Vector reconstructedElemBoundary = _diffPb.ReconstructHigherOrderOnBoundaryOnly(faceSolution, _zeroDirichlet, b_source);
				return _normalDerivativeMatrix * reconstructedElemBoundary;
			}
		}
		else
			return _diffPb.ReconstructHigherOrderApproximationFromFaceCoeffs(faceSolution, _zeroDirichlet, b_source);
	}

	double L2InnerProdOnBoundary(const Vector& v1, const Vector& v2) override
	{
		return _diffPb.BoundarySpace.L2InnerProd(v1, v2);
	}

	Vector ComputeSolution(const Vector& theta) override
	{
		// Solve problem 1 (f=source, Dirich=<theta>)
		Vector lambda = Solve1stDiffProblemWithFSource(theta);

		// Solve problem 2 (f=<lambda>, Dirich=0)
		return Solve2ndDiffProblem(lambda);
	}
};
