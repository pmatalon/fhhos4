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
	SparseMatrix _HPTranspose;

	SparseMatrix _CellStiff;
	SparseMatrix _CellMass;
	SparseMatrix _DomainSolveCellUknTranspose;
	SparseMatrix _SolveCellUknTranspose;

	SparseMatrix _A_T_T;
	SparseMatrix _A_T_T_Stab;
	
	SparseMatrix _Trace;
	SparseMatrix _normalDerivativeMatrix;

	SparseMatrix _Proj_HO_LO_Bndry;

	EigenSparseCholesky _PT_M_P_solver;

	HHOParameters<Dim>* HHO;

	Vector _b_fSource;
	Vector _zeroDirichlet;

	int _option = 0;
public:
	bool UseIntegrationByParts = true;

	HHOBoundarySpace<Dim>* ThetaSpace = nullptr;

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
		_diffPb.Assemble(diffActions, true);

		ThetaSpace = &_diffPb.BoundarySpace;

		_option = Utils::ProgramArgs.Actions.Option;

		if (UseIntegrationByParts)
		{
			_ReconstructStiff = _diffPb.ReconstructStiffnessMatrixOnBoundaryElements();
			_ReconstructMass = _diffPb.ReconstructMassMatrixOnBoundaryElements();
			if (_option == 4 || _option == 9)
				_HPTranspose = _diffPb.H0PTransposeOnBoundary();
			else
				_HPTranspose = _diffPb.HPTransposeOnBoundary();
			_DomainSolveCellUknTranspose = _diffPb.DomainSolveCellUknTransposeOnBoundary();
			_SolveCellUknTranspose = _diffPb.SolveCellUknTransposeOnBoundary();

			_A_T_T = _diffPb.A_T_T_Matrix();
			_A_T_T_Stab = _diffPb.A_T_T_Stab_Matrix();

			if (_option == 3)
				_Trace = _diffPb.LowerOrderTraceMatrixOnBoundary();

			if (_option == 1)
			{
				_higherOrderBoundary = HigherOrderBoundary<Dim>(&_diffPb);
				_higherOrderBoundary.Setup();

				SparseMatrix& hoTrace = _higherOrderBoundary.TraceMatrix();
				SparseMatrix M = _HPTranspose * hoTrace.transpose() * _higherOrderBoundary.BoundaryFaceMassMatrix() * hoTrace * _HPTranspose.transpose();
				_PT_M_P_solver.Setup(M);
			}
			if (_option == 5 || _option == 11 || _option == 13)
			{
				_CellStiff = _diffPb.CellStiffnessMatrixOnBoundaryElements();
				_CellMass = _diffPb.CellMassMatrixOnBoundaryElements();
			}
			if (_option == 7)
			{
				_higherOrderBoundary = HigherOrderBoundary<Dim>(&_diffPb);
				_higherOrderBoundary.Setup();

				ThetaSpace = &_higherOrderBoundary.BoundarySpace;

				_Proj_HO_LO_Bndry = LowerOrderProjectionOnBoundary(_higherOrderBoundary);
			}
			if (_option == 8)
			{
				_higherOrderBoundary = HigherOrderBoundary<Dim>(&_diffPb);
				_higherOrderBoundary.Setup();

				ThetaSpace = &_higherOrderBoundary.BoundarySpace;

				_Trace = _higherOrderBoundary.TraceMatrix();
			}
		}
		else
			_normalDerivativeMatrix = _diffPb.NormalDerivativeMatrix();

		_b_fSource = _diffPb.AssembleSourceTerm(_testCase->SourceFunction);
		_zeroDirichlet = Vector::Zero(HHO->nDirichletCoeffs);
	}

	Vector FindCompatibleTheta() override
	{
		if (UseIntegrationByParts)
		{
			if (_option == 2)
				return Vector::Zero(_mesh->NBoundaryElements() * HHO->nReconstructUnknowns);
			else
				return ThetaSpace->ZeroVector();
		}
		else
			return ThetaSpace->ZeroVector();
	}

private:
	Vector BuildDirichlet(const Vector& dirichletArg)
	{
		if (UseIntegrationByParts)
		{
			if (_option == 0)
				return dirichletArg;
			else if (_option == 1)
				return _diffPb.ReconstructAndAssembleDirichletTerm(dirichletArg);
			else if (_option == 2)
				return _diffPb.AssembleDirichletTermFromReconstructedBoundaryElem(dirichletArg);
			else if (_option == 6 || _option == 4)
				return ThetaSpace->SolveMassMatrix(dirichletArg);
			else if (_option == 7)
				return _higherOrderBoundary.AssembleDirichletTerm(dirichletArg);
			else if (_option == 8)
				return _higherOrderBoundary.AssembleDirichletTerm(dirichletArg);
			else
				return dirichletArg;
		}
		else
			return dirichletArg;
	}

public:
	// Solve problem 1 (f=source, Dirich=<dirichlet>)
	Vector Solve1stDiffProblemWithFSource(const Vector& dirichletArg) override
	{
		// Define problem
		Vector dirichlet = BuildDirichlet(dirichletArg);
		Vector b_T   = _diffPb.ComputeB_T(_b_fSource, dirichlet);
		Vector b_ndF = _diffPb.ComputeB_ndF_noNeumann(dirichlet);
		Vector rhs = _diffPb.CondensedRHS(b_T, b_ndF);

		// Solve
		Vector faceSolution = this->_diffSolver->Solve(rhs);
		this->CheckDiffSolverConvergence();

		if (_option == 5 || _option == 11 || _option == 13 || _option == 14)
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
		Vector dirichlet = BuildDirichlet(dirichletArg);
		Vector b_T = _diffPb.ComputeB_T_zeroSource(dirichlet);
		Vector b_ndF = _diffPb.ComputeB_ndF_noNeumann(dirichlet);
		Vector rhs = _diffPb.CondensedRHS(b_T, b_ndF);

		// Solve
		Vector faceSolution = this->_diffSolver->Solve(rhs);
		this->CheckDiffSolverConvergence();

		if (_option == 5 || _option == 11 || _option == 13 || _option == 14)
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
				if (_option == 0)
					return ThetaSpace->SolveMassMatrix(_HPTranspose * (_ReconstructStiff * reconstructedElemBoundary - _ReconstructMass * sourceElemBoundary));
				else if (_option == 1)
					return _PT_M_P_solver.Solve(_HPTranspose * (_ReconstructStiff * reconstructedElemBoundary - _ReconstructMass * sourceElemBoundary));
				else if (_option == 2)
					return _ReconstructStiff * reconstructedElemBoundary - _ReconstructMass * sourceElemBoundary;
				else if (_option == 3)
					return ThetaSpace->SolveMassMatrix(_Trace * (_ReconstructStiff * reconstructedElemBoundary - _ReconstructMass * sourceElemBoundary));
				else if (_option == 4)
				{
					// Second best
					Vector cellSolution = _diffPb.SolveCellUnknowns(faceSolution, b_source);
					Vector normalDerivative = _diffPb.A_T_dF.transpose() * cellSolution + _diffPb.A_ndF_dF.transpose() * faceSolution;
					normalDerivative -= _HPTranspose * _ReconstructMass * sourceElemBoundary;
					return normalDerivative;
				}
				else if (_option == 5)
				{
					Vector cellSolution = _diffPb.SolveCellUnknownsOnBoundaryOnly(faceSolution, b_source);
					return ThetaSpace->SolveMassMatrix(_SolveCellUknTranspose * (_CellStiff * cellSolution - _CellMass * sourceElemBoundary));
				}
				else if (_option == 6)
					return _HPTranspose * (_ReconstructStiff * reconstructedElemBoundary - _ReconstructMass * sourceElemBoundary);
				else if (_option == 7)
					return ThetaSpace->SolveMassMatrix(_Proj_HO_LO_Bndry.transpose() * _HPTranspose * (_ReconstructStiff * reconstructedElemBoundary - _ReconstructMass * sourceElemBoundary));
				else if (_option == 8)
					return ThetaSpace->SolveMassMatrix(_Trace * (_ReconstructStiff * reconstructedElemBoundary - _ReconstructMass * sourceElemBoundary));
				else if (_option == 9) // close to 4
				{
					//Vector cellSolution = _diffPb.SolveCellUnknownsOnBoundaryOnly(faceSolution, b_source);
					//_diffPb.StabilizationMatrixOnBoundaryFaces();

					Vector cellSolution = _diffPb.SolveCellUnknowns(faceSolution, b_source);
					Vector normalDerivative = _diffPb.A_T_dF.transpose()      * cellSolution + _diffPb.A_ndF_dF.transpose()      * faceSolution;
					normalDerivative       -= _diffPb.A_Stab_T_dF.transpose() * cellSolution + _diffPb.A_Stab_ndF_dF.transpose() * faceSolution;
					normalDerivative -= _HPTranspose * _ReconstructMass * sourceElemBoundary;
					return normalDerivative;
				}
				else if (_option == 10) // close to 0
				{
					Vector cellSolution = _diffPb.SolveCellUnknowns(faceSolution, b_source);
					return ThetaSpace->SolveMassMatrix(_HPTranspose * (_ReconstructStiff * reconstructedElemBoundary - _ReconstructMass * sourceElemBoundary) - _diffPb.A_Stab_T_dF.transpose() * cellSolution - _diffPb.A_Stab_ndF_dF.transpose() * faceSolution);
				}
				else if (_option == 11) // close to 5
				{
					Vector cellSolution = _diffPb.SolveCellUnknownsOnBoundaryOnly(faceSolution, b_source);
					Vector normalDerivative = _SolveCellUknTranspose * (_CellStiff * cellSolution - _CellMass * sourceElemBoundary);

					cellSolution = _diffPb.SolveCellUnknowns(faceSolution, b_source);
					normalDerivative -= _diffPb.A_Stab_T_dF.transpose() * cellSolution + _diffPb.A_Stab_ndF_dF.transpose() * faceSolution;
					return ThetaSpace->SolveMassMatrix(normalDerivative);
				}
				else if (_option == 12) // close to 9
				{
					Vector cellSolution = _diffPb.SolveCellUnknowns(faceSolution, b_source);

					Vector normalDerivative = _diffPb.A_T_dF.transpose() * cellSolution + _diffPb.A_ndF_dF.transpose() * faceSolution;
					normalDerivative += _DomainSolveCellUknTranspose * (_A_T_T * cellSolution + _diffPb.A_T_ndF * faceSolution);

					normalDerivative       -= _diffPb.A_Stab_T_dF.transpose() * cellSolution + _diffPb.A_Stab_ndF_dF.transpose() * faceSolution;
					normalDerivative -= _DomainSolveCellUknTranspose * (_A_T_T_Stab * cellSolution + _diffPb.A_Stab_T_ndF * faceSolution);
					
					normalDerivative -= _HPTranspose * _ReconstructMass * sourceElemBoundary;
					return ThetaSpace->SolveMassMatrix(normalDerivative);
				}
				else if (_option == 13) // close to 12
				{
					// Best (with -kc 1)

					Vector cellSolution = _diffPb.SolveCellUnknowns(faceSolution, b_source);

					Vector normalDerivative = _diffPb.A_T_dF.transpose() * cellSolution + _diffPb.A_ndF_dF.transpose() * faceSolution;
					normalDerivative += _DomainSolveCellUknTranspose * (_A_T_T * cellSolution + _diffPb.A_T_ndF * faceSolution);

					normalDerivative -= _SolveCellUknTranspose * _CellMass * sourceElemBoundary;

					//normalDerivative += _diffPb.A_Stab_T_dF.transpose() * cellSolution + _diffPb.A_Stab_ndF_dF.transpose() * faceSolution;
					//normalDerivative += _DomainSolveCellUknTranspose * (_A_T_T_Stab * cellSolution + _diffPb.A_Stab_T_ndF * faceSolution);
					return ThetaSpace->SolveMassMatrix(normalDerivative);
				}
				else if (_option == 14) // close to 13
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

private:
	SparseMatrix LowerOrderProjectionOnBoundary(HigherOrderBoundary<Dim>& hoBoundary)
	{
		int hoUnknowns = hoBoundary.HHO->nFaceUnknowns;
		int loUnknowns = _diffPb.HHO->nFaceUnknowns;

		FaceParallelLoop<Dim> parallelLoop(_mesh->BoundaryFaces);
		parallelLoop.ReserveChunkCoeffsSize(hoUnknowns * loUnknowns);

		parallelLoop.Execute([this, &hoBoundary](Face<Dim>* f, ParallelChunk<CoeffsChunk>* chunk)
			{
				Diff_HHOFace<Dim>* hoFace = hoBoundary.HHOFace(f);
				int hoUnknowns = hoBoundary.HHO->nFaceUnknowns;
				BigNumber j = f->Number - HHO->nInteriorFaces;

				Diff_HHOFace<Dim>* loFace = _diffPb.HHOFace(f);
				int loUnknowns = _diffPb.HHO->nFaceUnknowns;
				BigNumber i = f->Number - HHO->nInteriorFaces;

				DenseMatrix proj = loFace->ProjectOnBasis(hoFace->Basis);

				chunk->Results.Coeffs.Add(i * loUnknowns, j * hoUnknowns, proj);
				
			});

		SparseMatrix mat(HHO->nBoundaryFaces * loUnknowns, HHO->nBoundaryFaces * hoUnknowns);
		parallelLoop.Fill(mat);
		return mat;
	}
};
