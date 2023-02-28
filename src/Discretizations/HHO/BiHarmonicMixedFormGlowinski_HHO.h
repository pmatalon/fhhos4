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

	DiffusionField<Dim> _diffField;
	VirtualDiffusionTestCase<Dim> _diffPbTestCase;
	Diffusion_HHO<Dim> _diffPb;

	SparseMatrix _Theta_T_bF_transpose;
	SparseMatrix _S_iF_bF_transpose;
	SparseMatrix _S_bF_bF;
	SparseMatrix _Stab_iF_T;
	SparseMatrix _Stab_iF_F;
	SparseMatrix _Stab_bF_T;
	SparseMatrix _Stab_bF_F;

	HHOParameters<Dim>* HHO;

	Vector _g_D;
	Vector _b_fSource;
public:
	bool print = false;

	BiHarmonicMixedFormGlowinski_HHO(Mesh<Dim>* mesh, BiHarmonicTestCase<Dim>* testCase, HHOParameters<Dim>* hho, bool saveMatrixBlocks)
	{
		_mesh = mesh;
		_testCase = testCase;
		HHO = hho;

		map<string, Tensor<Dim>> tensors;
		for (PhysicalGroup<Dim>* phyPart : mesh->PhysicalParts)
		{
			tensors.insert({ phyPart->Name, Tensor<Dim>() });
		}
		_diffField = DiffusionField<Dim>(tensors);
		mesh->SetDiffusionField(&_diffField);

		_diffPbTestCase = VirtualDiffusionTestCase<Dim>(testCase->SourceFunction, _diffField);
		_diffPbTestCase.BC = BoundaryConditions::HomogeneousDirichletEverywhere();
		_diffPb = Diffusion_HHO<Dim>(mesh, &_diffPbTestCase, HHO, true, true);
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

		SparseMatrix Theta_T_F_transpose = _diffPb.Theta_T_F_transpose();
		SparseMatrix Theta_T_iF_transpose = Theta_T_F_transpose.topRows(HHO->nInteriorFaces * HHO->nFaceUnknowns);
		_Theta_T_bF_transpose = Theta_T_F_transpose.bottomRows(HHO->nBoundaryFaces * HHO->nFaceUnknowns);
		auto nBoundaryElemUnknowns = _mesh->NBoundaryElements() * HHO->nCellUnknowns;
		
		_S_iF_bF_transpose = _Theta_T_bF_transpose * _diffPb.A_T_ndF.topRows(nBoundaryElemUnknowns) + _diffPb.A_ndF_dF.transpose();

		SparseMatrix tmp = _diffPb.A_T_dF.topRows(nBoundaryElemUnknowns).transpose() * _Theta_T_bF_transpose.transpose();
		_S_bF_bF = tmp + _diffPb.A_dF_dF;

		SparseMatrix B_T_T = _diffPb.BiharStab_T_T();
		SparseMatrix B_T_F = _diffPb.BiharStab_T_F();
		SparseMatrix B_F_F = _diffPb.BiharStab_F_F();
		auto B_F_bF = B_F_F.rightCols(HHO->nBoundaryFaces * HHO->nFaceUnknowns);
		auto B_iF_F = B_F_F.topRows(HHO->nInteriorFaces * HHO->nFaceUnknowns);
		auto B_T_iF = B_T_F.leftCols(HHO->nInteriorFaces * HHO->nFaceUnknowns);
		auto B_T_bF = B_T_F.rightCols(HHO->nBoundaryFaces * HHO->nFaceUnknowns);
		_Stab_bF_T = _Theta_T_bF_transpose * B_T_T + B_T_bF.transpose();
		_Stab_bF_F = _Theta_T_bF_transpose * B_T_F + B_F_bF.transpose();

		_Stab_iF_T  = Theta_T_iF_transpose * B_T_T  + B_T_iF.transpose();
		_Stab_iF_F = Theta_T_iF_transpose * B_T_F;
		_Stab_iF_F += B_iF_F;


		/*ExportModule out(Utils::ProgramArgs.OutputDirectory, "", Utils::ProgramArgs.Actions.Export.ValueSeparator);

		cout << Utils::MatrixInfo(Theta_T_bF_transpose, "Theta_T_bF_transpose") << endl;
		out.ExportMatrix(_Theta_T_bF_transpose, "Theta_T_bF_transpose");

		cout << Utils::MatrixInfo(_diffPb.A_T_ndF, "A_T_ndF") << endl;
		out.ExportMatrix(_diffPb.A_T_ndF, "A_T_ndF");

		cout << Utils::MatrixInfo(_diffPb.A_ndF_dF, "A_ndF_dF") << endl;
		out.ExportMatrix(_diffPb.A_ndF_dF, "A_ndF_dF");

		cout << Utils::MatrixInfo(_diffPb.A_dF_dF, "A_dF_dF") << endl;
		out.ExportMatrix(_diffPb.A_dF_dF, "A_dF_dF");

		cout << Utils::MatrixInfo(_S_iF_bF_transpose, "NormalDerStiff_interior") << endl;
		out.ExportMatrix(_S_iF_bF_transpose, "NormalDerStiff_interior");

		cout << Utils::MatrixInfo(_S_bF_bF, "NormalDerStiff_boundary") << endl;
		out.ExportMatrix(_S_bF_bF, "NormalDerStiff_boundary");*/

		_b_fSource = _diffPb.AssembleSourceTerm(_testCase->SourceFunction);
		_g_D = _diffPb.AssembleDirichletTerm(_testCase->DirichletBC.DirichletFunction);
	}

	Vector FindCompatibleTheta() override
	{
		return _diffPb.BoundarySpace.ZeroVector();
	}

	// Solve problem 1 (source=f, Dirich=<dirichlet>)
	Vector Solve1stDiffProblem(const Vector& dirichlet) override
	{
		Vector dummy;
		return Solve1stDiffProblem(dirichlet, false, dummy);
	}

private:
	// Solve problem 1 (source=f, Dirich=<dirichlet>)
	Vector Solve1stDiffProblem(const Vector& dirichlet, bool reconstruct, Vector& reconstruction)
	{
		// Define problem
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

		//if (Utils::ProgramArgs.Actions.Option1 == 1)
			//return _diffPb.ReconstructHigherOrderApproximationFromFaceCoeffs(faceSolution, dirichlet, b_T);
		return cellSolution;
	}

public:
	// Solve problem 1 (source=0, Dirich=<dirichlet>)
	std::array<Vector, 2> Solve1stDiffProblem_Homogeneous(const Vector& dirichlet) //override
	{
		// Define problem
		Vector b_T = _diffPb.ComputeB_T_zeroSource(dirichlet);
		Vector b_ndF = _diffPb.ComputeB_ndF_noNeumann(dirichlet);
		Vector rhs = _diffPb.CondensedRHS(b_T, b_ndF);

		/*if (print)
		{
			Vector v1 = _Theta_T_bF_transpose.transpose()* dirichlet;
			Vector v2 = _diffPb.A_ndF_dF * dirichlet;
			Vector v3 = _S_iF_bF_transpose.transpose() * dirichlet;
			cout << "v1.norm(): " << v1.norm() << endl;
			cout << "v2.norm(): " << v2.norm() << endl;
			cout << "v3.norm(): " << v3.norm() << endl;

			ExportModule out(Utils::ProgramArgs.OutputDirectory, "", Utils::ProgramArgs.Actions.Export.ValueSeparator);

			cout << Utils::MatrixInfo(_diffPb.A_T_dF, "A_T_dF") << endl;
			out.ExportMatrix(_diffPb.A_T_dF, "A_T_dF");

			cout << Utils::MatrixInfo(_diffPb.A_ndF_dF, "A_ndF_dF") << endl;
			out.ExportMatrix(_diffPb.A_ndF_dF, "A_ndF_dF");
			

			cout << "dirichlet.norm(): " << dirichlet.norm() << endl;
			cout << "b_T.norm(): " << b_T.norm() << endl;
			cout << "b_ndF.norm(): " << b_ndF.norm() << endl;
			cout << "rhs.norm(): " << rhs.norm() << endl;
		}*/

		std::array<Vector, 2> hybridSolution;
		auto& cellSolution = hybridSolution[0];
		auto& faceSolution = hybridSolution[1];

		// Solve
		faceSolution = this->_diffSolver->Solve(rhs);
		this->CheckDiffSolverConvergence();

		cellSolution = _diffPb.SolveCellUnknowns(faceSolution, b_T);

		/*if (print)
		{
			//cout << "faceSolution: " << endl << faceSolution.transpose() << endl;
			cout << "faceSolution.norm(): " << _diffPb.NonDirichletFaceSpace.L2Norm(faceSolution) << endl;

			//cout << "cellSolution: " << endl << cellSolution.transpose() << endl;
			cout << "cellSolution.norm(): " << _diffPb.CellSpace.L2Norm(cellSolution) << endl;

			//Vector reconstruct = _diffPb.ReconstructHigherOrderApproximationFromFaceCoeffs(faceSolution, dirichlet, b_T);
			//cout << "reconstruct: " << endl << reconstruct.transpose() << endl;

			double a = _diffPb.a(cellSolution, faceSolution, dirichlet, cellSolution, faceSolution, dirichlet);
			cout << "a = " << a << endl;

			Vector reconstruction = _diffPb.ReconstructHigherOrderApproximationFromFaceCoeffs(faceSolution, dirichlet, b_T);
			//cout << "reconstruction: " << endl << reconstruction.transpose() << endl;
			cout << "L2 norm = " << _diffPb.ReconstructSpace.L2Norm(reconstruction) << endl;
			//ExportModule out(Utils::ProgramArgs.OutputDirectory, "", Utils::ProgramArgs.Actions.Export.ValueSeparator);
			//_diffPb.ExportReconstructedVectorToGMSH(reconstruction, out, "lambda_0");
		}*/

		//if (Utils::ProgramArgs.Actions.Option1 == 1)
			//return _diffPb.ReconstructHigherOrderApproximationFromFaceCoeffs(faceSolution, dirichlet, b_T);
		return hybridSolution;
	}

	// Solve problem 2 (source=<source>, Dirich=g_D)
	Vector Solve2ndDiffProblem(const Vector& source, bool returnBoundaryNormalDerivative = false) override
	{
		// Define problem
		Vector b_source = _diffPb.AssembleSourceTerm(source);
		Vector b_T = _diffPb.ComputeB_T(b_source, _g_D);
		Vector b_ndF = _diffPb.ComputeB_ndF_noNeumann(_g_D);
		Vector rhs = _diffPb.CondensedRHS(b_T, b_ndF);

		// Solve
		Vector faceSolution = this->_diffSolver->Solve(rhs);
		this->CheckDiffSolverConvergence();

		if (returnBoundaryNormalDerivative)
		{
			Vector sourceElemBoundary = _diffPb.ExtractElemBoundary(b_source);
			Vector normalDerivative = _S_iF_bF_transpose * faceSolution + _S_bF_bF * _g_D - _Theta_T_bF_transpose * sourceElemBoundary;
			//return _diffPb.BoundarySpace.SolveMassMatrix(normalDerivative);
			return normalDerivative;
		}
		else
			return _diffPb.ReconstructHigherOrderApproximationFromFaceCoeffs(faceSolution, _g_D, b_T);
	}

	// Solve problem 2 (source=<source>, Dirich=0)
	Vector Solve2ndDiffProblem_Homogeneous(const std::array<Vector, 2>& source, const Vector& boundarySource) //override
	{
		auto& cellSource = source[0];
		auto& faceSource = source[1];

		// Define problem
		Vector b_source = _diffPb.AssembleSourceTerm(cellSource);
		Vector rhs = _diffPb.CondensedRHS_noNeumannZeroDirichlet(b_source);
		// Stabilization
		rhs += _Stab_iF_T * _diffPb.ExtractElemBoundary(cellSource);
		rhs += _Stab_iF_F.leftCols(_diffPb.HHO->nInteriorFaces * _diffPb.HHO->nFaceUnknowns) * faceSource;
		rhs += _Stab_iF_F.rightCols(_diffPb.HHO->nBoundaryFaces * _diffPb.HHO->nFaceUnknowns) * boundarySource;

		// Solve
		Vector faceSolution = this->_diffSolver->Solve(rhs);
		this->CheckDiffSolverConvergence();

		// Normal derivative
		Vector sourceElemBoundary = _diffPb.ExtractElemBoundary(b_source);
		//cout << "_S_iF_bF_transpose * faceSolution: " << endl << (_S_iF_bF_transpose * faceSolution).transpose() << endl;
		//cout << "_Theta_T_bF_transpose * sourceElemBoundary: " << endl << (_Theta_T_bF_transpose * sourceElemBoundary).transpose() << endl;
		Vector normalDerivative = _S_iF_bF_transpose * faceSolution - _Theta_T_bF_transpose * sourceElemBoundary;
		// Stabilization
		normalDerivative -= _Stab_bF_T * _diffPb.ExtractElemBoundary(cellSource);
		normalDerivative -= _Stab_bF_F.leftCols(_diffPb.HHO->nInteriorFaces * _diffPb.HHO->nFaceUnknowns) * faceSource;
		normalDerivative -= _Stab_bF_F.rightCols(_diffPb.HHO->nBoundaryFaces * _diffPb.HHO->nFaceUnknowns) * boundarySource;
		//return _diffPb.BoundarySpace.SolveMassMatrix(normalDerivative);
		return normalDerivative;
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
		Vector lambdaCells = Solve1stDiffProblem(theta, true, lambda);

		// Solve problem 2 (source=<lambda>, Dirich=0)
		solution = Solve2ndDiffProblem(lambdaCells);

		return p;
	}
	
	Vector ProblemOperator(const Vector& x) override
	{
		auto delta = Solve1stDiffProblem_Homogeneous(x);
		return -Solve2ndDiffProblem_Homogeneous(delta, x);
	}

	DenseMatrix Matrix() override
	{
		Vector theta0 = FindCompatibleTheta();
		int n = theta0.rows();
		DenseMatrix A(n, n);

		int nThreads = 1;

		NumberParallelLoop<EmptyResultChunk> parallelLoop(n, nThreads);
		parallelLoop.Execute([this, n, &A](BigNumber i, ParallelChunk<EmptyResultChunk>* chunk)
			{
				Vector e_i = Vector::Zero(n);
				e_i[i] = 1;
				auto lambda = Solve1stDiffProblem_Homogeneous(e_i);
				A.col(i) = -Solve2ndDiffProblem_Homogeneous(lambda, e_i);
			});
		return A;
	}
};
