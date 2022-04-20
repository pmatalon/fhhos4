#pragma once
#include "../Preconditioner.h"
#include "../../Discretizations/HHO/BiHarmonicMixedForm_HHO.h"
#include "NeighbourhoodDiffusion_HHO.h"
using namespace std;

template <int Dim>
class BiharPatchPreconditioner : public Preconditioner
{
private:
	Diffusion_HHO<Dim>& _diffPb;
	bool _useIntegrationByParts = true;

	int _neighbourhoodDepth = 1;
	bool _blockDiagPrec = false;
	vector<Eigen::FullPivLU<DenseMatrix>> _invD;
	EigenSparseLU _solver;
public:
	BiharPatchPreconditioner(BiHarmonicMixedFormGlowinski_HHO<Dim>& biHarPb, int neighbourhoodDepth, bool blockDiagPrec = false) :
		_diffPb(biHarPb.DiffPb())
	{
		_useIntegrationByParts = biHarPb.UseIntegrationByParts;
		_neighbourhoodDepth = neighbourhoodDepth;
		_blockDiagPrec = blockDiagPrec;
	}

	void Setup()
	{
		FaceParallelLoop<Dim> parallelLoop(_diffPb._mesh->BoundaryFaces);
		parallelLoop.ReserveChunkCoeffsSize(_diffPb.HHO->nFaceUnknowns * _diffPb.HHO->nFaceUnknowns);

		parallelLoop.Execute([this](Face<Dim>* f, ParallelChunk<CoeffsChunk>* chunk)
			{
				int i = f->Number - _diffPb.HHO->nInteriorFaces;
				Element<Dim>* e = f->Element1;

				int nFaceUnknowns = _diffPb.HHO->nFaceUnknowns;
				int nCellUnknowns = _diffPb.HHO->nCellUnknowns;

				Neighbourhood<Dim> nbh(e, _neighbourhoodDepth);
				NeighbourhoodDiffusion_HHO<Dim> nbhDiff(nbh, _diffPb);


				SparseMatrix PTranspose;
				SparseMatrix ReconstructStiff;
				SparseMatrix ReconstructMass;
				SparseMatrix SolveCellUknTranspose;
				SparseMatrix CellStiff;
				SparseMatrix CellMass;
				SparseMatrix NormalDerivative;
				if (_useIntegrationByParts)
				{
					if (Utils::ProgramArgs.Actions.Option == 0)
					{
						PTranspose = nbhDiff.PTranspose();
						ReconstructStiff = nbhDiff.ReconstructStiffnessMatrix();
						ReconstructMass = nbhDiff.ReconstructMassMatrix();
					}
					else if (Utils::ProgramArgs.Actions.Option == 4)
					{
						PTranspose = nbhDiff.PTranspose();
						ReconstructMass = nbhDiff.ReconstructMassMatrix();
					}
					else if (Utils::ProgramArgs.Actions.Option == 5)
					{
						SolveCellUknTranspose = nbhDiff.SolveCellUknTranspose();
						CellStiff = nbhDiff.CellStiffnessMatrix();
						CellMass = nbhDiff.CellMassMatrix();
					}
				}
				else
					NormalDerivative = nbhDiff.NormalDerivative();
				
				for (int k = 0; k < nFaceUnknowns; k++)
				{
					// Problem 1 (Dirichlet 1 at the current unknown, zero source)
					Vector dirichlet = Vector::Zero(nbh.BoundaryFaces.size() * nFaceUnknowns);
					dirichlet[nbh.BoundaryFaceNumber(f) * nFaceUnknowns + k] = 1;
					// Reconstruct higher-order like in the Glowinski scheme?
					Vector b_T = nbhDiff.ComputeB_T_zeroSource(dirichlet); // TODO don't compute the whole matrix vector product, there is only one 1 in the vector
					Vector b_ndF = nbhDiff.ComputeB_ndF_noNeumann(dirichlet);
					Vector rhs = nbhDiff.CondensedRHS(b_T, b_ndF);

					Vector faceSolution = nbhDiff.FaceUnknownSolver.Solve(rhs);

					Vector lambda;
					if (Utils::ProgramArgs.Actions.Option == 5)
						lambda = nbhDiff.SolveCellUnknowns(faceSolution, b_T);
					else
						lambda = nbhDiff.ReconstructHigherOrder(faceSolution, dirichlet, b_T);

					// Problem 2 (Dirichlet 0, source)
					Vector b_source = nbhDiff.AssembleSourceTerm(lambda);
					rhs = nbhDiff.CondensedRHS_noNeumannZeroDirichlet(b_source);

					faceSolution = nbhDiff.FaceUnknownSolver.Solve(rhs);

					// Normal derivative
					Vector normalDerivative;
					if (_useIntegrationByParts)
					{
						if (Utils::ProgramArgs.Actions.Option == 0)
						{
							Vector reconstruction = nbhDiff.ReconstructHigherOrder(faceSolution, Vector::Zero(nbh.BoundaryFaces.size() * nFaceUnknowns), b_source);
							normalDerivative = nbhDiff.SolveFaceMassMatrixOnBoundary(PTranspose * (ReconstructStiff * reconstruction - ReconstructMass * lambda));
						}
						else if (Utils::ProgramArgs.Actions.Option == 4)
						{
							Vector cellSolution = nbhDiff.SolveCellUnknowns(faceSolution, b_source);

							normalDerivative = nbhDiff.A_T_dF.transpose() * cellSolution + nbhDiff.A_ndF_dF.transpose() * faceSolution;
							normalDerivative -= PTranspose * ReconstructMass * lambda;
						}
						else if (Utils::ProgramArgs.Actions.Option == 5)
						{
							Vector cellSolution = nbhDiff.SolveCellUnknowns(faceSolution, b_source);
							normalDerivative = nbhDiff.SolveFaceMassMatrixOnBoundary(SolveCellUknTranspose * (CellStiff * cellSolution - CellMass * lambda));
						}
						else
							Utils::FatalError("Preconditioner not managed for -opt " + to_string(Utils::ProgramArgs.Actions.Option));
					}
					else
					{
						Vector reconstruction = nbhDiff.ReconstructHigherOrder(faceSolution, Vector::Zero(nbh.BoundaryFaces.size() * nFaceUnknowns), b_source);
						normalDerivative = NormalDerivative * reconstruction;
					}

					for (int i2 = 0; i2 < nbh.BoundaryFaces.size(); i2++)
					{
						Face<Dim>* f2 = nbh.BoundaryFaces[i2];
						if (f2->IsDomainBoundary)
						{
							int j = f2->Number - _diffPb.HHO->nInteriorFaces;
							// Add minus sign
							chunk->Results.Coeffs.Add(j * nFaceUnknowns, i * nFaceUnknowns + k, -normalDerivative.segment(i2 * nFaceUnknowns, nFaceUnknowns));
						}
					}
				}
			});

		SparseMatrix mat(_diffPb.HHO->nBoundaryFaces * _diffPb.HHO->nFaceUnknowns, _diffPb.HHO->nBoundaryFaces * _diffPb.HHO->nFaceUnknowns);
		parallelLoop.Fill(mat);

		//cout << "Preconditioner matrix:" << endl << DenseMatrix(mat) << endl << endl;

		if (_blockDiagPrec)
		{
			int blockSize = _diffPb.HHO->nFaceUnknowns;
			_invD = vector<Eigen::FullPivLU<DenseMatrix>>(_diffPb.HHO->nBoundaryFaces);
			NumberParallelLoop<EmptyResultChunk> parallelLoop2(_diffPb.HHO->nBoundaryFaces);
			parallelLoop2.Execute([this, &mat, blockSize](BigNumber i, ParallelChunk<EmptyResultChunk>* chunk)
				{
					DenseMatrix Di = mat.block(i * blockSize, i * blockSize, blockSize, blockSize);
					_invD[i].compute(Di);
				});
		}
		else
			_solver.Setup(mat);
	}

	MFlops SetupComputationalWork() override
	{
		
	}

	Vector Apply(const Vector& r) override
	{
		if (_blockDiagPrec)
		{
			int nb = _invD.size();
			int blockSize = _diffPb.HHO->nFaceUnknowns;
			assert(r.rows() == nb * blockSize);

			Vector res(r.rows());
			for (int i = 0; i < nb; i++)
				res.segment(i * blockSize, blockSize) = _invD[i].solve(r.segment(i * blockSize, blockSize));
			return res;
		}
		else
			return _solver.Solve(r);
	}

	MFlops SolvingComputationalWork() override
	{
		Utils::FatalError("TODO: SolvingComputationalWork() not implemented");
		return 0;
	}
};