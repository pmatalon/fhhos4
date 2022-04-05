#pragma once
#include "../Preconditioner.h"
#include "../../Discretizations/HHO/BiHarmonicMixedForm_HHO.h"
#include "NeighbourhoodDiffusion_HHO.h"
using namespace std;

template <int Dim>
class NeighbourhoodBiHarmonicPreconditioner : public Preconditioner
{
private:
	Diffusion_HHO<Dim>& _diffPb;

	bool _blockDiagPrec = false;
	vector<Eigen::FullPivLU<DenseMatrix>> _invD;
	EigenSparseLU _solver;
public:
	NeighbourhoodBiHarmonicPreconditioner(BiHarmonicMixedForm_HHO<Dim>& biHarPb, bool blockDiagPrec = false) :
		_diffPb(biHarPb.DiffPb())
	{
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

				Neighbourhood<Dim> nbh(e);
				NeighbourhoodDiffusion_HHO<Dim> nbhDiff(nbh, _diffPb);
				
				for (int k = 0; k < nFaceUnknowns; k++)
				{
					// Problem 1 (Dirichlet 1 at the current unknown, zero source)
					Vector dirichlet = Vector::Zero(nbh.BoundaryFaces.size() * nFaceUnknowns);
					dirichlet[nbh.BoundaryFaceNumber(f) * nFaceUnknowns + k] = 1;
					Vector b_T = nbhDiff.ComputeB_T_zeroSource(dirichlet); // TODO don't compute the whole matrix vector product, there is only one 1 in the vector
					Vector b_ndF = nbhDiff.ComputeB_ndF_noNeumann(dirichlet);
					Vector rhs = nbhDiff.CondensedRHS(b_T, b_ndF);

					Vector faceSolution = nbhDiff.FaceUnknownSolver.Solve(rhs);

					Vector lambda = nbhDiff.ReconstructHigherOrder(faceSolution, dirichlet, b_T);

					// Problem 2 (Dirichlet 0, source)
					Vector b_source = nbhDiff.AssembleSourceTerm(lambda);
					rhs = nbhDiff.CondensedRHS_noNeumannZeroDirichlet(b_source);

					faceSolution = nbhDiff.FaceUnknownSolver.Solve(rhs);

					Vector cellSolution = nbhDiff.SolveCellUnknowns(faceSolution, b_source);

					// Normal derivative
					Vector normalDerivative = nbhDiff.A_T_dF.transpose() * cellSolution + nbhDiff.A_ndF_dF.transpose() * faceSolution;
					normalDerivative -= nbhDiff.PTranspose_Mass * lambda;
					for (int i2 = 0; i2 < nbh.BoundaryFaces.size(); i2++)
					{
						Face<Dim>* f2 = nbh.BoundaryFaces[i2];
						if (f2->IsDomainBoundary)
						{
							int j = f2->Number - _diffPb.HHO->nInteriorFaces;
							// Add minus sign
							chunk->Results.Coeffs.Add(j * nFaceUnknowns, i * nFaceUnknowns + k, -normalDerivative.segment(i2 * nFaceUnknowns, nFaceUnknowns));
							if (f2 != f)
								chunk->Results.Coeffs.Add(i * nFaceUnknowns + k, j * nFaceUnknowns, -normalDerivative.segment(i2 * nFaceUnknowns, nFaceUnknowns).transpose());
						}
					}
				}
			});

		SparseMatrix mat(_diffPb.HHO->nBoundaryFaces * _diffPb.HHO->nFaceUnknowns, _diffPb.HHO->nBoundaryFaces * _diffPb.HHO->nFaceUnknowns);
		parallelLoop.Fill(mat);

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