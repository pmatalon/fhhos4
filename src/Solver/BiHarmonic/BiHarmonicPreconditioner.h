#pragma once
#include "../Preconditioner.h"
#include "../../Discretizations/HHO/BiHarmonicMixedForm_HHO.h"
using namespace std;

template <int Dim>
class BiHarmonicPreconditioner : public Preconditioner
{
private:
	Diffusion_HHO<Dim>& _diffPb;
	vector<Eigen::FullPivLU<DenseMatrix>> _invD;
public:
	BiHarmonicPreconditioner(BiHarmonicMixedForm_HHO<Dim>& biHarPb) :
		_diffPb(biHarPb.DiffPb())
	{}

	void Setup()
	{
		FaceParallelLoop<Dim> parallelLoop(_diffPb._mesh->BoundaryFaces);
		parallelLoop.ReserveChunkCoeffsSize(_diffPb.HHO->nFaceUnknowns * _diffPb.HHO->nFaceUnknowns);

		parallelLoop.Execute([this](Face<Dim>* f, ParallelChunk<CoeffsChunk>* chunk)
			{
				int i = f->Number - _diffPb.HHO->nInteriorFaces;
				Diff_HHOFace<Dim>* face = _diffPb.HHOFace(f);

				Diff_HHOElement<Dim>* elem = _diffPb.HHOElement(f->Element1);

				int nFaceUnknowns = _diffPb.HHO->nFaceUnknowns;
				int nCellUnknowns = _diffPb.HHO->nCellUnknowns;

				int i_loc = elem->LocalNumberOf(face);

				int nElemFaceUnknowns = elem->Faces.size() * nFaceUnknowns;
				for (int k = 0; k < nFaceUnknowns; k++)
				{
					// Problem 1 (Dirichlet 1 at the current unknown, zero source)
					Vector faceCoeffs = Vector::Zero(nElemFaceUnknowns);
					faceCoeffs[i_loc * nFaceUnknowns + k] = 1;
					Vector lambda = elem->ReconstructFromFaces(faceCoeffs);

					// Problem 2 (Dirichlet 0, source)
					Vector b_source = elem->ApplyCellReconstructMassMatrix(lambda);
					Vector cellSolution = elem->AttSolver.solve(b_source);
					// 1st term
					DenseMatrix Atf = elem->A.topRightCorner(nCellUnknowns, nElemFaceUnknowns);
					Vector normalDerivative = Atf.middleCols(i_loc * nFaceUnknowns, nFaceUnknowns).transpose() * cellSolution;
					// 2nd term
					DenseMatrix P = elem->P.middleCols(nCellUnknowns + i_loc * nFaceUnknowns, nFaceUnknowns);
					DenseMatrix m = P.transpose() * elem->MassMatrix(elem->ReconstructionBasis);
					normalDerivative -= m * lambda; 

					// Add minus sign
					chunk->Results.Coeffs.Add(i * nFaceUnknowns, i * nFaceUnknowns + k, -normalDerivative);
				}
			});

		SparseMatrix mat(_diffPb.HHO->nBoundaryFaces * _diffPb.HHO->nFaceUnknowns, _diffPb.HHO->nBoundaryFaces * _diffPb.HHO->nFaceUnknowns);
		parallelLoop.Fill(mat);

		int blockSize = _diffPb.HHO->nFaceUnknowns;
		_invD = vector<Eigen::FullPivLU<DenseMatrix>>(_diffPb.HHO->nBoundaryFaces);
		NumberParallelLoop<EmptyResultChunk> parallelLoop2(_diffPb.HHO->nBoundaryFaces);
		parallelLoop2.Execute([this, &mat, blockSize](BigNumber i, ParallelChunk<EmptyResultChunk>* chunk)
			{
				DenseMatrix Di = mat.block(i * blockSize, i * blockSize, blockSize, blockSize);
				_invD[i].compute(Di);
			});
	}

	MFlops SetupComputationalWork() override
	{
		
	}

	Vector Apply(const Vector& r) override
	{
		int nb = _invD.size();
		int blockSize = _diffPb.HHO->nFaceUnknowns;
		assert(r.rows() == nb * blockSize);

		Vector res(r.rows());
		for (int i = 0; i < nb; i++)
			res.segment(i * blockSize, blockSize) = _invD[i].solve(r.segment(i * blockSize, blockSize));
		return res;
	}

	MFlops SolvingComputationalWork() override
	{
		Utils::FatalError("TODO: SolvingComputationalWork() not implemented");
		return 0;
	}
};