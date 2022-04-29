#pragma once
#include "../Preconditioner.h"
#include "../../Discretizations/HHO/BiHarmonicMixedForm_HHO.h"
#include "../../Discretizations/HHO/ZeroMeanEnforcer.h"
#include "../../Discretizations/HHO/NumericImageEnforcer.h"
#include "NeighbourhoodDiffusion_HHO.h"
using namespace std;

template <int Dim>
class BiharPatchPreconditionerFalk : public Preconditioner
{
private:
	Diffusion_HHO<Dim>& _diffPb;
	bool _reconstructHigherOrderBoundary = false;

	int _neighbourhoodDepth = 1;
	bool _blockDiagPrec = false;
	vector<Eigen::FullPivLU<DenseMatrix>> _invD;
	EigenSparseLU _solver;
public:
	BiharPatchPreconditionerFalk(BiHarmonicMixedFormFalk_HHO<Dim>& biHarPb, int neighbourhoodDepth, bool blockDiagPrec = false) :
		_diffPb(biHarPb.DiffPb())
	{
		_reconstructHigherOrderBoundary = biHarPb.ReconstructHigherOrderBoundary;
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

				Vector b_zeroSource = Vector::Zero(nbh.Elements.size() * nCellUnknowns);
				Vector noDirichlet = Vector();

				ZeroMeanEnforcer integralZeroOnBoundary(_boundarySpace);
				integralZeroOnBoundary.Setup();

				NumericImageEnforcer imageEnforcer(nbhDiff.SkeletonSpace);
				imageEnforcer.Setup();

				ZeroMeanEnforcer integralZeroOnDomain(nbhDiff.ReconstructSpace);
				integralZeroOnDomain.Setup();
				
				for (int k = 0; k < nFaceUnknowns; k++)
				{
					// Problem 1 (Neumann 1 at the current unknown, zero source)
					// ------
					// Define problem
					Vector neumann = Vector::Zero(nbh.BoundaryFaces.size() * nFaceUnknowns);
					neumann[nbh.BoundaryFaceNumber(f) * nFaceUnknowns + k] = 1;
					integralZeroOnBoundary.Enforce(neumann);

					Vector b_neumann = ReconstructHigherOrderBoundary ? _higherOrderBoundary.AssembleNeumannTerm(neumann) : nbhDiff.AssembleNeumannTerm(neumann);
					Vector& rhs = b_neumann;
					// Solve
					imageEnforcer.ProjectOntoImage(rhs); // enforce numerical compatibility
					Vector faceSolution = nbhDiff.FaceUnknownSolver.Solve(rhs);

					// Reconstruct the higher-order polynomial
					Vector lambda = nbhDiff.ReconstructHigherOrderApproximationFromFaceCoeffs(faceSolution, noDirichlet, b_zeroSource);

					// Enforce (lambda|1) = 0
					integralZeroOnDomain.Enforce(lambda);


					// Problem 2 (source=lambda, zero Neumann)
					// ------
					// Define problem
					Vector b_source = nbhDiff.AssembleSourceTerm(lambda);
					Vector rhs2 = nbhDiff.CondensedRHS_noDirichletZeroNeumann(b_source);

					// Solve
					imageEnforcer.ProjectOntoImage(rhs2); // enforce numerical compatibility
					faceSolution = nbhDiff.FaceUnknownSolver.Solve(rhs2);

					// Extract boundary
					Vector boundary;
					if (ReconstructHigherOrderBoundary)
					{
						Vector reconstructedElemBoundary = nbhDiff.ReconstructHigherOrderOnBoundaryOnly(faceSolution, noDirichlet, b_source);
						boundary = _higherOrderBoundary.Trace(reconstructedElemBoundary);
					}
					else
						boundary = faceSolution.tail(nbh.BoundaryFaces.size() * nFaceUnknowns); // keep only the boundary unknowns
					integralZeroOnBoundary.Enforce(boundary);



					for (int i2 = 0; i2 < nbh.BoundaryFaces.size(); i2++)
					{
						Face<Dim>* f2 = nbh.BoundaryFaces[i2];
						if (f2->IsDomainBoundary)
						{
							int j = f2->Number - _diffPb.HHO->nInteriorFaces;
							// Add minus sign
							chunk->Results.Coeffs.Add(j * nFaceUnknowns, i * nFaceUnknowns + k, -boundary.segment(i2 * nFaceUnknowns, nFaceUnknowns));
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