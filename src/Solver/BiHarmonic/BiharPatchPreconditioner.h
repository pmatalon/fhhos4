#pragma once
#include "../Preconditioner.h"
#include "../../Discretizations/HHO/BiHarmonicMixedForm_HHO.h"
#include "NeighbourhoodDiffusion_HHO.h"
#include "BoundaryFacePatch.h"
using namespace std;

template <int Dim>
class BiharPatchPreconditioner : public Preconditioner
{
public:
	enum class Type : unsigned
	{
		SingleFaceNeighbourhood,
		FacePatchNeighbourhood
	};
private:
	Diffusion_HHO<Dim>& _diffPb;

	Type _type;
	int _facePatchSize = 3;
	int _neighbourhoodDepth = 1;
	bool _blockDiagPrec = false;
	vector<Eigen::FullPivLU<DenseMatrix>> _invD;
	EigenSparseLU _solver;
public:
	BiharPatchPreconditioner(BiHarmonicMixedFormGlowinski_HHO<Dim>& biHarPb, Type type, int facePatchSize, int neighbourhoodDepth, bool blockDiagPrec = false) :
		_diffPb(biHarPb.DiffPb())
	{
		_type = type;
		_facePatchSize = facePatchSize;
		_neighbourhoodDepth = neighbourhoodDepth;
		_blockDiagPrec = blockDiagPrec;
	}

	void Serialize(ostream& os) const override
	{
		if (_type == Type::SingleFaceNeighbourhood)
		{
			if (_blockDiagPrec)
				os << "diagonal ";
			os << "single face neighbourhood (neighbourhood depth = " << _neighbourhoodDepth << ")";
		}
		else if (_type == Type::FacePatchNeighbourhood)
		{
			if (_blockDiagPrec)
				os << "diagonal ";
			os << "face patch neighbourhood (patch size = " << _facePatchSize << ", neighbourhood depth = " << _neighbourhoodDepth << ")";
		}
	}

	void Setup(const DenseMatrix& A) override
	{
		if (_type == Type::SingleFaceNeighbourhood)
			Setup1(A);
		else if (_type == Type::FacePatchNeighbourhood)
			Setup2(A);
	}

private:
	void Setup1(const DenseMatrix& A)
	{
		ExportModule out(Utils::ProgramArgs.OutputDirectory, "", Utils::ProgramArgs.Actions.Export.ValueSeparator);

		FaceParallelLoop<Dim> parallelLoop(_diffPb._mesh->BoundaryFaces);
		parallelLoop.ReserveChunkCoeffsSize(_diffPb.HHO->nFaceUnknowns * _diffPb.HHO->nFaceUnknowns);

		parallelLoop.Execute([this, &out](Face<Dim>* f, ParallelChunk<CoeffsChunk>* chunk)
			{
				int i = f->Number - _diffPb.HHO->nInteriorFaces;
				Element<Dim>* e = f->Element1;

				int nFaceUnknowns = _diffPb.HHO->nFaceUnknowns;

				Neighbourhood<Dim> nbh(e, _neighbourhoodDepth);
				NeighbourhoodDiffusion_HHO<Dim> nbhDiff(nbh, _diffPb);

				SparseMatrix SolveCellUknTranspose = nbhDiff.SolveCellUknTranspose();
				SparseMatrix CellMass = nbhDiff.CellMassMatrix();

				SparseMatrix normalDerStiff = SolveCellUknTranspose * nbhDiff.A_T_ndF + nbhDiff.A_ndF_dF.transpose();
				SparseMatrix normalDerMass = SolveCellUknTranspose * CellMass;
				
				for (int k = 0; k < nFaceUnknowns; k++)
				{
					// Dirichlet 1 at the current unknown
					Vector dirichlet = Vector::Zero(nbh.BoundaryFaces.size() * nFaceUnknowns);
					dirichlet[nbh.BoundaryFaceNumber(f) * nFaceUnknowns + k] = 1;

					Vector approxColumn = BiharOperator(dirichlet, nbhDiff, normalDerStiff, normalDerMass);

					for (int i2 = 0; i2 < nbh.BoundaryFaces.size(); i2++)
					{
						Face<Dim>* f2 = nbh.BoundaryFaces[i2];
						if (f2->IsDomainBoundary)
						{
							int j = f2->Number - _diffPb.HHO->nInteriorFaces;
							chunk->Results.Coeffs.Add(j * nFaceUnknowns, i * nFaceUnknowns + k, approxColumn.segment(i2 * nFaceUnknowns, nFaceUnknowns));
						}
					}
				}
			});

		SparseMatrix mat(_diffPb.HHO->nBoundaryFaces * _diffPb.HHO->nFaceUnknowns, _diffPb.HHO->nBoundaryFaces * _diffPb.HHO->nFaceUnknowns);
		parallelLoop.Fill(mat);

		//SparseMatrix matT = mat.triangularView<Eigen::StrictlyLower>().transpose();
		//mat = mat.triangularView<Eigen::Lower>() + matT;

		if (Utils::ProgramArgs.Actions.PrintDebug)
		{
			if (Dim <= 2)
				cout << "Preconditioner matrix:" << endl << DenseMatrix(mat) << endl << endl;
			else
				cout << "Preconditioner matrix:" << endl << DenseMatrix(mat.topLeftCorner(3 * _diffPb.HHO->nFaceUnknowns, 3 * _diffPb.HHO->nFaceUnknowns)) << endl << endl;
			cout << "||A-mat|| = " << (A - DenseMatrix(mat)).norm() << endl;
		}

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

	void Setup2(const DenseMatrix& A)
	{
		vector<BoundaryFacePatch<Dim>> patches;
		patches.reserve(_diffPb._mesh->BoundaryFaces.size() / _facePatchSize);

		auto nInteriorFaces = _diffPb.HHO->nInteriorFaces;
		vector<bool> isFaceInAPatch(_diffPb.HHO->nBoundaryFaces);

		//vector<Element<Dim>*> remainingBoundaryFaces = _diffPb._mesh->BoundaryFaces;
		//random_shuffle(remainingBoundaryFaces.begin(), remainingBoundaryFaces.end());

		for (Face<Dim>* f : _diffPb._mesh->BoundaryFaces)
		{
			if (!isFaceInAPatch[f->Number - nInteriorFaces])
				patches.emplace_back(f, isFaceInAPatch, nInteriorFaces, _facePatchSize);
		}

		if (Utils::ProgramArgs.Actions.PrintDebug)
		{
			for (auto p : patches)
			{
				cout << "[";
				for (auto f : p.Faces)
					cout << f->Number - _diffPb.HHO->nInteriorFaces << " ";
				cout << "]" << endl;
			}
		}

		ParallelLoop<BoundaryFacePatch<Dim>, CoeffsChunk> parallelLoop(patches);
		parallelLoop.Execute([this](BoundaryFacePatch<Dim> patch, ParallelChunk<CoeffsChunk>* chunk)
			{
				int nFaceUnknowns = _diffPb.HHO->nFaceUnknowns;

				vector<Element<Dim>*> boundaryElements;
				for (Face<Dim>* f : patch.Faces)
					boundaryElements.push_back(f->Element1);

				std::sort(boundaryElements.begin(), boundaryElements.end());
				std::unique(boundaryElements.begin(), boundaryElements.end());

				Neighbourhood<Dim> nbh(boundaryElements, _neighbourhoodDepth);
				NeighbourhoodDiffusion_HHO<Dim> nbhDiff(nbh, _diffPb);

				SparseMatrix SolveCellUknTranspose = nbhDiff.SolveCellUknTranspose();
				SparseMatrix CellMass = nbhDiff.CellMassMatrix();

				SparseMatrix normalDerStiff = SolveCellUknTranspose * nbhDiff.A_T_ndF + nbhDiff.A_ndF_dF.transpose();
				SparseMatrix normalDerMass = SolveCellUknTranspose * CellMass;

				for (int i1 = 0; i1 < patch.Faces.size(); ++i1)
				{
					Face<Dim>* f1 = patch.Faces[i1];
					int globNum1 = f1->Number - _diffPb.HHO->nInteriorFaces;
					int locNum1 = nbh.BoundaryFaceNumber(f1);

					for (int k = 0; k < nFaceUnknowns; k++)
					{
						// Dirichlet 1 at the current unknown
						Vector dirichlet = Vector::Zero(nbh.BoundaryFaces.size() * nFaceUnknowns);
						dirichlet[locNum1 * nFaceUnknowns + k] = 1;

						Vector approxColumn = BiharOperator(dirichlet, nbhDiff, normalDerStiff, normalDerMass);

						//cout << approxColumn.segment(locNum1 * nFaceUnknowns, nFaceUnknowns) << endl << endl;

						chunk->Results.Coeffs.Add(globNum1 * nFaceUnknowns, globNum1 * nFaceUnknowns + k, approxColumn.segment(locNum1 * nFaceUnknowns, nFaceUnknowns));

						/*for (int i2 = i1 + 1; i2 < patch.Faces.size(); ++i2)
						{
							Face<Dim>* f2 = patch.Faces[i2];
							int globNum2 = f2->Number - _diffPb.HHO->nInteriorFaces;
							int locNum2 = nbh.BoundaryFaceNumber(f2);

							chunk->Results.Coeffs.Add(globNum2 * nFaceUnknowns, globNum1 * nFaceUnknowns + k, approxColumn.segment(locNum2 * nFaceUnknowns, nFaceUnknowns));
							//chunk->Results.Coeffs.Add(globNum1 * nFaceUnknowns + k, globNum2 * nFaceUnknowns, approxColumn.segment(locNum2 * nFaceUnknowns, nFaceUnknowns).transpose());
							chunk->Results.Coeffs.Add(globNum1 * nFaceUnknowns, globNum2 * nFaceUnknowns + k, approxColumn.segment(locNum2 * nFaceUnknowns, nFaceUnknowns));
						}*/
						for (int i2 = 0; i2 < patch.Faces.size(); ++i2)
						{
							if (i1 == i2)
								continue;

							Face<Dim>* f2 = patch.Faces[i2];
							int globNum2 = f2->Number - _diffPb.HHO->nInteriorFaces;
							int locNum2 = nbh.BoundaryFaceNumber(f2);

							chunk->Results.Coeffs.Add(globNum2 * nFaceUnknowns, globNum1 * nFaceUnknowns + k, approxColumn.segment(locNum2 * nFaceUnknowns, nFaceUnknowns));
							//chunk->Results.Coeffs.Add(globNum1 * nFaceUnknowns + k, globNum2 * nFaceUnknowns, approxColumn.segment(locNum2 * nFaceUnknowns, nFaceUnknowns).transpose());
							//chunk->Results.Coeffs.Add(globNum1 * nFaceUnknowns, globNum2 * nFaceUnknowns + k, approxColumn.segment(locNum2 * nFaceUnknowns, nFaceUnknowns));
						}
					}
				}
			});

		SparseMatrix mat(_diffPb.HHO->nBoundaryFaces * _diffPb.HHO->nFaceUnknowns, _diffPb.HHO->nBoundaryFaces * _diffPb.HHO->nFaceUnknowns);
		parallelLoop.Fill(mat);

		if (Utils::ProgramArgs.Actions.PrintDebug)
		{
			if (Dim <= 2)
				cout << "Preconditioner matrix:" << endl << DenseMatrix(mat) << endl << endl;
			else
				cout << "Preconditioner matrix:" << endl << DenseMatrix(mat.topLeftCorner(3 * _diffPb.HHO->nFaceUnknowns, 3 * _diffPb.HHO->nFaceUnknowns)) << endl << endl;
			cout << "Condition number = " << Utils::Cond(DenseMatrix(mat)) << endl;
			cout << "||A-mat|| = " << (A - DenseMatrix(mat)).norm() << endl;
		}

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

	Vector BiharOperator(const Vector& dirichlet, NeighbourhoodDiffusion_HHO<Dim>& nbhDiff, const SparseMatrix& normalDerStiff, const SparseMatrix& normalDerMass)
	{
		// Problem 1 (Dirichlet, zero source)
		Vector b_T = nbhDiff.ComputeB_T_zeroSource(dirichlet); // TODO don't compute the whole matrix vector product, there is only one 1 in the vector
		Vector b_ndF = nbhDiff.ComputeB_ndF_noNeumann(dirichlet);
		Vector rhs = nbhDiff.CondensedRHS(b_T, b_ndF);

		Vector faceSolutionLambda = nbhDiff.FaceUnknownSolver.Solve(rhs);

		Vector lambda = nbhDiff.SolveCellUnknowns(faceSolutionLambda, b_T);

		// Problem 2 (Dirichlet 0, source)
		Vector b_source = nbhDiff.AssembleSourceTerm(lambda);
		rhs = nbhDiff.CondensedRHS_noNeumannZeroDirichlet(b_source);

		Vector faceSolution = nbhDiff.FaceUnknownSolver.Solve(rhs);


		// To generate the figures:
		// Neighbouhood (nbh): -pb bihar -geo square -source poly -mesh tri -n 8 -k 3 -nbh-depth 3 
		// Domain       (dom): idem with -nbh-depth 30 
		//
		// In GMSH, 
		// Open the square Geometry, go to Tools->Options->Geometry
		// - Visibility tab: untick Points
		// - Aspect tab: set 1.0 for curve width
		// - Color tab: set Curves to black
		// Open all 4 .pos file, go to Tools->Options. For each view, one by one:
		// - General tab: choose the "Filled iso-values" interval type
		// - Map tab: type 9 to select the B/W color map, then type i to inverse the colors
		// - In order to get the same scale for u_nbh and u_dom, select the view of u_nbh, tab General,
		//   select the "Custom" range mode, and enter -4.85e-5 (Min) and 0.00394 (Max)
		// - Visibility tab, untick "show value scale"
		/*if (k == 0 && f->Center().X == 0)
		{
			lambda = nbhDiff.ReconstructHigherOrder(faceSolutionLambda, dirichlet, b_T);
			Vector u = nbhDiff.ReconstructHigherOrder(faceSolution, Vector::Zero(nbh.BoundaryFaces.size() * nFaceUnknowns), b_source);
			_diffPb.ExportReconstructedVectorToGMSH(NeighbourhoodToDomain(nbh, lambda), out, "lambda_prec");
			_diffPb.ExportReconstructedVectorToGMSH(NeighbourhoodToDomain(nbh, u), out, "u_prec");
		}*/

		// Normal derivative
		Vector normalDerivative = normalDerStiff * faceSolution - normalDerMass * lambda;
		//return -nbhDiff.SolveFaceMassMatrixOnBoundary(normalDerivative);
		return -normalDerivative;
	}

public:
	MFlops SetupComputationalWork() override
	{
		return 0;
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

private:
	Vector NeighbourhoodToDomain(const Neighbourhood<Dim>& nbh, const Vector& v_nbh)
	{
		Vector v_domain = Vector::Zero(_diffPb.HHO->nTotalReconstructUnknowns);
		int nReconstructUnknowns = _diffPb.HHO->nReconstructUnknowns;
		for (int i = 0; i < nbh.Elements.size(); i++)
		{
			Element<Dim>* e = nbh.Elements[i];
			v_domain.segment(e->Number * nReconstructUnknowns, nReconstructUnknowns) = v_nbh.segment(i * nReconstructUnknowns, nReconstructUnknowns);
		}
		return v_domain;
	}
};