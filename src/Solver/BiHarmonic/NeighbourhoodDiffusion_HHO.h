#pragma once
#include "../../Discretizations/HHO/Diffusion_HHO.h"
#include "../../Mesh/Neighbourhood.h"
#include "../Direct/EigenSparseCholesky.h"
using namespace std;

template <int Dim>
class NeighbourhoodDiffusion_HHO
{
private:
	const Neighbourhood<Dim>& _nbh;
	Diffusion_HHO<Dim>& _diffPb;
	bool dirichlet = true;
public:
	SparseMatrix A;
	SparseMatrix A_T_dF;
	SparseMatrix A_T_ndF;
	SparseMatrix A_ndF_dF;

	HHOSkeletonSpace<Dim> SkeletonSpace;
	HHOReconstructSpace<Dim> ReconstructSpace;
	HHOBoundarySpace<Dim> BoundarySpace;

	EigenSparseCholesky FaceUnknownSolver;

	NeighbourhoodDiffusion_HHO(const Neighbourhood<Dim>& nbh, Diffusion_HHO<Dim>& diffPb) :
		_nbh(nbh),
		_diffPb(diffPb)
	{
		Init();
	}

private:
	void Init()
	{
		int nFaceUnknowns = _diffPb.HHO->nFaceUnknowns;
		int nCellUnknowns = _diffPb.HHO->nCellUnknowns;

		// Non-Dirichlet faces
		const vector<Face<Dim>*>& nDFaces = dirichlet ? _nbh.InteriorFaces : _nbh.Faces;

		// A

		A = SparseMatrix(nDFaces.size() * nFaceUnknowns, nDFaces.size() * nFaceUnknowns);
		NonZeroCoefficients coeffs(nDFaces.size() * nFaceUnknowns * nFaceUnknowns);
		for (int i = 0; i < nDFaces.size(); i++)
		{
			Face<Dim>* fi = nDFaces[i];
			DenseMatrix block = _diffPb.A.block(fi->Number * nFaceUnknowns, fi->Number * nFaceUnknowns, nFaceUnknowns, nFaceUnknowns);
			coeffs.Add(i * nFaceUnknowns, i * nFaceUnknowns, block);
			for (int j = i+1; j < nDFaces.size(); j++)
			{
				Face<Dim>* fj = nDFaces[j];
				block = _diffPb.A.block(fi->Number * nFaceUnknowns, fj->Number * nFaceUnknowns, nFaceUnknowns, nFaceUnknowns);
				coeffs.Add(i * nFaceUnknowns, j * nFaceUnknowns, block);
				coeffs.Add(j * nFaceUnknowns, i * nFaceUnknowns, block.transpose());
			}
		}
		coeffs.Fill(A);

		FaceUnknownSolver.Setup(A);

		if (dirichlet)
		{
			// A_T_dF
			A_T_dF = this->A_T_dF_Matrix();
		}

		// A_T_ndF
		A_T_ndF = this->A_T_ndF_Matrix();


		if (dirichlet)
		{
			// A_ndF_dF
			A_ndF_dF = this->A_ndF_dF_Matrix();
		}

		SkeletonSpace    = HHOSkeletonSpace   (&_nbh, _diffPb.HHO, _diffPb._hhoFaces);
		ReconstructSpace = HHOReconstructSpace(&_nbh, _diffPb.HHO, _diffPb._hhoElements);
		BoundarySpace    = HHOBoundarySpace   (&_nbh, _diffPb.HHO, _diffPb._hhoFaces);
	}


public:
	const Neighbourhood<Dim>& Nbh() const
	{
		return _nbh;
	}

	SparseMatrix A_T_dF_Matrix()
	{
		int nCellUnknowns = _diffPb.HHO->nCellUnknowns;
		int nFaceUnknowns = _diffPb.HHO->nFaceUnknowns;

		SparseMatrix A_T_dF(_nbh.Elements.size() * nCellUnknowns, _nbh.BoundaryFaces.size() * nFaceUnknowns);
		NonZeroCoefficients coeffs;
		for (int i = 0; i < _nbh.Elements.size(); i++)
		{
			Element<Dim>* e = _nbh.Elements[i];
			for (int j = 0; j < _nbh.BoundaryFaces.size(); j++)
			{
				Face<Dim>* f = _nbh.BoundaryFaces[j];
				if (f->IsIn(e->Faces))
				{
					DenseMatrix block;
					if (f->IsDomainBoundary)
						block = _diffPb.A_T_dF.block(e->Number * nCellUnknowns, (f->Number - _diffPb.HHO->nInteriorFaces) * nFaceUnknowns, nCellUnknowns, nFaceUnknowns);
					else
						block = _diffPb.A_T_ndF.block(e->Number * nCellUnknowns, f->Number * nFaceUnknowns, nCellUnknowns, nFaceUnknowns);
					coeffs.Add(i * nCellUnknowns, j * nFaceUnknowns, block);
				}
			}
		}
		coeffs.Fill(A_T_dF);
		return A_T_dF;
	}

	SparseMatrix A_T_ndF_Matrix()
	{
		int nCellUnknowns = _diffPb.HHO->nCellUnknowns;
		int nFaceUnknowns = _diffPb.HHO->nFaceUnknowns;

		const vector<Face<Dim>*>& nDFaces = dirichlet ? _nbh.InteriorFaces : _nbh.Faces;

		SparseMatrix A_T_ndF(_nbh.Elements.size() * nCellUnknowns, nDFaces.size() * nFaceUnknowns);
		NonZeroCoefficients coeffs;
		for (int i = 0; i < _nbh.Elements.size(); i++)
		{
			Element<Dim>* e = _nbh.Elements[i];
			for (int j = 0; j < nDFaces.size(); j++)
			{
				Face<Dim>* f = nDFaces[j];
				if (f->IsIn(e->Faces))
				{
					DenseMatrix block = _diffPb.A_T_ndF.block(e->Number * nCellUnknowns, f->Number * nFaceUnknowns, nCellUnknowns, nFaceUnknowns);
					coeffs.Add(i * nCellUnknowns, j * nFaceUnknowns, block);
				}
			}
		}
		coeffs.Fill(A_T_ndF);
		return A_T_ndF;
	}

	SparseMatrix A_ndF_dF_Matrix()
	{
		int nFaceUnknowns = _diffPb.HHO->nFaceUnknowns;

		SparseMatrix A_ndF_dF(_nbh.InteriorFaces.size() * nFaceUnknowns, _nbh.BoundaryFaces.size() * nFaceUnknowns);
		NonZeroCoefficients coeffs;
		for (int i = 0; i < _nbh.InteriorFaces.size(); i++)
		{
			Face<Dim>* fi = _nbh.InteriorFaces[i];
			for (int j = 0; j < _nbh.BoundaryFaces.size(); j++)
			{
				Face<Dim>* fj = _nbh.BoundaryFaces[j];

				DenseMatrix block;
				if (fj->IsDomainBoundary)
					block = _diffPb.A_ndF_dF.block(fi->Number * nFaceUnknowns, (fj->Number - _diffPb.HHO->nInteriorFaces) * nFaceUnknowns, nFaceUnknowns, nFaceUnknowns);
				else
					block = _diffPb.A_ndF_ndF.block(fi->Number * nFaceUnknowns, fj->Number * nFaceUnknowns, nFaceUnknowns, nFaceUnknowns);
				coeffs.Add(i * nFaceUnknowns, j * nFaceUnknowns, block);
			}
		}
		coeffs.Fill(A_ndF_dF);
		return A_ndF_dF;
	}

	//         [    0    ]   <---- interior faces
	// returns [---------]
	//         [ b_neumF ]   <---- Neumann faces
	// 
	// where b_neumF = [ M * coeffs ]
	/*Vector AssembleNeumannTerm(const Vector& neumannFuncCoeffs)
	{
		int nFaceUnknowns = _diffPb.HHO->nFaceUnknowns;
		assert(neumannFuncCoeffs.rows() == _nbh->BoundaryFaces.size() * nFaceUnknowns);

		Vector b_neumann = Vector::Zero(_nbh.Faces.size() * nFaceUnknowns);
		for (int i = 0; i < _nbh.BoundaryFaces.size(); i++)
		{
			Face<Dim>* f = _nbh.BoundaryFaces[i];
			Diff_HHOFace<Dim>* face = _diffPb.HHOFace(f);
			BigNumber i = face->Number() * nFaceUnknowns;
			assert(i >= _nbh->InteriorFaces.size() * nFaceUnknowns);
			BigNumber j = (face->Number() - _nbh->InteriorFaces.size()) * nFaceUnknowns;
			b_neumann.segment(i, nFaceUnknowns) = face->ApplyMassMatrix(neumannFuncCoeffs.segment(j, nFaceUnknowns));
		}
		return b_neumann;
	}*/

	SparseMatrix Theta_T_bF_transpose()
	{
		int nFaceUnknowns = _diffPb.HHO->nFaceUnknowns;
		int nCellUnknowns = _diffPb.HHO->nCellUnknowns;

		SparseMatrix mat(_nbh.BoundaryFaces.size() * nFaceUnknowns, _nbh.Elements.size() * nCellUnknowns);
		NonZeroCoefficients coeffs;
		for (int i = 0; i < _nbh.BoundaryFaces.size(); i++)
		{
			Face<Dim>* f = _nbh.BoundaryFaces[i];
			if (f->IsDomainBoundary)
			{
				Element<Dim>* e = f->Element1;
				DenseMatrix P = _diffPb.HHOElement(e)->SolveCellUnknownsMatrix().middleCols(e->LocalNumberOf(f) * nFaceUnknowns, nFaceUnknowns);

				coeffs.Add(i * nFaceUnknowns, _nbh.ElementNumber(e) * nCellUnknowns, P.transpose());
			}
		}
		coeffs.Fill(mat);
		return mat;
	}

	SparseMatrix Theta_T_F_transpose()
	{
		int nFaceUnknowns = _diffPb.HHO->nFaceUnknowns;
		int nCellUnknowns = _diffPb.HHO->nCellUnknowns;

		SparseMatrix mat(_nbh.Faces.size() * nFaceUnknowns, _nbh.Elements.size() * nCellUnknowns);
		NonZeroCoefficients coeffs;
		for (int j = 0; j < _nbh.Elements.size(); j++)
		{
			Element<Dim>* e = _nbh.Elements[j];
			if (!e->IsOnBoundary())
				continue;
			DenseMatrix Theta = _diffPb.HHOElement(e)->SolveCellUnknownsMatrix();
			for (auto f : e->Faces)
			{
				DenseMatrix Theta_Tf = Theta.middleCols(e->LocalNumberOf(f) * nFaceUnknowns, nFaceUnknowns);
				coeffs.Add(_nbh.FaceNumber(f) * nFaceUnknowns, j * nCellUnknowns, Theta_Tf.transpose());
			}
		}
		coeffs.Fill(mat);
		return mat;
	}

	SparseMatrix BiharStab_T_T()
	{
		int nFaceUnknowns = _diffPb.HHO->nFaceUnknowns;
		int nCellUnknowns = _diffPb.HHO->nCellUnknowns;

		SparseMatrix mat(_nbh.Elements.size() * nCellUnknowns, _nbh.Elements.size() * nCellUnknowns);
		NonZeroCoefficients coeffs;
		for (int i = 0; i < _nbh.Elements.size(); i++)
		{
			Element<Dim>* e = _nbh.Elements[i];
			if (!e->IsOnBoundary())
				continue;
			auto hhoElem = _diffPb.HHOElement(e);
			DenseMatrix Stab_T_T = hhoElem->Astab.topLeftCorner(nCellUnknowns, nCellUnknowns);
			coeffs.Add(i * nCellUnknowns, i * nCellUnknowns, Stab_T_T);
		}
		coeffs.Fill(mat);
		return mat;
	}

	SparseMatrix BiharStab_T_F()
	{
		int nFaceUnknowns = _diffPb.HHO->nFaceUnknowns;
		int nCellUnknowns = _diffPb.HHO->nCellUnknowns;

		SparseMatrix mat(_nbh.Elements.size() * nCellUnknowns, _nbh.Faces.size() * nFaceUnknowns);
		NonZeroCoefficients coeffs;
		for (int i = 0; i < _nbh.Elements.size(); i++)
		{
			Element<Dim>* e = _nbh.Elements[i];
			if (!e->IsOnBoundary())
				continue;
			auto hhoElem = _diffPb.HHOElement(e);
			for (auto f : hhoElem->Faces)
			{
				coeffs.AddBlock(i * nCellUnknowns, _nbh.FaceNumber(f->MeshFace) * nFaceUnknowns, hhoElem->Astab, 0, hhoElem->FirstDOFNumber(f), nCellUnknowns, nFaceUnknowns);
			}
		}
		coeffs.Fill(mat);
		return mat;
	}

	SparseMatrix BiharStab_F_F()
	{
		int nFaceUnknowns = _diffPb.HHO->nFaceUnknowns;
		int nCellUnknowns = _diffPb.HHO->nCellUnknowns;

		SparseMatrix mat(_nbh.Faces.size() * nFaceUnknowns, _nbh.Faces.size() * nFaceUnknowns);
		NonZeroCoefficients coeffs;
		for (int i = 0; i < _nbh.Elements.size(); i++)
		{
			Element<Dim>* e = _nbh.Elements[i];
			if (!e->IsOnBoundary())
				continue;
			auto hhoElem = _diffPb.HHOElement(e);
			for (auto f1 : hhoElem->Faces)
			{
				for (auto f2 : hhoElem->Faces)
				{
					coeffs.AddBlock(_nbh.FaceNumber(f1->MeshFace) * nFaceUnknowns, _nbh.FaceNumber(f2->MeshFace) * nFaceUnknowns, hhoElem->Astab, hhoElem->FirstDOFNumber(f1), hhoElem->FirstDOFNumber(f2), nFaceUnknowns, nFaceUnknowns);
				}
			}
		}
		coeffs.Fill(mat);
		return mat;
	}

	/*SparseMatrix CellMassMatrix()
	{
		int nCellUnknowns = _diffPb.HHO->nCellUnknowns;

		SparseMatrix mat(_nbh.Elements.size() * nCellUnknowns, _nbh.Elements.size() * nCellUnknowns);
		NonZeroCoefficients coeffs;
		for (int i = 0; i < _nbh.Elements.size(); i++)
		{
			Diff_HHOElement<Dim>* e = _diffPb.HHOElement(_nbh.Elements[i]);
			coeffs.Add(i * nCellUnknowns, i * nCellUnknowns, e->MassMatrix(e->CellBasis));
		}
		coeffs.Fill(mat);
		return mat;
	}*/


	Vector AssembleSourceTerm(const Vector& sourceFuncCoeffs)
	{
		Vector b_source = Vector(_nbh.Elements.size() * _diffPb.HHO->nCellUnknowns);

		if (sourceFuncCoeffs.rows() == _nbh.Elements.size() * _diffPb.HHO->nReconstructUnknowns) // degree k+1
		{
			for (int i = 0; i < _nbh.Elements.size(); i++)
			{
				Diff_HHOElement<Dim>* e = _diffPb.HHOElement(_nbh.Elements[i]);
				b_source.segment(i * _diffPb.HHO->nCellUnknowns, _diffPb.HHO->nCellUnknowns) = e->ApplyCellReconstructMassMatrix(sourceFuncCoeffs.segment(i * _diffPb.HHO->nReconstructUnknowns, _diffPb.HHO->nReconstructUnknowns));
			}
		}
		else if (sourceFuncCoeffs.rows() == _nbh.Elements.size() * _diffPb.HHO->nCellUnknowns)
		{
			for (int i = 0; i < _nbh.Elements.size(); i++)
			{
				Diff_HHOElement<Dim>* e = _diffPb.HHOElement(_nbh.Elements[i]);
				b_source.segment(i * _diffPb.HHO->nCellUnknowns, _diffPb.HHO->nCellUnknowns) = e->ApplyCellMassMatrix(sourceFuncCoeffs.segment(i * _diffPb.HHO->nCellUnknowns, _diffPb.HHO->nCellUnknowns));
			}
		}
		else
			Utils::FatalError("sourceFuncCoeffs has not the right size");

		return b_source;
	}

	// Right-hand side

	Vector ComputeB_T_zeroSource(const Vector& x_dF) const
	{
		assert(x_dF.rows() == _nbh.BoundaryFaces.size() * _diffPb.HHO->nFaceUnknowns);

		return -A_T_dF * x_dF;
	}

	Vector ComputeB_ndF_noNeumann(const Vector& x_dF) const
	{
		assert(x_dF.rows() == _nbh.BoundaryFaces.size() * _diffPb.HHO->nFaceUnknowns);

		return -A_ndF_dF * x_dF;
	}

	Vector CondensedRHS(const Vector& b_T, const Vector& b_ndF) const
	{
		return b_ndF - A_T_ndF.transpose() * Solve_A_T_T(b_T);
	}

	Vector CondensedRHS_noNeumannZeroDirichlet(const Vector& b_T)
	{
		assert(b_T.rows() == _nbh.Elements.size() * _diffPb.HHO->nCellUnknowns);

		return -A_T_ndF.transpose() * Solve_A_T_T(b_T);
	}



	Vector SolveCellUnknowns(const Vector& faceUnknowns, const Vector& b_T) const
	{
		int nCellUnknowns = _diffPb.HHO->nCellUnknowns;

		assert(faceUnknowns.rows() == _nbh.InteriorFaces.size() * _diffPb.HHO->nFaceUnknowns);
		assert(b_T.rows() == _nbh.Elements.size() * nCellUnknowns);

		return Solve_A_T_T(b_T - A_T_ndF * faceUnknowns);
	}

	Vector Solve_A_T_T(const Vector& v) const
	{
		int nCellUnknowns = _diffPb.HHO->nCellUnknowns;

		Vector result(_nbh.Elements.size() * nCellUnknowns);
		for (int i = 0; i < _nbh.Elements.size(); i++)
		{
			Diff_HHOElement<Dim>* e = _diffPb.HHOElement(_nbh.Elements[i]);
			result.segment(i * nCellUnknowns, nCellUnknowns) = e->AttSolver.solve(v.segment(i * nCellUnknowns, nCellUnknowns));
		}
		return result;
	}

	/*SparseMatrix A_T_T()
	{
		int nCellUnknowns = _diffPb.HHO->nCellUnknowns;

		SparseMatrix mat(_nbh.Elements.size() * nCellUnknowns, _nbh.Elements.size() * nCellUnknowns);
		NonZeroCoefficients coeffs;
		for (int i = 0; i < _nbh.Elements.size(); i++)
		{
			Diff_HHOElement<Dim>* e = _diffPb.HHOElement(_nbh.Elements[i]);
			coeffs.Add(i * nCellUnknowns, i * nCellUnknowns, e->A.topLeftCorner(nCellUnknowns, nCellUnknowns));
		}
		coeffs.Fill(mat);
		return mat;
	}*/

	// Reconstruction

	Vector ReconstructHigherOrder(const Vector& faceUnknowns, const Vector& dirichletCoeffs, const Vector& b_T)
	{
		int nFaceUnknowns = _diffPb.HHO->nFaceUnknowns;
		int nCellUnknowns = _diffPb.HHO->nCellUnknowns;
		int nReconstructUnknowns = _diffPb.HHO->nReconstructUnknowns;

		Vector v = b_T - A_T_ndF * faceUnknowns;

		Vector faceCoeffs((_nbh.InteriorFaces.size() + _nbh.BoundaryFaces.size()) * nFaceUnknowns);
		faceCoeffs.head(_nbh.InteriorFaces.size() * nFaceUnknowns) = faceUnknowns;
		faceCoeffs.tail(_nbh.BoundaryFaces.size() * nFaceUnknowns) = dirichletCoeffs;

		Vector reconstruction(_nbh.Elements.size() * nReconstructUnknowns);

		for (int i = 0; i < _nbh.Elements.size(); i++)
		{
			Diff_HHOElement<Dim>* e = _diffPb.HHOElement(_nbh.Elements[i]);

			Vector localHybrid(nCellUnknowns + nFaceUnknowns * e->Faces.size());
			localHybrid.head(nCellUnknowns) = e->AttSolver.solve(v.segment(i * nCellUnknowns, nCellUnknowns));
			for (auto f : e->Faces)
				localHybrid.segment(e->FirstDOFNumber(f), nFaceUnknowns) = faceCoeffs.segment(_nbh.FaceNumber(f->MeshFace) * nFaceUnknowns, nFaceUnknowns);

			reconstruction.segment(i * nReconstructUnknowns, nReconstructUnknowns) = e->Reconstruct(localHybrid);

		}
		return reconstruction;
	}

	/*Vector SolveFaceMassMatrixOnBoundary(const Vector& v)
	{
		int nFaceUnknowns = _diffPb.HHO->nFaceUnknowns;
		Vector res(v.rows());
		for (int i = 0; i < _nbh.BoundaryFaces.size(); i++)
		{
			Face<Dim>* f = _nbh.BoundaryFaces[i];
			Diff_HHOFace<Dim>* face = _diffPb.HHOFace(f);

			res.segment(i * nFaceUnknowns, nFaceUnknowns) = face->SolveMassMatrix(v.segment(i * nFaceUnknowns, nFaceUnknowns));
		}
		return res;
	}*/
};