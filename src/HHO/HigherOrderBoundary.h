#pragma once
#include "Diffusion_HHO.h"


template <int Dim>
class HigherOrderBoundary
{
private:
	vector<Diff_HHOFace<Dim>> _hhoFaces;
	Diffusion_HHO<Dim>* _diffPb;
	SparseMatrix _trace;
	SparseMatrix _normalDerivative;
public:
	HHOParameters<Dim>* HHO;
	Mesh<Dim>* _mesh;
	HHOBoundarySpace<Dim> BoundarySpace;

	HigherOrderBoundary() {}

	HigherOrderBoundary(Diffusion_HHO<Dim>* diffPb)
	{
		_diffPb = diffPb;
		FunctionalBasis<Dim - 1>* higherOrderFaceBasis = new FunctionalBasis<Dim - 1>(diffPb->HHO->FaceBasis->CreateSameBasisForDegree(diffPb->HHO->FaceBasis->GetDegree()+1));
		this->HHO = new HHOParameters<Dim>(diffPb->_mesh, "", diffPb->HHO->ReconstructionBasis, diffPb->HHO->CellBasis, higherOrderFaceBasis, diffPb->HHO->OrthogonalizeElemBasesCode, diffPb->HHO->OrthogonalizeFaceBasesCode);
		_mesh = diffPb->_mesh;
	}

	void Setup(bool setupTraceMatrix, bool setupNormalDerivativeMatrix)
	{
		_hhoFaces = vector<Diff_HHOFace<Dim>>(this->_mesh->BoundaryFaces.size());

		Diffusion_HHO<Dim>::InitReferenceShapes(HHO, nullptr);

		ParallelLoop<Face<Dim>*>::Execute(this->_mesh->BoundaryFaces, [this](Face<Dim>* f)
			{
				int i = f->Number - HHO->nInteriorFaces;
				_hhoFaces[i].MeshFace = f;
				_hhoFaces[i].InitHHO(HHO);
			}
		);

		BoundarySpace = HHOBoundarySpace(_mesh, HHO, _hhoFaces);

		_mesh->AssignNumberToBoundaryElements();

		if (setupTraceMatrix)
			SetupTraceMatrix();
		if (setupNormalDerivativeMatrix)
			SetupNormalDerivativeMatrix();
	}

private:
	Diff_HHOFace<Dim>* HHOFace(Face<Dim>* f)
	{
		return &this->_hhoFaces[f->Number - HHO->nInteriorFaces];
	}

	void SetupTraceMatrix()
	{
		FaceParallelLoop<Dim> parallelLoop(_mesh->BoundaryFaces);
		parallelLoop.ReserveChunkCoeffsSize(HHO->nReconstructUnknowns * HHO->nFaceUnknowns);

		parallelLoop.Execute([this](Face<Dim>* f, ParallelChunk<CoeffsChunk>* chunk)
			{
				int i = f->Number - HHO->nInteriorFaces;
				Diff_HHOFace<Dim>* face = HHOFace(f);

				int j = _mesh->BoundaryElementNumber(f->Element1);
				Diff_HHOElement<Dim>* elem = _diffPb->HHOElement(f->Element1);
				
				chunk->Results.Coeffs.Add(i * HHO->nFaceUnknowns, j * HHO->nReconstructUnknowns, face->Trace(elem->MeshElement, elem->ReconstructionBasis));
			});

		_trace = SparseMatrix(HHO->nBoundaryFaces * HHO->nFaceUnknowns, _mesh->NBoundaryElements() * HHO->nReconstructUnknowns);
		parallelLoop.Fill(_trace);
	}

	void SetupNormalDerivativeMatrix()
	{
		FaceParallelLoop<Dim> parallelLoop(_mesh->BoundaryFaces);
		parallelLoop.ReserveChunkCoeffsSize(HHO->nReconstructUnknowns * HHO->nFaceUnknowns);

		parallelLoop.Execute([this](Face<Dim>* f, ParallelChunk<CoeffsChunk>* chunk)
			{
				int i = f->Number - HHO->nInteriorFaces;
				Diff_HHOFace<Dim>* face = HHOFace(f);

				int j = _mesh->BoundaryElementNumber(f->Element1);
				Diff_HHOElement<Dim>* elem = _diffPb->HHOElement(f->Element1);

				chunk->Results.Coeffs.Add(i * HHO->nFaceUnknowns, j * HHO->nReconstructUnknowns, face->NormalDerivative(elem->MeshElement, elem->ReconstructionBasis));
			});

		_normalDerivative = SparseMatrix(HHO->nBoundaryFaces * HHO->nFaceUnknowns, _mesh->NBoundaryElements() * HHO->nReconstructUnknowns);
		parallelLoop.Fill(_normalDerivative);
	}

public:
	Vector Trace(const Vector& v)
	{
		assert(v.rows() == _mesh->NBoundaryElements() * HHO->nReconstructUnknowns);
		return _trace * v;
	}

	Vector NormalDerivative(const Vector& v)
	{
		assert(v.rows() == _mesh->NBoundaryElements() * HHO->nReconstructUnknowns);
		return _normalDerivative * v;
	}

	Vector AssembleNeumannTerm(const Vector& neumannHigherOrderCoeffs)
	{
		assert(neumannHigherOrderCoeffs.rows() == HHO->nBoundaryFaces * HHO->nFaceUnknowns);
		assert(_mesh->NeumannFaces.size() == _mesh->BoundaryFaces.size());

		Vector b_ndF = Vector::Zero(_diffPb->HHO->nTotalFaceUnknowns);
		ParallelLoop<Face<Dim>*>::Execute(this->_mesh->NeumannFaces, [this, &b_ndF, &neumannHigherOrderCoeffs](Face<Dim>* f)
			{
				Diff_HHOFace<Dim>* loFace = _diffPb->HHOFace(f);
				int loUnknowns = _diffPb->HHO->nFaceUnknowns;
				BigNumber i = f->Number * loUnknowns;
				assert(i >= _diffPb->HHO->nInteriorFaces * loUnknowns);

				Diff_HHOFace<Dim>* hoFace = this->HHOFace(f);
				int hoUnknowns = this->HHO->nFaceUnknowns;
				BigNumber j = (f->Number - HHO->nInteriorFaces) * hoUnknowns;

				b_ndF.segment(i, loUnknowns) = loFace->MassMatrix(hoFace->Basis) * neumannHigherOrderCoeffs.segment(j, hoUnknowns);
			}
		);
		return b_ndF;
	}

	Vector AssembleDirichletTerm(const Vector& dirichletHigherOrderCoeffs)
	{
		assert(dirichletHigherOrderCoeffs.rows() == HHO->nDirichletCoeffs);

		Vector x_dF = Vector::Zero(_diffPb->HHO->nDirichletCoeffs);
		ParallelLoop<Face<Dim>*>::Execute(this->_mesh->DirichletFaces, [this, &x_dF, &dirichletHigherOrderCoeffs](Face<Dim>* f)
			{
				Diff_HHOFace<Dim>* loFace = _diffPb->HHOFace(f);
				int loUnknowns = _diffPb->HHO->nFaceUnknowns;
				BigNumber i = (f->Number - HHO->nInteriorFaces - HHO->nNeumannFaces) * loUnknowns;

				Diff_HHOFace<Dim>* hoFace = this->HHOFace(f);
				int hoUnknowns = this->HHO->nFaceUnknowns;
				BigNumber j = (f->Number - HHO->nInteriorFaces - HHO->nNeumannFaces) * hoUnknowns;

				x_dF.segment(i, loUnknowns) = loFace->SolveMassMatrix(loFace->MassMatrix(hoFace->Basis) * dirichletHigherOrderCoeffs.segment(j, hoUnknowns));
			}
		);
		return x_dF;
	}

	~HigherOrderBoundary()
	{
		//delete HHO->FaceBasis;
		//delete HHO;
	}
};
