#pragma once
#include "Diffusion_HHO.h"


template <int Dim>
class HigherOrderBoundary
{
private:
	vector<Diff_HHOFace<Dim>> _hhoFaces;
	Diffusion_HHO<Dim>* _diffPb;
	SparseMatrix _trace;
	SparseMatrix _boundaryFaceMass;
	//SparseMatrix _normalDerivative;
public:
	HHOParameters<Dim>* HHO;
	Mesh<Dim>* _mesh;
	HHOBoundarySpace<Dim> BoundarySpace;

	HigherOrderBoundary() {}

	HigherOrderBoundary(Diffusion_HHO<Dim>* diffPb)
	{
		_diffPb = diffPb;
		FunctionalBasis<Dim - 1>* higherOrderFaceBasis = diffPb->HHO->FaceBasis->CreateSameBasisForDegree(diffPb->HHO->FaceBasis->GetDegree()+1);
		this->HHO = new HHOParameters<Dim>(diffPb->_mesh, "", diffPb->HHO->ReconstructionBasis, diffPb->HHO->CellBasis, higherOrderFaceBasis, diffPb->HHO->OrthogonalizeElemBasesCode, diffPb->HHO->OrthogonalizeFaceBasesCode);
		_mesh = diffPb->_mesh;
	}

	void Setup()
	{
		_hhoFaces = vector<Diff_HHOFace<Dim>>(this->_mesh->BoundaryFaces.size());

		Diffusion_HHO<Dim>::InitReferenceShapes(HHO, nullptr);

		#pragma omp parallel for
		for (Face<Dim>* f : this->_mesh->BoundaryFaces)
		{
			int i = f->Number - HHO->nInteriorFaces;
			_hhoFaces[i].MeshFace = f;
			_hhoFaces[i].InitHHO(HHO);
		}

		BoundarySpace = HHOBoundarySpace(_mesh, HHO, _hhoFaces);

		SetupTraceMatrix();
		SetupBoundaryFaceMassMatrix();
	}

	Diff_HHOFace<Dim>* HHOFace(Face<Dim>* f)
	{
		return &this->_hhoFaces[f->Number - HHO->nInteriorFaces];
	}

private:
	void SetupTraceMatrix()
	{
		ThreadLocalCoeffs coeffs(_mesh->BoundaryFaces.size(), HHO->nReconstructUnknowns * HHO->nFaceUnknowns);
		#pragma omp parallel for
		for (Face<Dim>* f : _mesh->BoundaryFaces)
		{
			Diff_HHOFace<Dim>* face = HHOFace(f);
			int i = f->Number - HHO->nInteriorFaces;

			Diff_HHOElement<Dim>* elem = _diffPb->HHOElement(f->Element1);
			int j = _mesh->BoundaryElementNumber(f->Element1);
			
			coeffs.Local().Add(i * HHO->nFaceUnknowns, j * HHO->nReconstructUnknowns, face->Trace(elem->MeshElement, elem->ReconstructionBasis));
		}

		_trace = SparseMatrix(HHO->nBoundaryFaces * HHO->nFaceUnknowns, _mesh->NBoundaryElements() * HHO->nReconstructUnknowns);
		coeffs.Fill(_trace);
	}

	void SetupBoundaryFaceMassMatrix()
	{
		ThreadLocalCoeffs coeffs(_mesh->BoundaryFaces.size(), HHO->nFaceUnknowns * HHO->nFaceUnknowns);
		#pragma omp parallel for
		for (Face<Dim>* f : _mesh->BoundaryFaces)
		{
			Diff_HHOFace<Dim>* face = HHOFace(f);
			int i = f->Number - HHO->nInteriorFaces;

			coeffs.Local().Add(i * HHO->nFaceUnknowns, i * HHO->nFaceUnknowns, face->MassMatrix());
		}

		_boundaryFaceMass = SparseMatrix(HHO->nBoundaryFaces * HHO->nFaceUnknowns, HHO->nBoundaryFaces * HHO->nFaceUnknowns);
		coeffs.Fill(_boundaryFaceMass);
	}

	/*void SetupNormalDerivativeMatrix()
	{
		ThreadLocalCoeffs coeffs(_mesh->BoundaryFaces.size(), HHO->nReconstructUnknowns * HHO->nFaceUnknowns);
		#pragma omp parallel for
		for (Face<Dim>* f : _mesh->BoundaryFaces)
		{
			int i = f->Number - HHO->nInteriorFaces;
			Diff_HHOFace<Dim>* face = HHOFace(f);

			int j = _mesh->BoundaryElementNumber(f->Element1);
			Diff_HHOElement<Dim>* elem = _diffPb->HHOElement(f->Element1);

			coeffs.Local().Add(i * HHO->nFaceUnknowns, j * HHO->nReconstructUnknowns, face->NormalDerivative(elem->MeshElement, elem->ReconstructionBasis));
		}

		_normalDerivative = SparseMatrix(HHO->nBoundaryFaces * HHO->nFaceUnknowns, _mesh->NBoundaryElements() * HHO->nReconstructUnknowns);
		coeffs.Fill(_normalDerivative);
	}*/

public:
	Vector Trace(const Vector& v)
	{
		assert(v.rows() == _mesh->NBoundaryElements() * HHO->nReconstructUnknowns);
		return _trace * v;
	}

	SparseMatrix& TraceMatrix()
	{
		return _trace;
	}

	SparseMatrix& BoundaryFaceMassMatrix()
	{
		return _boundaryFaceMass;
	}

	/*Vector NormalDerivative(const Vector& v)
	{
		assert(v.rows() == _mesh->NBoundaryElements() * HHO->nReconstructUnknowns);
		return _normalDerivative * v;
	}*/

	Vector AssembleNeumannTerm(const Vector& neumannHigherOrderCoeffs)
	{
		assert(neumannHigherOrderCoeffs.rows() == HHO->nBoundaryFaces * HHO->nFaceUnknowns);
		assert(_mesh->NeumannFaces.size() == _mesh->BoundaryFaces.size());

		Vector b_ndF = Vector::Zero(_diffPb->HHO->nTotalFaceUnknowns);
		#pragma omp parallel for
		for (Face<Dim>* f : this->_mesh->NeumannFaces)
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
		return b_ndF;
	}

	Vector AssembleDirichletTerm(const Vector& dirichletHigherOrderCoeffs)
	{
		assert(dirichletHigherOrderCoeffs.rows() == HHO->nDirichletCoeffs);
		assert(_mesh->DirichletFaces.size() == _mesh->BoundaryFaces.size());

		Vector x_dF = Vector::Zero(_diffPb->HHO->nDirichletCoeffs);
		#pragma omp parallel for
		for (Face<Dim>* f : this->_mesh->DirichletFaces)
		{
			Diff_HHOFace<Dim>* loFace = _diffPb->HHOFace(f);
			int loUnknowns = _diffPb->HHO->nFaceUnknowns;

			Diff_HHOFace<Dim>* hoFace = this->HHOFace(f);
			int hoUnknowns = this->HHO->nFaceUnknowns;

			BigNumber i = f->Number - HHO->nInteriorFaces - HHO->nNeumannFaces;

			x_dF.segment(i * loUnknowns, loUnknowns) = loFace->ProjectOnBasis(hoFace->Basis) * dirichletHigherOrderCoeffs.segment(i * hoUnknowns, hoUnknowns);
		}
		return x_dF;
	}

	Vector ExtractBoundaryElements(const Vector& v)
	{
		assert(v.rows() == _mesh->Elements.size() * HHO->nReconstructUnknowns);

		Vector boundary(_mesh->NBoundaryElements() * HHO->nReconstructUnknowns);
		#pragma omp parallel for
		for (Element<Dim>* e : _mesh->Elements)
		{
			if (e->IsOnBoundary())
			{
				int boundaryElemNumber = _mesh->BoundaryElementNumber(e);
				boundary.segment(boundaryElemNumber * HHO->nReconstructUnknowns, HHO->nReconstructUnknowns) = v.segment(e->Number * HHO->nReconstructUnknowns, HHO->nReconstructUnknowns);
			}
		}
		return boundary;
	}

	~HigherOrderBoundary()
	{
		//delete HHO->FaceBasis;
		//delete HHO;
	}
};
