#pragma once
#include "Diffusion_HHO.h"


template <int Dim>
class HigherOrderBoundary
{
private:
	vector<Diff_HHOFace<Dim>> _hhoFaces;
	Diffusion_HHO<Dim>* _diffPb;
	SparseMatrix _trace;
public:
	HHOParameters<Dim>* HHO;
	Mesh<Dim>* _mesh;

	HigherOrderBoundary() {}

	HigherOrderBoundary(Diffusion_HHO<Dim>* diffPb)
	{
		_diffPb = diffPb;
		FunctionalBasis<Dim - 1>* higherOrderFaceBasis = new FunctionalBasis<Dim - 1>(diffPb->HHO->FaceBasis->CreateSameBasisForDegree(diffPb->HHO->FaceBasis->GetDegree()+1));
		this->HHO = new HHOParameters<Dim>(diffPb->_mesh, "", diffPb->HHO->ReconstructionBasis, diffPb->HHO->CellBasis, higherOrderFaceBasis, diffPb->HHO->OrthogonalizeElemBasesCode, diffPb->HHO->OrthogonalizeFaceBasesCode);
		_mesh = diffPb->_mesh;
	}

	void Setup()
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

		_mesh->AssignNumberToBoundaryElements();
		SetupTraceMatrix();
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

public:
	Vector Trace(const Vector& v)
	{
		assert(v.rows() == _mesh->NBoundaryElements() * HHO->nReconstructUnknowns);
		return _trace * v;
	}

	Vector InnerProdWithBoundaryFaceBasis(DomFunction func)
	{
		Vector innerProds = Vector(HHO->nBoundaryFaces * HHO->nFaceUnknowns);
		ParallelLoop<Face<Dim>*>::Execute(this->_mesh->BoundaryFaces, [this, &innerProds, func](Face<Dim>* f)
			{
				Diff_HHOFace<Dim>* face = HHOFace(f);
				BigNumber i = (f->Number - HHO->nInteriorFaces) * HHO->nFaceUnknowns;

				innerProds.segment(i, HHO->nFaceUnknowns) = face->InnerProductWithBasis(func);
			}
		);
		return innerProds;
	}

	Vector ProjectOnBoundaryDiscreteSpace(DomFunction func)
	{
		Vector vectorOfDoFs = Vector(HHO->nBoundaryFaces * HHO->nFaceUnknowns);
		ParallelLoop<Face<Dim>*>::Execute(this->_mesh->BoundaryFaces, [this, &vectorOfDoFs, func](Face<Dim>* f)
			{
				Diff_HHOFace<Dim>* face = HHOFace(f);
				BigNumber i = (face->Number() - HHO->nInteriorFaces) * HHO->nFaceUnknowns;

				vectorOfDoFs.segment(i, HHO->nFaceUnknowns) = face->ProjectOnBasis(func);
			}
		);
		return vectorOfDoFs;
	}

	Vector SolveBoundaryFaceMassMatrix(const Vector& v)
	{
		assert(v.rows() == HHO->nBoundaryFaces * HHO->nFaceUnknowns);
		Vector res(v.rows());
		ParallelLoop<Face<Dim>*>::Execute(this->_mesh->BoundaryFaces, [this, &v, &res](Face<Dim>* f)
			{
				Diff_HHOFace<Dim>* face = HHOFace(f);
				BigNumber i = (f->Number - HHO->nInteriorFaces) * HHO->nFaceUnknowns;

				res.segment(i, HHO->nFaceUnknowns) = face->SolveMassMatrix(v.segment(i, HHO->nFaceUnknowns));
			});
		return res;
	}

	double IntegralOverBoundaryFromFaceCoeffs(const Vector& boundaryFaceCoeffs)
	{
		assert(boundaryFaceCoeffs.rows() == HHO->nBoundaryFaces * HHO->nFaceUnknowns);

		struct ChunkResult { double total = 0; };

		ParallelLoop<Face<Dim>*, ChunkResult> parallelLoop(_mesh->BoundaryFaces);
		parallelLoop.Execute([this, &boundaryFaceCoeffs](Face<Dim>* f, ParallelChunk<ChunkResult>* chunk)
			{
				Diff_HHOFace<Dim>* hhoFace = HHOFace(f);
				auto i = (f->Number - HHO->nInteriorFaces) * HHO->nFaceUnknowns;
				chunk->Results.total += hhoFace->Integral(boundaryFaceCoeffs.segment(i, HHO->nFaceUnknowns));
			});

		double total = 0;
		parallelLoop.AggregateChunkResults([&total](ChunkResult chunkResult)
			{
				total += chunkResult.total;
			});
		return total;
	}

	double L2InnerProdOnBoundary(const Vector& v1, const Vector& v2)
	{
		assert(v1.rows() == HHO->nBoundaryFaces * HHO->nFaceUnknowns);
		assert(v2.rows() == HHO->nBoundaryFaces * HHO->nFaceUnknowns);

		struct ChunkResult { double total = 0; };

		ParallelLoop<Face<Dim>*, ChunkResult> parallelLoop(_mesh->BoundaryFaces);
		parallelLoop.Execute([this, &v1, &v2](Face<Dim>* f, ParallelChunk<ChunkResult>* chunk)
			{
				Diff_HHOFace<Dim>* face = HHOFace(f);
				BigNumber i = (face->Number() - HHO->nInteriorFaces) * HHO->nFaceUnknowns;

				chunk->Results.total += face->InnerProd(v1.segment(i, HHO->nFaceUnknowns), v2.segment(i, HHO->nFaceUnknowns));
			});

		double total = 0;
		parallelLoop.AggregateChunkResults([&total](ChunkResult chunkResult)
			{
				total += chunkResult.total;
			});
		return total;
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

	~HigherOrderBoundary()
	{
		//delete HHO->FaceBasis;
		//delete HHO;
	}
};
