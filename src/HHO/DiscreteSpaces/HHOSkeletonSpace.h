#pragma once
#include "IDiscreteSpace.h"
#include "../Diff_HHOElement.h"
#include "../../Utils/ElementParallelLoop.h"


template<int Dim>
class HHOSkeletonSpace : public IDiscreteSpace
{
private:
	Mesh<Dim>* _mesh = nullptr;
	vector<Diff_HHOFace<Dim>>* _hhoFaces = nullptr;
	HHOParameters<Dim>* HHO = nullptr;

public:
	HHOSkeletonSpace()
	{}

	HHOSkeletonSpace(Mesh<Dim>* mesh, HHOParameters<Dim>* hho, vector<Diff_HHOFace<Dim>>& hhoFaces)
	{
		HHO = hho;
		_mesh = mesh;
		_hhoFaces = &hhoFaces;
	}

	Diff_HHOFace<Dim>* HHOFace(Face<Dim>* f)
	{
		return &(*_hhoFaces)[f->Number];
	}

	//--------------------------------------------//
	// Implementation of interface IDiscreteSpace //
	//--------------------------------------------//

	Vector InnerProdWithBasis(DomFunction func) override
	{
		Vector innerProds = Vector(HHO->nTotalFaceUnknowns);
		ParallelLoop<Face<Dim>*>::Execute(this->_mesh->Faces, [this, &innerProds, func](Face<Dim>* f)
			{
				Diff_HHOFace<Dim>* face = HHOFace(f);
				BigNumber i = f->Number * HHO->nFaceUnknowns;

				innerProds.segment(i, HHO->nFaceUnknowns) = face->InnerProductWithBasis(func);
			}
		);
		return innerProds;
	}

	Vector SolveMassMatrix(const Vector& v) override
	{
		assert(v.rows() == HHO->nFaces * HHO->nFaceUnknowns);
		Vector res(v.rows());
		ParallelLoop<Face<Dim>*>::Execute(this->_mesh->Faces, [this, &v, &res](Face<Dim>* f)
			{
				Diff_HHOFace<Dim>* face = HHOFace(f);
				BigNumber i = f->Number * HHO->nFaceUnknowns;

				res.segment(i, HHO->nFaceUnknowns) = face->SolveMassMatrix(v.segment(i, HHO->nFaceUnknowns));
			});
		return res;
	}

	double Measure() override
	{
		return _mesh->SkeletonMeasure();
	}
};
