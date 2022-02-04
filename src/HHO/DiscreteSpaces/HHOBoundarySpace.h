#pragma once
#include "IDiscreteSpace.h"
#include "../Diff_HHOElement.h"
#include "../../Utils/ElementParallelLoop.h"


template<int Dim>
class HHOBoundarySpace : public IDiscreteSpace
{
private:
	Mesh<Dim>* _mesh = nullptr;
	vector<Diff_HHOFace<Dim>>* _hhoFaces = nullptr;
	HHOParameters<Dim>* HHO = nullptr;
	bool _onlyBoundaryHHOFacesProvided = false;

public:
	HHOBoundarySpace()
	{}

	HHOBoundarySpace(Mesh<Dim>* mesh, HHOParameters<Dim>* hho, vector<Diff_HHOFace<Dim>>& hhoFaces)
	{
		HHO = hho;
		_mesh = mesh;
		_hhoFaces = &hhoFaces;
		_onlyBoundaryHHOFacesProvided = hhoFaces.size() == hho->nBoundaryFaces;
	}

	Diff_HHOFace<Dim>* HHOFace(Face<Dim>* f)
	{
		if (_onlyBoundaryHHOFacesProvided)
			return &(*_hhoFaces)[f->Number - HHO->nInteriorFaces];
		else
			return &(*_hhoFaces)[f->Number];
	}

	//--------------------------------------------//
	// Implementation of interface IDiscreteSpace //
	//--------------------------------------------//

	Vector InnerProdWithBasis(DomFunction func) override
	{
		Vector innerProds = Vector(HHO->nBoundaryFaces * HHO->nFaceUnknowns);
		ParallelLoop<Face<Dim>*>::Execute(this->_mesh->BoundaryFaces, [this, &innerProds, func](Face<Dim>* f)
			{
				Diff_HHOFace<Dim>* face = HHOFace(f);
				BigNumber i = (face->Number() - HHO->nInteriorFaces) * HHO->nFaceUnknowns;

				innerProds.segment(i, HHO->nFaceUnknowns) = face->InnerProductWithBasis(func);
			}
		);
		return innerProds;
	}

	Vector SolveMassMatrix(const Vector& v) override
	{
		assert(v.rows() == HHO->nBoundaryFaces * HHO->nFaceUnknowns);
		Vector res(v.rows());
		ParallelLoop<Face<Dim>*>::Execute(this->_mesh->BoundaryFaces, [this, &v, &res](Face<Dim>* f)
			{
				Diff_HHOFace<Dim>* face = HHOFace(f);
				BigNumber i = (face->Number() - HHO->nInteriorFaces) * HHO->nFaceUnknowns;

				res.segment(i, HHO->nFaceUnknowns) = face->SolveMassMatrix(v.segment(i, HHO->nFaceUnknowns));
			});
		return res;
	}

	double Measure() override
	{
		return _mesh->BoundaryMeasure();
	}
};
