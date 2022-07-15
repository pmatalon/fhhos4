#pragma once
#include "HHOFaceListSpace.h"


template<int Dim>
class HHOBoundarySpace : public HHOFaceListSpace<Dim>
{
private:
	const Mesh<Dim>* _mesh = nullptr;
	vector<Diff_HHOFace<Dim>>* _hhoFaces = nullptr;
	HHOParameters<Dim>* HHO = nullptr;
	bool _onlyBoundaryHHOFacesProvided = false;

public:
	HHOBoundarySpace()
	{}

	HHOBoundarySpace(const Mesh<Dim>* mesh, HHOParameters<Dim>* hho, vector<Diff_HHOFace<Dim>>& hhoFaces)
		: HHOFaceListSpace<Dim>(hho->nFaceUnknowns)
	{
		_mesh = mesh;
		HHO = hho;
		_hhoFaces = &hhoFaces;
		_onlyBoundaryHHOFacesProvided = hhoFaces.size() == hho->nBoundaryFaces;
	}

private:
	const vector<Face<Dim>*>& ListFaces() override
	{
		return _mesh->BoundaryFaces;
	}

	Diff_HHOFace<Dim>* HHOFace(Face<Dim>* f) override
	{
		if (_onlyBoundaryHHOFacesProvided)
			return &(*_hhoFaces)[f->Number - HHO->nInteriorFaces];
		else
			return &(*_hhoFaces)[f->Number];
	}

	BigNumber Number(Face<Dim>* f) override
	{
		return f->Number - HHO->nInteriorFaces;
	}

public:
	double Measure() override
	{
		return _mesh->BoundaryMeasure();
	}
};
