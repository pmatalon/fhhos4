#pragma once
#include "HHOFaceListSpace.h"


template<int Dim>
class HHONeumannSpace : public HHOFaceListSpace<Dim>
{
private:
	const Mesh<Dim>* _mesh = nullptr;
	vector<Diff_HHOFace<Dim>>* _hhoFaces = nullptr;
	HHOParameters<Dim>* HHO = nullptr;

public:
	HHONeumannSpace() {}

	HHONeumannSpace(const Mesh<Dim>* mesh, HHOParameters<Dim>* hho, vector<Diff_HHOFace<Dim>>& hhoFaces)
		: HHOFaceListSpace<Dim>(hho->nFaceUnknowns)
	{
		_mesh = mesh;
		HHO = hho;
		_hhoFaces = &hhoFaces;
		assert(HHO->nNeumannFaces == _mesh->NeumannFaces.size());
		assert(HHO->nNeumannFaces * HHO->nFaceUnknowns == HHO->nNeumannUnknowns);
	}

private:
	const vector<Face<Dim>*>& ListFaces() override
	{
		return _mesh->NeumannFaces;
	}

	Diff_HHOFace<Dim>* HHOFace(Face<Dim>* f) override
	{
		return &(*_hhoFaces)[f->Number];
	}

	BigNumber Number(Face<Dim>* f) override
	{
		return f->Number - HHO->nInteriorFaces;
	}

public:
	double Measure() override
	{
		Utils::FatalError("HHONeumannSpace::Measure() not implemented");
	}
};
