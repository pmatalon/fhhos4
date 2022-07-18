#pragma once
#include "HHOFaceListSpace.h"


template<int Dim>
class HHONeumannSpace : public HHOFaceListSpace<Dim>
{
private:
	const Mesh<Dim>* _mesh = nullptr;
	vector<Diff_HHOFace<Dim>>* _hhoFaces = nullptr;

public:
	HHONeumannSpace() {}

	HHONeumannSpace(const Mesh<Dim>* mesh, HHOParameters<Dim>* hho, vector<Diff_HHOFace<Dim>>& hhoFaces)
		: HHOFaceListSpace<Dim>(hho)
	{
		_mesh = mesh;
		_hhoFaces = &hhoFaces;
		assert(hho->nNeumannFaces == _mesh->NeumannFaces.size());
		assert(hho->nNeumannFaces * hho->nFaceUnknowns == hho->nNeumannUnknowns);
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
		return f->Number - this->HHO->nInteriorFaces;
	}

public:
	double Measure() override
	{
		Utils::FatalError("HHONeumannSpace::Measure() not implemented");
	}
};
