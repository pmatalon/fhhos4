#pragma once
#include "HHOFaceListSpace.h"


template<int Dim>
class HHODirichletSpace : public HHOFaceListSpace<Dim>
{
private:
	const Mesh<Dim>* _mesh = nullptr;
	vector<Diff_HHOFace<Dim>>* _hhoFaces = nullptr;
	HHOParameters<Dim>* HHO = nullptr;

public:
	HHODirichletSpace() {}

	HHODirichletSpace(const Mesh<Dim>* mesh, HHOParameters<Dim>* hho, vector<Diff_HHOFace<Dim>>& hhoFaces)
		: HHOFaceListSpace<Dim>(hho->nFaceUnknowns)
	{
		_mesh = mesh;
		HHO = hho;
		_hhoFaces = &hhoFaces;
		assert(HHO->nDirichletFaces == _mesh->DirichletFaces.size());
		assert(HHO->nDirichletFaces * HHO->nFaceUnknowns == HHO->nDirichletCoeffs);
	}

private:
	const vector<Face<Dim>*>& ListFaces() override
	{
		return _mesh->DirichletFaces;
	}

	Diff_HHOFace<Dim>* HHOFace(Face<Dim>* f) override
	{
		return &(*_hhoFaces)[f->Number];
	}

	BigNumber Number(Face<Dim>* f) override
	{
		return f->Number - HHO->nInteriorAndNeumannFaces;
	}

public:
	double Measure() override
	{
		Utils::FatalError("HHODirichletSpace::Measure() not implemented");
	}
};
