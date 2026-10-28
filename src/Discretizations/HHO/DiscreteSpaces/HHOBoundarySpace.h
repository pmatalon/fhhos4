#pragma once
#include "HHOFaceListSpace.h"
#include "../../../Utils/Communicator.h"


template<int Dim>
class HHOBoundarySpace : public HHOFaceListSpace<Dim>
{
private:
	const Mesh<Dim>* _mesh = nullptr;
	vector<Diff_HHOFace<Dim>>* _hhoFaces = nullptr;
	bool _onlyBoundaryHHOFacesProvided = false;
	const vector<Face<Dim>*>* _ownFaces = nullptr; // those of Integral(func), reduced over the processes

public:
	HHOBoundarySpace()
	{}

	HHOBoundarySpace(const Mesh<Dim>* mesh, HHOParameters<Dim>* hho, vector<Diff_HHOFace<Dim>>& hhoFaces)
		: HHOBoundarySpace(mesh, hho, hhoFaces, mesh->BoundaryFaces)
	{}

	HHOBoundarySpace(const Mesh<Dim>* mesh, HHOParameters<Dim>* hho, vector<Diff_HHOFace<Dim>>& hhoFaces, const vector<Face<Dim>*>& ownFaces)
		: HHOFaceListSpace<Dim>(hho)
	{
		_mesh = mesh;
		_hhoFaces = &hhoFaces;
		_onlyBoundaryHHOFacesProvided = hhoFaces.size() == hho->nBoundaryFaces;
		_ownFaces = &ownFaces;
	}

	double Integral(DomFunction func) override
	{
		double total = 0;
		#pragma omp parallel for reduction(+:total)
		for (Face<Dim>* f : *_ownFaces)
			total += f->Integral(func);
		return Communicator::Current().Sum(total);
	}
	using HHOFaceListSpace<Dim>::Integral;

private:
	const vector<Face<Dim>*>& ListFaces() override
	{
		return _mesh->BoundaryFaces;
	}

	Diff_HHOFace<Dim>* HHOFace(Face<Dim>* f) override
	{
		if (_onlyBoundaryHHOFacesProvided)
			return &(*_hhoFaces)[f->Number - this->HHO->nInteriorFaces];
		else
			return &(*_hhoFaces)[f->Number];
	}

	BigNumber Number(Face<Dim>* f) override
	{
		return f->Number - this->HHO->nInteriorFaces;
	}

public:
	double Measure() override
	{
		return _mesh->BoundaryMeasure();
	}
};
