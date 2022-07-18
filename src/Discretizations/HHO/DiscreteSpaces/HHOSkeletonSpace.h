#pragma once
#include "HHOFaceListSpace.h"


template<int Dim>
class HHOSkeletonSpace : public HHOFaceListSpace<Dim>
{
protected:
	const Mesh<Dim>* _mesh = nullptr;
	vector<Diff_HHOFace<Dim>>* _hhoFaces = nullptr;

public:
	HHOSkeletonSpace()
	{}

	HHOSkeletonSpace(const Mesh<Dim>* mesh, HHOParameters<Dim>* hho, vector<Diff_HHOFace<Dim>>& hhoFaces)
		: HHOFaceListSpace<Dim>(hho)
	{
		_mesh = mesh;
		_hhoFaces = &hhoFaces;
	}

protected:
	const vector<Face<Dim>*>& ListFaces() override
	{
		return _mesh->Faces;
	}

	Diff_HHOFace<Dim>* HHOFace(Face<Dim>* f) override
	{
		return &(*_hhoFaces)[f->Number];
	}

	BigNumber Number(Face<Dim>* f) override
	{
		return f->Number;
	}

public:
	double Measure() override
	{
		return _mesh->SkeletonMeasure();
	}
};
