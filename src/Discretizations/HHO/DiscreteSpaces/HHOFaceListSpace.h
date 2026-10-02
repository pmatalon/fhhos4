#pragma once
#include "IDiscreteSpace.h"
#include "../Diff_HHOElement.h"
#include "../../../Utils/Parallelism.h"


template<int Dim>
class HHOFaceListSpace : public IDiscreteSpace
{
protected:
	int _nUnknowns;
	HHOParameters<Dim>* HHO = nullptr;

	HHOFaceListSpace() {}

	HHOFaceListSpace(HHOParameters<Dim>* hho)
	{
		_nUnknowns = hho->nFaceUnknowns;
		HHO = hho;
	}
private:
	virtual const vector<Face<Dim>*>& ListFaces() = 0;

	virtual Diff_HHOFace<Dim>* HHOFace(Face<Dim>* f) = 0;

	virtual BigNumber Number(Face<Dim>* f) = 0;

public:
	//--------------------------------------------//
	// Implementation of interface IDiscreteSpace //
	//--------------------------------------------//

	BigNumber Dimension() override
	{
		return ListFaces().size() * _nUnknowns;
	}

	Vector InnerProdWithBasis(DomFunction func) override
	{
		Vector innerProds = Vector(Dimension());
		#pragma omp parallel for
		for (Face<Dim>* f : ListFaces())
		{
			BigNumber i = Number(f) * _nUnknowns;
			innerProds.segment(i, _nUnknowns) = HHOFace(f)->InnerProductWithBasis(func);
		}
		return innerProds;
	}

	Vector ApplyMassMatrix(const Vector& v) override
	{
		assert(v.rows() == Dimension());

		if (HHO->OrthonormalizeFaceBases())
			return v;

		Vector res(v.rows());
		#pragma omp parallel for
		for (Face<Dim>* f : ListFaces())
		{
			BigNumber i = Number(f) * _nUnknowns;
			res.segment(i, _nUnknowns) = HHOFace(f)->ApplyMassMatrix(v.segment(i, _nUnknowns));
		}
		return res;
	}

	Vector SolveMassMatrix(const Vector& v) override
	{
		assert(v.rows() == Dimension());

		if (HHO->OrthonormalizeFaceBases())
			return v;

		Vector res(v.rows());
		#pragma omp parallel for
		for (Face<Dim>* f : ListFaces())
		{
			BigNumber i = Number(f) * _nUnknowns;
			res.segment(i, _nUnknowns) = HHOFace(f)->SolveMassMatrix(v.segment(i, _nUnknowns));
		}
		return res;
	}

	Vector Project(DomFunction func) override
	{
		Vector vectorOfDoFs = Vector(Dimension());
		#pragma omp parallel for
		for (Face<Dim>* f : ListFaces())
		{
			BigNumber i = Number(f) * _nUnknowns;
			vectorOfDoFs.segment(i, _nUnknowns) = HHOFace(f)->ProjectOnBasis(func);
		}
		return vectorOfDoFs;
	}

	double L2InnerProd(const Vector& v1, const Vector& v2) override
	{
		assert(v1.rows() == Dimension());
		assert(v2.rows() == Dimension());

		if (HHO->OrthonormalizeFaceBases())
			return v1.dot(v2);

		double total = 0;
		#pragma omp parallel for reduction(+:total)
		for (Face<Dim>* f : ListFaces())
		{
			BigNumber i = Number(f) * _nUnknowns;
			total += HHOFace(f)->InnerProd(v1.segment(i, _nUnknowns), v2.segment(i, _nUnknowns));
		}
		return total;
	}

	double Integral(const Vector& v) override
	{
		assert(v.rows() == Dimension());

		double total = 0;
		#pragma omp parallel for reduction(+:total)
		for (Face<Dim>* f : ListFaces())
		{
			auto i = Number(f) * _nUnknowns;
			total += HHOFace(f)->Integral(v.segment(i, _nUnknowns));
		}
		return total;
	}

	double Integral(DomFunction func) override
	{
		double total = 0;
		#pragma omp parallel for reduction(+:total)
		for (Face<Dim>* f : ListFaces())
			total += f->Integral(func);
		return total;
	}
};