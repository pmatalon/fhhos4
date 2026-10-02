#pragma once
#include "IDiscreteSpace.h"
#include "../Diff_HHOElement.h"
#include "../../../Utils/Parallelism.h"


template<int Dim>
class HHOReconstructSpace : public IDiscreteSpace
{
private:
	const Mesh<Dim>* _mesh = nullptr;
	vector<Diff_HHOElement<Dim>>* _hhoElements = nullptr;
	HHOParameters<Dim>* HHO = nullptr;

public:
	HHOReconstructSpace()
	{}

	HHOReconstructSpace(const Mesh<Dim>* mesh, HHOParameters<Dim>* hho, vector<Diff_HHOElement<Dim>>& hhoElements)
	{
		HHO = hho;
		_mesh = mesh;
		_hhoElements = &hhoElements;
	}

	Diff_HHOElement<Dim>* HHOElement(Element<Dim>* e)
	{
		return &(*_hhoElements)[e->Number];
	}

	//--------------------------------------------//
	// Implementation of interface IDiscreteSpace //
	//--------------------------------------------//
	
	BigNumber Dimension() override
	{
		return _mesh->Elements.size() * HHO->nReconstructUnknowns;
	}

	double Measure() override
	{
		return _mesh->Measure();
	}

	Vector InnerProdWithBasis(DomFunction func) override
	{
		Vector innerProds = Vector(Dimension());
		#pragma omp parallel for
		for (Element<Dim>* e : this->_mesh->Elements)
		{
			Diff_HHOElement<Dim>* elem = HHOElement(e);
			BigNumber i = e->Number * HHO->nReconstructUnknowns;

			innerProds.segment(i, HHO->nReconstructUnknowns) = elem->InnerProductWithBasis(elem->ReconstructionBasis, func);
		}
		return innerProds;
	}

	Vector ApplyMassMatrix(const Vector& v) override
	{
		assert(v.rows() == Dimension());

		if (HHO->OrthonormalizeElemBases())
			return v;

		Vector res(v.rows());
		#pragma omp parallel for
		for (Element<Dim>* e : this->_mesh->Elements)
		{
			BigNumber i = e->Number * HHO->nReconstructUnknowns;
			res.segment(i, HHO->nReconstructUnknowns) = HHOElement(e)->ApplyReconstructMassMatrix(v.segment(i, HHO->nReconstructUnknowns));
		}
		return res;
	}

	Vector SolveMassMatrix(const Vector& v) override
	{
		assert(v.rows() == Dimension());

		if (HHO->OrthonormalizeElemBases())
			return v;

		Vector res(v.rows());
		#pragma omp parallel for
		for (Element<Dim>* e : this->_mesh->Elements)
		{
			BigNumber i = e->Number * HHO->nReconstructUnknowns;
			res.segment(i, HHO->nReconstructUnknowns) = HHOElement(e)->SolveReconstructMassMatrix(v.segment(i, HHO->nReconstructUnknowns));
		}
		return res;
	}
	
	Vector Project(DomFunction func) override
	{
		Vector vectorOfDoFs = Vector(Dimension());
		#pragma omp parallel for
		for (Element<Dim>* e : this->_mesh->Elements)
		{
			BigNumber i = e->Number * HHO->nReconstructUnknowns;
			vectorOfDoFs.segment(i, HHO->nReconstructUnknowns) = HHOElement(e)->ProjectOnReconstructBasis(func);
		}
		return vectorOfDoFs;
	}

	double L2InnerProd(const Vector& v1, const Vector& v2) override
	{
		assert(v1.rows() == Dimension());
		assert(v2.rows() == Dimension());

		if (HHO->OrthonormalizeElemBases())
			return v1.dot(v2);

		double total = 0;
		#pragma omp parallel for reduction(+:total)
		for (Element<Dim>* e : _mesh->Elements)
		{
			BigNumber i = e->Number * HHO->nReconstructUnknowns;
			total += v1.segment(i, HHO->nReconstructUnknowns).dot(HHOElement(e)->ApplyReconstructMassMatrix(v2.segment(i, HHO->nReconstructUnknowns)));
		}
		return total;
	}

	double Integral(const Vector& reconstructedCoeffs) override
	{
		assert(reconstructedCoeffs.rows() == Dimension());

		double total = 0;
		#pragma omp parallel for reduction(+:total)
		for (Element<Dim>* e : _mesh->Elements)
		{
			auto i = e->Number * HHO->nReconstructUnknowns;
			total += HHOElement(e)->IntegralReconstruct(reconstructedCoeffs.segment(i, HHO->nReconstructUnknowns));
		}
		return total;
	}

	double Integral(DomFunction func) override
	{
		double total = 0;
		#pragma omp parallel for reduction(+:total)
		for (Element<Dim>* e : _mesh->Elements)
			total += e->Integral(func);
		return total;
	}

	/*
	SparseMatrix StiffnessMatrix()
	{
		ThreadLocalCoeffs coeffs(_mesh->Elements.size(), HHO->nReconstructUnknowns * HHO->nReconstructUnknowns);
		#pragma omp parallel for
		for (Element<Dim>* e : _mesh->Elements)
		{
			Diff_HHOElement<Dim>* elem = this->HHOElement(e);
			coeffs.Local().Add(e->Number * HHO->nReconstructUnknowns, e->Number * HHO->nReconstructUnknowns, e->IntegralGradGradMatrix(elem->ReconstructionBasis));
		}

		SparseMatrix stiff(HHO->nTotalReconstructUnknowns, HHO->nTotalReconstructUnknowns);
		coeffs.Fill(stiff);
		return stiff;
	}*/
};