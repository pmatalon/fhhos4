#pragma once
#include "IDiscreteSpace.h"
#include "../Diff_HHOElement.h"
#include "../../../Utils/Communicator.h"


template<int Dim>
class HHOCellSpace : public IDiscreteSpace
{
private:
	const Mesh<Dim>* _mesh = nullptr;
	vector<Diff_HHOElement<Dim>>* _hhoElements = nullptr;
	HHOParameters<Dim>* HHO = nullptr;
	const vector<Element<Dim>*>* _elements = nullptr;    // those of the computations of the process (the vectors are 0 elsewhere)
	const vector<Element<Dim>*>* _ownElements = nullptr; // those of the sums, reduced over the processes

public:
	HHOCellSpace()
	{}

	HHOCellSpace(const Mesh<Dim>* mesh, HHOParameters<Dim>* hho, vector<Diff_HHOElement<Dim>>& hhoElements)
		: HHOCellSpace(mesh, hho, hhoElements, mesh->Elements, mesh->Elements)
	{}

	HHOCellSpace(const Mesh<Dim>* mesh, HHOParameters<Dim>* hho, vector<Diff_HHOElement<Dim>>& hhoElements, const vector<Element<Dim>*>& elements, const vector<Element<Dim>*>& ownElements)
	{
		HHO = hho;
		_mesh = mesh;
		_hhoElements = &hhoElements;
		_elements = &elements;
		_ownElements = &ownElements;
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
		return _mesh->Elements.size() * HHO->nCellUnknowns;
	}

	double Measure() override
	{
		return _mesh->Measure();
	}

	Vector InnerProdWithBasis(DomFunction func) override
	{
		Vector innerProds = Vector::Zero(Dimension());
		#pragma omp parallel for
		for (Element<Dim>* e : *_elements)
		{
			Diff_HHOElement<Dim>* elem = HHOElement(e);
			BigNumber i = e->Number * HHO->nCellUnknowns;

			innerProds.segment(i, HHO->nCellUnknowns) = HHOElement(e)->InnerProductWithBasis(elem->CellBasis, func);
		}
		return innerProds;
	}

	Vector ApplyMassMatrix(const Vector& v) override
	{
		assert(v.rows() == Dimension());

		if (HHO->OrthonormalizeElemBases())
			return v;
		
		Vector res = Vector::Zero(v.rows());
		#pragma omp parallel for
		for (Element<Dim>* e : *_elements)
		{
			BigNumber i = e->Number * HHO->nCellUnknowns;
			res.segment(i, HHO->nCellUnknowns) = HHOElement(e)->ApplyCellMassMatrix(v.segment(i, HHO->nCellUnknowns));
		}
		return res;
	}

	Vector SolveMassMatrix(const Vector& v) override
	{
		assert(v.rows() == Dimension());

		if (HHO->OrthonormalizeElemBases())
			return v;

		Vector res = Vector::Zero(v.rows());
		#pragma omp parallel for
		for (Element<Dim>* e : *_elements)
		{
			BigNumber i = e->Number * HHO->nCellUnknowns;
			res.segment(i, HHO->nCellUnknowns) = HHOElement(e)->SolveCellMassMatrix(v.segment(i, HHO->nCellUnknowns));
		}
		return res;
	}
	
	Vector Project(DomFunction func) override
	{
		Vector vectorOfDoFs = Vector::Zero(Dimension());
		#pragma omp parallel for
		for (Element<Dim>* e : *_elements)
			vectorOfDoFs.segment(e->Number * HHO->nCellUnknowns, HHO->nCellUnknowns) = HHOElement(e)->ProjectOnCellBasis(func);
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
		for (Element<Dim>* e : *_ownElements)
		{
			BigNumber i = e->Number * HHO->nCellUnknowns;
			total += v1.segment(i, HHO->nCellUnknowns).dot(HHOElement(e)->ApplyCellMassMatrix(v2.segment(i, HHO->nCellUnknowns)));
		}
		return Communicator::Current().Sum(total);
	}

	double Integral(const Vector& cellCoeffs) override
	{
		assert(cellCoeffs.rows() == Dimension());

		double total = 0;
		#pragma omp parallel for reduction(+:total)
		for (Element<Dim>* e : *_ownElements)
		{
			auto i = e->Number * HHO->nCellUnknowns;
			total += HHOElement(e)->IntegralCell(cellCoeffs.segment(i, HHO->nCellUnknowns));
		}
		return Communicator::Current().Sum(total);
	}

	double Integral(DomFunction func) override
	{
		double total = 0;
		#pragma omp parallel for reduction(+:total)
		for (Element<Dim>* e : *_ownElements)
			total += e->Integral(func);
		return Communicator::Current().Sum(total);
	}

	//----------------------------------//
	// Functions specific to this space //
	//----------------------------------//

	Vector Solve_A_T_T(const Vector& b_T)
	{
		assert(b_T.rows() == Dimension());

		Vector result = Vector::Zero(b_T.rows());
		#pragma omp parallel for
		for (Element<Dim>* e : *_elements)
		{
			auto i = e->Number * HHO->nCellUnknowns;
			result.segment(i, HHO->nCellUnknowns) = HHOElement(e)->AttSolver.solve(b_T.segment(i, HHO->nCellUnknowns));
		}
		return result;
	}

	SparseMatrix Solve_A_T_T(const SparseMatrix& A_T_ndF)
	{
		ThreadLocalCoeffs coeffs(_elements->size(), HHO->nCellUnknowns * HHO->nCellUnknowns);
		#pragma omp parallel for
		for (Element<Dim>* e : *_elements)
		{
			auto i = e->Number * HHO->nCellUnknowns;
			for (auto f : e->Faces)
			{
				if (f->HasDirichletBC())
					continue;
				auto j = f->Number * HHO->nFaceUnknowns;
				DenseMatrix block_A_T_ndF = A_T_ndF.block(i, j, HHO->nCellUnknowns, HHO->nFaceUnknowns);
				coeffs.Local().Add(i, j, HHOElement(e)->AttSolver.solve(block_A_T_ndF));
			}
		}
		SparseMatrix result = SparseMatrix(HHO->nTotalCellUnknowns, HHO->nTotalFaceUnknowns);
		coeffs.Fill(result);
		return result;
	}
};