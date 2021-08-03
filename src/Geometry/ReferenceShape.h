#pragma once
#include "GeometricShape.h"

template <int Dim>
class ReferenceShape : public GeometricShape<Dim>
{
private:
	map<FunctionalBasis<Dim>*, DenseMatrix> _massMatrices;
	map<FunctionalBasis<Dim>*, DenseMatrix> _cellReconstructMatrices;

public:
	ReferenceShape() {}

	double InRadius() const override
	{
		assert(false);
	}

	virtual vector<RefPoint> QuadraturePoints() const = 0;
	virtual vector<RefPoint> QuadraturePoints(int polynomialDegree) const = 0;

	virtual double Integral(RefFunction func) const = 0;
	virtual double Integral(RefFunction func, int polynomialDegree) const = 0;

	double Integral(DomFunction globalFunction) const override
	{
		assert(false);
	}
	double Integral(DomFunction globalFunction, int polynomialDegree) const override
	{
		assert(false);
	}

	virtual void Serialize(ostream& os) const
	{
		assert(false);
	}

	//--------//
	//   DG   //
	//--------//

	double MassTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) const
	{
		assert(_massMatrices.size() == 1);
		_massMatrices.begin()->second(phi1->LocalNumber, phi2->LocalNumber);
	}

	//---------//
	//   HHO   //
	//---------//

	const DenseMatrix& StoredMassMatrix(FunctionalBasis<Dim>* basis) const
	{
		auto it = _massMatrices.find(basis);
		if (it != _massMatrices.end())
			return it->second;
		Utils::FatalError("The mass matrix for this basis should be computed once and stored for this reference element.");
	}

	const DenseMatrix& StoredCellReconstructMassMatrix(FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim>* reconstructBasis) const
	{
		auto it = _cellReconstructMatrices.find(cellBasis);
		if (it != _cellReconstructMatrices.end())
			return it->second;
		Utils::FatalError("The cell-reconstruct mass matrix should be computed once and stored for this reference element.");
	}

	void ComputeAndStoreMassMatrix(FunctionalBasis<Dim>* basis)
	{
		if (_massMatrices.find(basis) == _massMatrices.end())
			_massMatrices[basis] = this->ComputeAndReturnMassMatrix(basis);
	}
	
	void ComputeAndStoreCellReconstructMassMatrix(FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim>* reconstructBasis)
	{
		if (_cellReconstructMatrices.find(cellBasis) == _cellReconstructMatrices.end())
			_cellReconstructMatrices[cellBasis] = this->ComputeAndReturnMassMatrix(cellBasis, reconstructBasis);
	}

};