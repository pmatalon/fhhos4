#pragma once
#include "GeometricShape.h"
#include "../FunctionalBasis/Orthogonal/OrthogonalBasis.h"
using namespace std;

struct StiffnessMatrices
{
	DenseMatrix tt;
	DenseMatrix uu;
	DenseMatrix tu;
	DenseMatrix vv;
	DenseMatrix tv;
	DenseMatrix uv;
};

template <int Dim>
class ReferenceShape : public GeometricShape<Dim>
{
private:
	// Keyed by FunctionalBasis::Id, not by the address of the basis: the reference shapes are static (they live as long as
	// the process), and a basis allocated at the address of a deleted one would get its matrices (e.g. of another degree).
	map<size_t, DenseMatrix> _massMatrices;
	map<size_t, DenseMatrix> _cellReconstructMatrices;

	map<size_t, StiffnessMatrices> _stiffMatrices;

	map<size_t, OrthogonalBasis<Dim>*> _orthogBases;

	map<size_t, Vector> _integralVectors;
public:
	ReferenceShape() {}

	double InRadius() const override
	{
		assert(false);
		return 0;
	}

	virtual string Name() const = 0;

	virtual vector<RefPoint> QuadraturePoints() const = 0;
	virtual vector<RefPoint> QuadraturePoints(int polynomialDegree) const = 0;

	virtual double Integral(RefFunction func) const = 0;
	virtual double Integral(RefFunction func, int polynomialDegree) const = 0;

	double Integral(DomFunction globalFunction) const override
	{
		assert(false);
		return 0;
	}
	double Integral(DomFunction globalFunction, int polynomialDegree) const override
	{
		assert(false);
		return 0;
	}

	void Serialize(ostream& os) const override
	{
		assert(false);
	}

	//--------//
	//   DG   //
	//--------//

	double MassTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) const
	{
		assert(_massMatrices.size() == 1);
		return _massMatrices.begin()->second(phi1->LocalNumber, phi2->LocalNumber);
	}

	//---------//
	//   HHO   //
	//---------//

	const DenseMatrix& StoredMassMatrix(FunctionalBasis<Dim>* basis) const
	{
		auto it = _massMatrices.find(basis->Id);
		if (it != _massMatrices.end())
			return it->second;
		Utils::FatalError("The mass matrix for this basis should be computed once and stored for this reference element ('" + Name() + "').");
		return it->second; // to avoid warning
	}

	const DenseMatrix& StoredCellReconstructMassMatrix(FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim>* reconstructBasis) const
	{
		auto it = _cellReconstructMatrices.find(cellBasis->Id);
		if (it != _cellReconstructMatrices.end())
			return it->second;
		Utils::FatalError("The cell-reconstruct mass matrix should be computed once and stored for this reference element ('" + Name() + "').");
		return it->second; // to avoid warning
	}

	const StiffnessMatrices& StoredStiffnessMatrices(FunctionalBasis<Dim>* basis) const
	{
		auto it = _stiffMatrices.find(basis->Id);
		if (it != _stiffMatrices.end())
			return it->second;
		Utils::FatalError("The stiffness matrices should be computed once and stored for this reference element ('" + Name() + "').");
		return it->second; // to avoid warning
	}

	const Vector& StoredIntegralVector(FunctionalBasis<Dim>* basis) const
	{
		auto it = _integralVectors.find(basis->Id);
		if (it != _integralVectors.end())
			return it->second;
		Utils::FatalError("The vector of integrals for this basis should be computed once and stored for this reference element ('" + Name() + "').");
		return it->second; // to avoid warning
	}

	void ComputeAndStoreMassMatrix(FunctionalBasis<Dim>* basis)
	{
		if (_massMatrices.find(basis->Id) == _massMatrices.end())
			_massMatrices[basis->Id] = this->ComputeAndReturnMassMatrix(basis);
	}
	
	void ComputeAndStoreCellReconstructMassMatrix(FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim>* reconstructBasis)
	{
		if (_cellReconstructMatrices.find(cellBasis->Id) == _cellReconstructMatrices.end())
			_cellReconstructMatrices[cellBasis->Id] = this->ComputeAndReturnMassMatrix(cellBasis, reconstructBasis);
	}

	// Refer to http://arturo.imati.cnr.it/~marini/didattica/Metodi-engl/Intro2FEM.pdf (page 30)
	void ComputeAndStoreStiffnessMatrices(FunctionalBasis<Dim>* basis)
	{
		if (_stiffMatrices.find(basis->Id) == _stiffMatrices.end())
		{
			StiffnessMatrices stiff;
			int t = 0;
			int u = 1;
			int v = 2;
			stiff.tt = this->ComputeAndReturnGradGradMatrix(basis, t, t);
			if (Dim >= 2)
			{
				stiff.tu = this->ComputeAndReturnGradGradMatrix(basis, t, u);
				stiff.uu = this->ComputeAndReturnGradGradMatrix(basis, u, u);
#ifdef ENABLE_3D
				if (Dim == 3)
				{
					stiff.tv = this->ComputeAndReturnGradGradMatrix(basis, t, v);
					stiff.uv = this->ComputeAndReturnGradGradMatrix(basis, u, v);
					stiff.vv = this->ComputeAndReturnGradGradMatrix(basis, v, v);
				}
#endif
			}
			_stiffMatrices[basis->Id] = stiff;
		}
	}

	void ComputeAndStoreIntegralVector(FunctionalBasis<Dim>* basis)
	{
		if (_integralVectors.find(basis->Id) == _integralVectors.end())
			_integralVectors[basis->Id] = this->ComputeIntegral(basis);
	}

private:
	DenseMatrix ComputeAndReturnGradGradMatrix(FunctionalBasis<Dim>* basis, int tuv1, int tuv2) const
	{
		DenseMatrix M = DenseMatrix(basis->Size(), basis->Size());
		auto locaFunctions = basis->LocalFunctions();
		for (BasisFunction<Dim>* phi1 : locaFunctions)
		{
			for (BasisFunction<Dim>* phi2 : locaFunctions)
			{
				if (tuv1 == tuv2 && phi2->LocalNumber > phi1->LocalNumber)
					break;
				double term = ComputeGradGradTerm(phi1, phi2, tuv1, tuv2);
				M(phi1->LocalNumber, phi2->LocalNumber) = term;
				if (tuv1 == tuv2)
					M(phi2->LocalNumber, phi1->LocalNumber) = term;
			}
		}
		return M;
	}
	double ComputeGradGradTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2, int tuv1, int tuv2) const
	{
		if (phi1->GetDegree() == 0 || phi2->GetDegree() == 0)
			return 0;

		RefFunction functionToIntegrate = [phi1, phi2, tuv1, tuv2](const RefPoint& p) {
			return phi1->Grad(p)[tuv1] * phi2->Grad(p)[tuv2];
		};

		int polynomialDegree = max(0, phi1->GetDegree() + phi2->GetDegree() - 2);
		return Integral(functionToIntegrate, polynomialDegree);
	}
public:

	OrthogonalBasis<Dim>* Orthogonalize(FunctionalBasis<Dim>* basis, int orthogonalizationSweeps = 1, bool normalize = true)
	{
		OrthogonalBasis<Dim>* orthoBasis = new OrthogonalBasis<Dim>(basis, this, orthogonalizationSweeps, normalize);
		_orthogBases[basis->Id] = orthoBasis;
		return orthoBasis;
	}

	OrthogonalBasis<Dim>* OrthogonalizedBasis(FunctionalBasis<Dim>* basis) const
	{
		auto it = _orthogBases.find(basis->Id);
		if (it != _orthogBases.end())
			return it->second;
		else
			return nullptr;
	}

};