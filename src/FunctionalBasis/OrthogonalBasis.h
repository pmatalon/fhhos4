#pragma once
#include "FunctionalBasis.h"
#include "../Geometry/GeometricShape.h"


template <int Dim>
class OrthogonalBasisFunction : public BasisFunction<Dim>
{
private:
	vector<BasisFunction<Dim>*> _functions;
	vector<double> _coeffs;

public:
	double NormSquare = -1;

	OrthogonalBasisFunction(BasisFunction<Dim>* phi)
	{
		this->LocalNumber = phi->LocalNumber;
		_functions.push_back(phi);
		_coeffs.push_back(1.0);
	}

	void Minus(double coeff, OrthogonalBasisFunction<Dim>* previousPhi)
	{
		for (int i = 0; i < previousPhi->_functions.size(); i++)
		{
			BasisFunction<Dim>* component = previousPhi->_functions[i];
			double coeff_component = previousPhi->_coeffs[i];
			bool found = false;
			for (int j = 1; j < this->_functions.size(); j++) // j=0 is excluded
			{
				if (component == this->_functions[j])
				{
					found = true;
					this->_coeffs[j] -= coeff * coeff_component;
					break;
				}
			}
			if (!found)
			{
				_functions.push_back(component);
				_coeffs.push_back(-coeff * coeff_component);
			}
		}
	}

	void DivideBy(double coeff)
	{
		for (int i = 0; i < _coeffs.size(); i++)
			_coeffs[i] = _coeffs[i] / coeff;
	}

	double Eval(const RefPoint& p) override
	{
		double res = 0;
		for (int i = 0; i < _functions.size(); i++)
			res += _coeffs[i] * _functions[i]->Eval(p);
		return res;
	}

	virtual DimVector<Dim> Grad(const RefPoint& p) override
	{
		DimVector<Dim> res = DimVector<Dim>::Zero();
		for (int i = 0; i < _functions.size(); i++)
			res += _coeffs[i] * _functions[i]->Grad(p);
		return res;
	}

	virtual int GetDegree() const override
	{
		return _functions[0]->GetDegree();
	}

	virtual string ToString() override
	{
		return "orthonormalized " + _functions[0]->ToString();
	}
};



template <int Dim>
class OrthogonalBasis : public FunctionalBasis<Dim>
{
private:
	vector<OrthogonalBasisFunction<Dim>> _localFunctions;
	FunctionalBasis<Dim>* _originalBasis;
public:
	OrthogonalBasis() {}

	OrthogonalBasis(FunctionalBasis<Dim>* basis, GeometricShape<Dim>* shape, int orthogonalizationSweeps = 1, bool normalize = true) :
		FunctionalBasis<Dim>()
	{
		_originalBasis = basis;
		assert(orthogonalizationSweeps > 0);
		if (normalize)
			Orthonormalize(basis, shape, orthogonalizationSweeps);
		else
			Orthogonalize(basis, shape, orthogonalizationSweeps);
	}

	vector<BasisFunction<Dim>*> LocalFunctions() override
	{
		vector<BasisFunction<Dim>*> list;
		list.reserve(_localFunctions.size());
		for (OrthogonalBasisFunction<Dim>& phi : _localFunctions)
			list.push_back(&phi);
		return list;
	}

	string BasisCode() const override
	{
		bool normalized = true;
		for (const OrthogonalBasisFunction<Dim>& phi : _localFunctions)
		{
			if (phi.NormSquare != 1)
			{
				normalized = false;
				break;
			}
		}
		if (normalized)
			return "orthonorm_" + _originalBasis->BasisCode();
		else
			return "orthogon_" + _originalBasis->BasisCode();
	}

	FunctionalBasis<Dim>* CreateSameBasisForDegree(int degree) override
	{
		assert(false && "TO IMPLEMENT");
		return nullptr;
	}

	FunctionalBasis<Dim>* CreateLowerDegreeBasis(int degree) override
	{
		OrthogonalBasis* lowerBasis = new OrthogonalBasis();
		lowerBasis->_originalBasis = _originalBasis;
		for (auto& phi : _localFunctions)
		{
			if (phi.GetDegree() <= degree)
				lowerBasis->_localFunctions.push_back(phi);
		}
		return lowerBasis;
	}

	bool IsHierarchical()      const override { return _originalBasis->IsHierarchical(); }
	int  GetDegree()           const override { return _originalBasis->GetDegree(); }
	int  Size()                const override { return (int)_localFunctions.size(); }
	bool UsePolynomialSpaceQ() const override { return _originalBasis->UsePolynomialSpaceQ(); }

private:
	// Modified Gram-Schmitt algorithm with reorthogonalization
	// (Giraud et al., The loss of orthogonality in the Gram-Schmidt orthogonalization process, 2003)
	void Orthonormalize(FunctionalBasis<Dim>* basis, GeometricShape<Dim>* shape, int orthogonalizationSweeps = 1)
	{
		auto originalFunctions = basis->LocalFunctions();
		_localFunctions.reserve(originalFunctions.size());
		for (int i = 0; i < originalFunctions.size(); i++)
		{
			_localFunctions.emplace_back(originalFunctions[i]);
			OrthogonalBasisFunction<Dim>& phi = _localFunctions.back();
			for (int nOrthogonalization = 0; nOrthogonalization < orthogonalizationSweeps; nOrthogonalization++) // possibly 2 passes of orthogonalization
			{
				for (int j = 0; j < i; j++)
				{
					OrthogonalBasisFunction<Dim>& previousPhi = _localFunctions[j];
					double innerprod = shape->ComputeMassTerm(&phi, &previousPhi);
					phi.Minus(innerprod, &previousPhi);
				}
			}
			double norm = shape->L2Norm(&phi);
			phi.DivideBy(norm);
			phi.NormSquare = 1;
		}
	}

	void Orthogonalize(FunctionalBasis<Dim>* basis, GeometricShape<Dim>* shape, int orthogonalizationSweeps = 1)
	{
		auto originalFunctions = basis->LocalFunctions();
		_localFunctions.reserve(originalFunctions.size());
		for (int i = 0; i < originalFunctions.size(); i++)
		{
			_localFunctions.emplace_back(originalFunctions[i]);
			OrthogonalBasisFunction<Dim>& phi = _localFunctions.back();
			for (int nOrthogonalization = 0; nOrthogonalization < orthogonalizationSweeps; nOrthogonalization++) // possibly 2 passes of orthogonalization
			{
				for (int j = 0; j < i; j++)
				{
					OrthogonalBasisFunction<Dim>& previousPhi = _localFunctions[j];
					double innerprod = shape->ComputeMassTerm(&phi, &previousPhi);
					phi.Minus(innerprod / previousPhi.NormSquare, &previousPhi);
				}
			}
			phi.NormSquare = shape->L2NormSquare(&phi);
		}
	}
};

/*
template <int Dim>
class ShapeWithConstantJacobianOrthogonalBasis : public FunctionalBasis<Dim>
{
private:
	vector<OrthogonalBasisFunction<Dim>> _localFunctions;

public:
	ShapeWithConstantJacobianOrthogonalBasis(const OrthogonalBasis<Dim>& refShapeBasis, double detJacobian, bool normalize)
	{
		_originalBasis = refShapeBasis->_originalBasis;
		this->_basisCode = refShapeBasis.BasisCode();
		_localFunctions.insert(this->LocalFunctions.end(), refShapeBasis.LocalFunctions.begin(), refShapeBasis.LocalFunctions.end());
		if (normalize)
		{
			for (BasisFunction<Dim>* phi : this->LocalFunctions)
			{
				OrthogonalBasisFunction<Dim>* refphi = dynamic_cast<OrthogonalBasisFunction<Dim>*>(phi);
				refphi->DivideBy(sqrt(detJacobian));
			}
		}
		else
		{
			for (BasisFunction<Dim>* phi : this->LocalFunctions)
			{
				OrthogonalBasisFunction<Dim>* refphi = dynamic_cast<OrthogonalBasisFunction<Dim>*>(phi);
				refphi->NormSquare *= detJacobian;
			}
		}
	}
};*/