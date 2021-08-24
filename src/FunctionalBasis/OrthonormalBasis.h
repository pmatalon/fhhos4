#pragma once
#include "FunctionalBasis.h"
#include "../Geometry/PhysicalShape.h"


template <int Dim>
class OrthonormalBasisFunction : public BasisFunction<Dim>
{
	vector<BasisFunction<Dim>*> _functions;
	vector<double> _coeffs;

public:
	OrthonormalBasisFunction(BasisFunction<Dim>* phi)
	{
		this->LocalNumber = phi->LocalNumber;
		_functions.push_back(phi);
		_coeffs.push_back(1.0);
	}

	void Minus(double coeff, OrthonormalBasisFunction<Dim>* previousPhi)
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

	virtual int GetDegree() override
	{
		return _functions[0]->GetDegree();
	}

	virtual string ToString() override
	{
		return "orthonormalized " + _functions[0]->ToString();
	}
};



template <int Dim>
class OrthonormalBasis : public FunctionalBasis<Dim>
{
public:
	OrthonormalBasis(FunctionalBasis<Dim>* basis, PhysicalShape<Dim>* shape) :
		FunctionalBasis<Dim>()
	{
		this->_maxPolynomialDegree = basis->GetDegree();
		this->_basisCode = "orthonorm_" + basis->BasisCode();
		this->IsHierarchical = basis->IsHierarchical;
		this->IsOrthogonal = true;
		Orthonormalize(basis, shape);
	}

private:
	// Modified Gram-Schmitt algorithm with reorthogonalization
	// (Giraud et al., The loss of orthogonality in the Gram-Schmidt orthogonalization process, 2003)
	void Orthonormalize(FunctionalBasis<Dim>* basis, PhysicalShape<Dim>* shape)
	{
		for (int i = 0; i < basis->LocalFunctions.size(); i++)
		{
			OrthonormalBasisFunction<Dim>* phi = new OrthonormalBasisFunction<Dim>(basis->LocalFunctions[i]);
			for (int nOrthogonalization = 0; nOrthogonalization < 2; nOrthogonalization++) // 2 passes of orthogonalization
			{
				for (int j = 0; j < i; j++)
				{
					OrthonormalBasisFunction<Dim>* previousPhi = dynamic_cast<OrthonormalBasisFunction<Dim>*>(this->LocalFunctions[j]);
					double innerprod = shape->ComputeMassTerm(phi, previousPhi);
					phi->Minus(innerprod, previousPhi);
				}
			}
			double norm = shape->L2Norm(phi);
			phi->DivideBy(norm);
			this->LocalFunctions.push_back(phi);
		}
	}
};