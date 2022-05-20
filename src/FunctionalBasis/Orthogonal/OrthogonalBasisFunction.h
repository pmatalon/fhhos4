#pragma once
#include "../BasisFunction.h"


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