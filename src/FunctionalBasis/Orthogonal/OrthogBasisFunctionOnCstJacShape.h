#pragma once
#include "OrthogonalBasisFunction.h"


template <int Dim>
class OrthogBasisFunctionOnCstJacShape : public OrthogonalBasisFunction<Dim>
{
private:
	const OrthogBasisFunctionOnGeoShape<Dim>* _refShapeFunction;
	double _detJacobian;
	bool _normalize;

public:
	OrthogBasisFunctionOnCstJacShape(const OrthogBasisFunctionOnGeoShape<Dim>* phi, double detJacobian, bool normalize)
	{
		this->LocalNumber = phi->LocalNumber;
		_refShapeFunction = phi;
		_detJacobian = detJacobian;
		_normalize = normalize;
	}

	bool IsNormalized() const
	{
		return _normalize;
	}

	double NormSquare() const override 
	{ 
		return _normalize ? 1 : _refShapeFunction->NormSquare() * _detJacobian;
	}

	double Eval(const RefPoint& p) const override
	{
		double valueOnRefShape = _refShapeFunction->Eval(p);
		if (_normalize)
			valueOnRefShape /= sqrt(_detJacobian);
		return valueOnRefShape;
	}

	DimVector<Dim> Grad(const RefPoint& p) const override
	{
		DimVector<Dim> valueOnRefShape = _refShapeFunction->Grad(p);
		if (_normalize)
			valueOnRefShape /= sqrt(_detJacobian);
		return valueOnRefShape;
	}

	int GetDegree() const override { return _refShapeFunction->GetDegree(); }

	string ToString() override { return ""; }
};