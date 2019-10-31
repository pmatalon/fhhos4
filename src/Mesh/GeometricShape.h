#pragma once
#include "../Utils/Utils.h"
#include "../FunctionalBasis/FunctionalBasis.h"

template <int Dim>
class GeometricShape
{
public:
	GeometricShape() {}

	//-----------------------//
	//   Virtual functions   //
	//-----------------------//

	// Geometric information
	virtual double Diameter() const = 0;
	virtual double Measure() const = 0;
	virtual DomPoint Center() const = 0;

	// Integrals
	virtual double Integral(RefFunction func) const = 0;
	virtual double Integral(RefFunction func, int polynomialDegree) const = 0;
	virtual double Integral(DomFunction globalFunction) const = 0;
	virtual double Integral(DomFunction globalFunction, int polynomialDegree) const = 0;

	virtual double Integral(BasisFunction<Dim>* phi) const
	{
		RefFunction func = [phi](RefPoint p) {
			return phi->Eval(p);
		};
		return Integral(func, phi->GetDegree());
	}

	virtual void Serialize(ostream& os) const = 0;
};