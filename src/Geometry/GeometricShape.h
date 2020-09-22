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

	virtual void Serialize(ostream& os) const = 0;

	// Geometric information
	virtual double Diameter() const = 0;
	virtual double Measure() const = 0;
	virtual DomPoint Center() const = 0;
	virtual double InRadius() const = 0;

	// Integrals
	virtual double Integral(RefFunction func) const = 0;
	virtual double Integral(RefFunction func, int polynomialDegree) const = 0;
	virtual double Integral(DomFunction globalFunction) const = 0;
	virtual double Integral(DomFunction globalFunction, int polynomialDegree) const = 0;

	//--------------------------------------------------------------------------------//

	double Regularity() const
	{
		return 2 * InRadius() / Diameter();
	}

	virtual ~GeometricShape() {}

	virtual double Integral(BasisFunction<Dim>* phi) const
	{
		RefFunction func = [phi](const RefPoint& p) {
			return phi->Eval(p);
		};
		return Integral(func, phi->GetDegree());
	}

protected:
	DenseMatrix ComputeAndReturnMassMatrix(FunctionalBasis<Dim>* basis)
	{
		DenseMatrix M = DenseMatrix(basis->Size(), basis->Size());
		for (BasisFunction<Dim>* phi1 : basis->LocalFunctions)
		{
			for (BasisFunction<Dim>* phi2 : basis->LocalFunctions)
			{
				if (phi2->LocalNumber > phi1->LocalNumber)
					break;
				double term = ComputeMassTerm(phi1, phi2);
				M(phi1->LocalNumber, phi2->LocalNumber) = term;
				M(phi2->LocalNumber, phi1->LocalNumber) = term;
			}
		}
		return M;
	}

	DenseMatrix ComputeAndReturnMassMatrix(FunctionalBasis<Dim>* basis1, FunctionalBasis<Dim>* basis2)
	{
		DenseMatrix M(basis1->LocalFunctions.size(), basis2->LocalFunctions.size());
		for (BasisFunction<Dim>* phi1 : basis1->LocalFunctions)
		{
			for (BasisFunction<Dim>* phi2 : basis2->LocalFunctions)
			{
				double term = ComputeMassTerm(phi1, phi2);
				M(phi1->LocalNumber, phi2->LocalNumber) = term;
			}
		}
		return M;
	}

	double ComputeMassTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) const
	{
		RefFunction functionToIntegrate = [phi1, phi2](const RefPoint& p) {
			return phi1->Eval(p)*phi2->Eval(p);
		};

		int polynomialDegree = phi1->GetDegree() + phi2->GetDegree();
		return Integral(functionToIntegrate, polynomialDegree);
	}

	//-------------------------------------------------------------------//
	//                            Unit tests                             //
	//-------------------------------------------------------------------//

public:
	virtual void UnitTests() const
	{
		RefFunction refOne = [](const RefPoint& p) { return 1; };
		DomFunction domOne = [](const DomPoint& p) { return 1; };

		double measure = Measure();

		for (int degree = 0; degree < 5; degree++)
		{
			double integral = Integral(refOne, degree);
			assert(abs(integral - measure) < 1e-4*measure);
		}
		double integral = Integral(domOne);
		assert(abs(integral - measure) < 1e-4*measure);
	}
};