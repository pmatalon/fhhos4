#pragma once
#include "../Mesh/Element.h"
#include "../Utils/SourceFunction.h"

template <int Dim>
class Poisson_DG_Element : virtual public Element<Dim>
{
public:
	Poisson_DG_Element(BigNumber number) : Element<Dim>(number) {}

	virtual double VolumicTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		return this->Kappa * this->StiffnessTerm(phi1, phi2);
	}

	virtual double StiffnessTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) = 0;
	virtual double MassTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) = 0;
	virtual double SourceTerm(BasisFunction<Dim>* phi, SourceFunction* f) = 0;

	virtual ~Poisson_DG_Element() {}
};