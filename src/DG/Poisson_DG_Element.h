#pragma once
#include "../Mesh/Element.h"
#include "../Problem/SourceFunction.h"

template <int Dim>
class Poisson_DG_Element : virtual public Element<Dim>
{
public:
	Poisson_DG_Element(BigNumber number) : Element<Dim>(number) {}

	double VolumicTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		return this->Kappa * this->StiffnessTerm(phi1, phi2);
	}

	double StiffnessTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		return this->Shape()->StiffnessTerm(phi1, phi2);
	}

	double MassTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		return this->Shape()->MassTerm(phi1, phi2);
	}

	virtual ~Poisson_DG_Element() {}
};