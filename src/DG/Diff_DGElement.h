#pragma once
#include "../Mesh/Element.h"

template <int Dim>
class Diff_DGElement : virtual public Element<Dim>
{
public:
	Diff_DGElement(BigNumber number) : Element<Dim>(number) {}

	double VolumicTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) const
	{
		return this->Kappa() * this->StiffnessTerm(phi1, phi2);
	}

	double StiffnessTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) const
	{
		return this->Shape()->StiffnessTerm(phi1, phi2);
	}

	double MassTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) const
	{
		return this->Shape()->MassTerm(phi1, phi2);
	}

	virtual ~Diff_DGElement() {}
};