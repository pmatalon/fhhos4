#pragma once
#include <map>
#include "../Mesh/Element.h"
#include "../Mesh/Face.h"
#include "../FunctionalBasis/FunctionalBasis.h"
#include "Poisson_DG_ReferenceElement.h"
#include "Poisson_DG_ReferenceInterval.h"
#include "Poisson_DG_ReferenceSquare.h"
#include "Poisson_DG_ReferenceCube.h"
#include "Poisson_DG_Element.h"
#include "Poisson_DG_Face.h"

template <short Dim>
class Poisson_DGTerms
{
private:
	DiffusionPartition _diffusionPartition;
public:
	std::map<StandardElementCode, Poisson_DG_ReferenceElement<Dim>*> ReferenceElements;

	Poisson_DGTerms(FunctionalBasis<Dim>* basis, DiffusionPartition diffusionPartition)
		: _diffusionPartition(diffusionPartition)
	{
		if (Dim == 1)
		{
			Poisson_DG_ReferenceElement<Dim>* refInterval = dynamic_cast<Poisson_DG_ReferenceElement<Dim>*>(new Poisson_DG_ReferenceInterval(basis->NumberOfLocalFunctionsInElement(NULL)));
			this->ComputeReferenceTerms(basis, refInterval);
			this->ReferenceElements.insert(std::make_pair(StandardElementCode::Interval, refInterval));
		}
		else if (Dim == 2)
		{
			Poisson_DG_ReferenceElement<Dim>* refSquare = dynamic_cast<Poisson_DG_ReferenceElement<Dim>*>(new Poisson_DG_ReferenceSquare(basis->NumberOfLocalFunctionsInElement(NULL)));
			this->ComputeReferenceTerms(basis, refSquare);
			this->ReferenceElements.insert(std::make_pair(StandardElementCode::Square, refSquare));
		}
		else if (Dim == 3)
		{
			Poisson_DG_ReferenceElement<Dim>* refCube = dynamic_cast<Poisson_DG_ReferenceElement<Dim>*>(new Poisson_DG_ReferenceCube(basis->NumberOfLocalFunctionsInElement(NULL)));
			this->ComputeReferenceTerms(basis, refCube);
			this->ReferenceElements.insert(std::make_pair(StandardElementCode::Cube, refCube));
		}
	}

	virtual bool IsGlobalBasis() { return false; }

	virtual double VolumicTerm(Poisson_DG_Element<Dim>* element, BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		Poisson_DG_ReferenceElement<Dim>* referenceElement = this->ReferenceElements[element->StdElementCode()];
		return element->VolumicTerm(phi1, phi2, referenceElement, this->_diffusionPartition);
	}

	virtual double MassTerm(Poisson_DG_Element<Dim>* element, BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		Poisson_DG_ReferenceElement<Dim>* referenceElement = this->ReferenceElements[element->StdElementCode()];
		return element->MassTerm(phi1, phi2, referenceElement);
	}

	virtual ~Poisson_DGTerms() {}

protected:

	void ComputeReferenceTerms(FunctionalBasis<Dim>* basis, Poisson_DG_ReferenceElement<Dim>* element)
	{
		for (BasisFunction<Dim>* phi1 : basis->LocalFunctions)
		{
			for (BasisFunction<Dim>* phi2 : basis->LocalFunctions)
			{
				element->ComputeVolumicTerm(phi1, phi2);
				element->ComputeMassTerm(phi1, phi2);
			}
		}
	}
};