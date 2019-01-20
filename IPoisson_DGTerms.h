#pragma once
#include <map>
#include "FunctionalBasisWithObjects.h"
#include "Element.h"
#include "Face.h"
#include "Poisson_DG_ReferenceElement.h"
#include "Poisson_DG_Element.h"
#include "Poisson_DG_Face.h"

template <class IBasisFunction>
class IPoisson_DGTerms
{
public:
	std::map<StandardElementCode, Poisson_DG_ReferenceElement*> ReferenceElements;

	virtual bool IsGlobalBasis() = 0;

	virtual double VolumicTerm(Element* element, BasisFunction* phi1, BasisFunction* phi2)
	{
		Poisson_DG_ReferenceElement* referenceElement = this->ReferenceElements[element->StdElementCode()];
		Poisson_DG_Element* dgElement = dynamic_cast<Poisson_DG_Element*>(element);

		return dgElement->VolumicTerm(phi1, phi2, referenceElement);
	}

	virtual double CouplingTerm(Face* face, Element* element1, BasisFunction* phi1, Element* element2, BasisFunction* phi2)
	{
		Poisson_DG_Face* dgFace = dynamic_cast<Poisson_DG_Face*>(face);
		Poisson_DG_Element* dgElement1 = dynamic_cast<Poisson_DG_Element*>(element1);
		Poisson_DG_Element* dgElement2 = dynamic_cast<Poisson_DG_Element*>(element2);
		return dgFace->CouplingTerm(dgElement1, phi1, dgElement2, phi2);
	}

	virtual double PenalizationTerm(Face* face, Element* element1, BasisFunction* phi1, Element* element2, BasisFunction* phi2, double penalizationCoefficient)
	{
		Poisson_DG_Face* dgFace = dynamic_cast<Poisson_DG_Face*>(face);
		Poisson_DG_Element* dgElement1 = dynamic_cast<Poisson_DG_Element*>(element1);
		Poisson_DG_Element* dgElement2 = dynamic_cast<Poisson_DG_Element*>(element2);
		return dgFace->PenalizationTerm(dgElement1, phi1, dgElement2, phi2, penalizationCoefficient);
	}

	virtual double RightHandSide(Element* element, IBasisFunction* func) = 0;

protected:

	void ComputeVolumicTerms(FunctionalBasisWithObjects<IBasisFunction>* basis, Poisson_DG_ReferenceElement* element)
	{
		for (int localFunctionNumber1 = 0; localFunctionNumber1 < basis->NumberOfLocalFunctionsInElement(NULL); localFunctionNumber1++)
		{
			BasisFunction* phi1 = basis->GetLocalBasisFunction(NULL, localFunctionNumber1);
			for (int localFunctionNumber2 = 0; localFunctionNumber2 < basis->NumberOfLocalFunctionsInElement(NULL); localFunctionNumber2++)
			{
				BasisFunction* phi2 = basis->GetLocalBasisFunction(NULL, localFunctionNumber2);
				element->ComputeVolumicTerm(phi1, phi2);
			}
		}
	}

	/*Poisson_DG_Element* Create(Element* element)
	{
		if (typeid(*element) == typeid(Interval))
			return new Poisson_DG_Interval(static_cast<Interval*>(element));
		if (typeid(*element) == typeid(Square))
			return new Poisson_DG_Square(static_cast<Square*>(element));
		if (typeid(*element) == typeid(Cube))
			return new Poisson_DG_Cube(static_cast<Cube*>(element));
		return NULL;
	}*/
};