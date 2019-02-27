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
	SourceFunction* _sourceFunction;
	DiffusionPartition _diffusionPartition;
public:
	std::map<StandardElementCode, Poisson_DG_ReferenceElement<Dim>*> ReferenceElements;

	Poisson_DGTerms(SourceFunction* sourceFunction, FunctionalBasis<Dim>* basis, DiffusionPartition diffusionPartition)
		: _diffusionPartition(diffusionPartition)
	{
		this->_sourceFunction = sourceFunction;

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

	virtual double VolumicTerm(Element<Dim>* element, BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		Poisson_DG_ReferenceElement<Dim>* referenceElement = this->ReferenceElements[element->StdElementCode()];
		Poisson_DG_Element<Dim>* dgElement = dynamic_cast<Poisson_DG_Element<Dim>*>(element);

		return dgElement->VolumicTerm(phi1, phi2, referenceElement, this->_diffusionPartition);
	}

	virtual double MassTerm(Element<Dim>* element, BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		Poisson_DG_ReferenceElement<Dim>* referenceElement = this->ReferenceElements[element->StdElementCode()];
		Poisson_DG_Element<Dim>* dgElement = dynamic_cast<Poisson_DG_Element<Dim>*>(element);

		return dgElement->MassTerm(phi1, phi2, referenceElement);
	}

	virtual double CouplingTerm(Face<Dim>* face, Element<Dim>* element1, BasisFunction<Dim>* phi1, Element<Dim>* element2, BasisFunction<Dim>* phi2)
	{
		Poisson_DG_Face<Dim>* dgFace = dynamic_cast<Poisson_DG_Face<Dim>*>(face);
		//Poisson_DG_Element<Dim>* dgElement1 = dynamic_cast<Poisson_DG_Element<Dim>*>(element1);
		//Poisson_DG_Element<Dim>* dgElement2 = dynamic_cast<Poisson_DG_Element<Dim>*>(element2);
		return dgFace->CouplingTerm(element1, phi1, element2, phi2, this->_diffusionPartition);
	}

	virtual double PenalizationTerm(Face<Dim>* face, Element<Dim>* element1, BasisFunction<Dim>* phi1, Element<Dim>* element2, BasisFunction<Dim>* phi2, double penalizationCoefficient)
	{
		Poisson_DG_Face<Dim>* dgFace = dynamic_cast<Poisson_DG_Face<Dim>*>(face);
		Poisson_DG_Element<Dim>* dgElement1 = dynamic_cast<Poisson_DG_Element<Dim>*>(element1);
		Poisson_DG_Element<Dim>* dgElement2 = dynamic_cast<Poisson_DG_Element<Dim>*>(element2);
		return dgFace->PenalizationTerm(dgElement1, phi1, dgElement2, phi2, penalizationCoefficient, this->_diffusionPartition);
	}

	virtual double RightHandSide(Element<Dim>* element, BasisFunction<Dim>* phi)
	{
		Poisson_DG_Element<Dim>* dgElement = dynamic_cast<Poisson_DG_Element<Dim>*>(element);
		return dgElement->SourceTerm(phi, this->_sourceFunction);
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