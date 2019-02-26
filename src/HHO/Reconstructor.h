#pragma once
#include <Eigen/Sparse>
#include "Poisson_HHO_Element.h"

template <short Dim>
class Reconstructor
{
private:
	Poisson_HHO_Element<Dim>* _element;
	FunctionalBasis<Dim>* _reconstructionBasis;
	FunctionalBasis<Dim>* _elementBasis;
	FunctionalBasis<Dim - 1>* _faceBasis;
public:
	Eigen::MatrixXd S;
	Eigen::MatrixXd Bt;
	Eigen::MatrixXd Bf;

	Reconstructor(Poisson_HHO_Element<Dim>* element, FunctionalBasis<Dim>* reconstructionBasis, FunctionalBasis<Dim>* elementBasis, FunctionalBasis<Dim - 1>* faceBasis) :
		S(reconstructionBasis->LocalFunctions.size(), reconstructionBasis->LocalFunctions.size()), 
		Bt(reconstructionBasis->LocalFunctions.size(), elementBasis->LocalFunctions.size()),
		Bf(reconstructionBasis->LocalFunctions.size(), faceBasis->LocalFunctions.size())
	{
		this->element = element;
		this->_reconstructionBasis = reconstructionBasis;
		this->_elementBasis = elementBasis;
		this->_faceBasis = faceBasis;
	}

	void AssembleS()
	{
		for (BasisFunction<Dim>* phi1 : this->_reconstructionBasis->LocalFunctions)
		{
			for (BasisFunction<Dim>* phi2 : this->_reconstructionBasis->LocalFunctions)
			{
				this->S(phi1->LocalNumber, phi2->LocalNumber) = _element->S(phi1, phi2);
			}
		}
	}

	void AssembleBt()
	{
		for (BasisFunction<Dim>* phi1 : this->_reconstructionBasis->LocalFunctions)
		{
			for (BasisFunction<Dim>* phi2 : this->_elementBasis->LocalFunctions)
			{
				this->Bt(phi1->LocalNumber, phi2->LocalNumber) = _element->Bt(phi1, phi2);
			}
		}
	}

	void AssembleBf()
	{
		for (BasisFunction<Dim>* phi1 : this->_reconstructionBasis->LocalFunctions)
		{
			for (BasisFunction<Dim-1>* phi2 : this->_faceBasis->LocalFunctions)
			{
				this->Bt(phi1->LocalNumber, phi2->LocalNumber) = _element->Bf(phi1, phi2);
			}
		}
	}
};