#include "Element.h"
#pragma once

template <int Dim>
class Face
{
public:
	BigNumber Number;
	bool IsDomainBoundary;
	Element<Dim>* Element1;
	Element<Dim>* Element2;

public:
	Face(BigNumber number, Element<Dim>* element1, Element<Dim>* element2)
	{
		this->Number = number;
		this->Element1 = element1;
		this->Element2 = element2;
		this->IsDomainBoundary = element2 == NULL;
	}
	Face(BigNumber number, Element<Dim>* element1)
		:Face(number, element1, NULL)
	{
		this->IsDomainBoundary = true;
	}

	bool IsBetween(Element<Dim>* element1, Element<Dim>* element2)
	{
		if (element1 == element2 && (element1 == this->Element1 || element1 == this->Element2))
			return true;
		if (element1 != this->Element1 && element1 != this->Element2)
			return false;
		if (element2 != this->Element1 && element2 != this->Element2)
			return false;
		return true;
	}

	Element<Dim>* GetNeighbour(Element<Dim>* element)
	{
		if (element == this->Element1)
			return this->Element2;
		else if (element == this->Element2)
			return this->Element1;
		return NULL;
	}

	virtual double MassTerm(BasisFunction<Dim - 1>* phi1, BasisFunction<Dim - 1>* phi2) = 0;

	virtual double MassTerm(BasisFunction<Dim - 1>* facePhi, Element<Dim>* element, BasisFunction<Dim>* reconstructPhi) = 0;

	virtual double GetDiameter() = 0;
	virtual double Measure() = 0;

	virtual DomPoint ConvertToDomain(RefPoint refPoint) = 0;

	virtual Eigen::MatrixXd MassMatrix(FunctionalBasis<Dim-1>* basis)
	{
		Eigen::MatrixXd M(basis->LocalFunctions.size(), basis->LocalFunctions.size());
		for (BasisFunction<Dim-1>* phi1 : basis->LocalFunctions)
		{
			for (BasisFunction<Dim-1>* phi2 : basis->LocalFunctions)
			{
				if (phi2->LocalNumber > phi1->LocalNumber)
					break;
				double term = this->MassTerm(phi1, phi2);
				M(phi1->LocalNumber, phi2->LocalNumber) = term;
				M(phi2->LocalNumber, phi1->LocalNumber) = term;
			}
		}
		return M;
	}

	Eigen::MatrixXd MassMatrix(FunctionalBasis<Dim-1>* basis, Element<Dim>* element, FunctionalBasis<Dim>* cellBasis)
	{
		Eigen::MatrixXd M(basis->LocalFunctions.size(), cellBasis->LocalFunctions.size());
		for (BasisFunction<Dim-1>* phi1 : basis->LocalFunctions)
		{
			for (BasisFunction<Dim>* phi2 : cellBasis->LocalFunctions)
			{
				double term = this->MassTerm(phi1, element, phi2);
				M(phi1->LocalNumber, phi2->LocalNumber) = term;
			}
		}
		return M;
	}

	virtual vector<RefPoint> GetNodalPoints(FunctionalBasis<Dim-1>* basis) = 0;

	virtual string ToString()
	{
		return "Interface " + to_string(this->Number);
	}

	virtual ~Face() {}
};