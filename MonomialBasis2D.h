#pragma once
#include "FunctionalBasis2D.h"
#include "BasisFunction2D.h"
//#include "Monomial2D.h"

using namespace std;

class Monomial2D : public BasisFunction2D
{
public:
	int DegreeX;
	int DegreeY;

	Monomial2D(int degreeX, int degreeY)
	{
		this->DegreeX = degreeX;
		this->DegreeY = degreeY;
	}

	double Eval(double x, double y)
	{
		return pow(x, this->DegreeX)*pow(y, this->DegreeY);
	}

	double EvalGradX(double x, double y)
	{
		if (this->DegreeX == 0)
			return 0;
		return this->DegreeX*pow(x, this->DegreeX - 1)*pow(y, this->DegreeY);
	}

	double EvalGradY(double x, double y)
	{
		if (this->DegreeY == 0)
			return 0;
		return this->DegreeY*pow(y, this->DegreeY - 1)*pow(x, this->DegreeX);
	}
};

class MonomialBasis2D : public FunctionalBasis2D
{
private:
	int _maxPolynomialDegree;

public:
	MonomialBasis2D(int maxPolynomialDegree, CartesianGrid2D* grid, int penalizationCoefficient, function<double(double, double)> sourceFunction)
		:FunctionalBasis2D(grid, penalizationCoefficient, sourceFunction)
	{
		this->_maxPolynomialDegree = maxPolynomialDegree;

		int functionNumber = 0;
		for (int i = 0; i <= this->_maxPolynomialDegree; i++)
		{
			for (int j = 0; i + j <= this->_maxPolynomialDegree; j++)
				this->_localFunctions[functionNumber++] = new Monomial2D(i, j);
		}
	}

	std::string Name()
	{
		return "monomials_p" + std::to_string(this->_maxPolynomialDegree);
	}

	/*double VolumicTerm(Element* element, int localFunctionNumber1, int localFunctionNumber2)
	{
		Square* square = static_cast<Square*>(element);
		return this->VolumicTerm(square, localFunctionNumber1, localFunctionNumber2);
	}

	double VolumicTerm(Square* element, int localFunctionNumber1, int localFunctionNumber2)
	{
		Monomial2D* func1 = (Monomial2D*)this->_localFunctions[localFunctionNumber1];
		Monomial2D* func2 = (Monomial2D*)this->_localFunctions[localFunctionNumber2];

		int i1 = func1->DegreeX;
		int j1 = func1->DegreeY;
		int i2 = func2->DegreeX;
		int j2 = func2->DegreeY;

		double integralDerivX = 0;
		if (i1 == 0 || i2 == 0)
			integralDerivX = 0;
		else
		{
			double integralOverX = (double)(i1 * i2) / (double)(i1 + i2 - 1) * (pow(element->X + element->Width, i1 + i2 - 1) - pow(element->X, i1 + i2 - 1));
			double integralOverY = 1 / (double)(j1 + j2 + 1) * (pow(element->Y + element->Width, j1 + j2 + 1) - pow(element->Y, j1 + j2 + 1));
			integralDerivX = integralOverX * integralOverY;
		}

		double integralDerivY = 0;
		if (j1 == 0 || j2 == 0)
			integralDerivY = 0;
		else
		{
			double integralOverY = (double)(j1 * j2) / (double)(j1 + j2 - 1) * (pow(element->Y + element->Width, j1 + j2 - 1) - pow(element->Y, j1 + j2 - 1));
			double integralOverX = 1 / (double)(i1 + i2 + 1) * (pow(element->X + element->Width, i1 + i2 + 1) - pow(element->X, i1 + i2 + 1));
			integralDerivY = integralOverY * integralOverX;
		}

		return integralDerivX + integralDerivY;
	}*/

	/*double CouplingTerm(ElementInterface* interface, Element* element1, int localFunctionNumber1, Element* element2, int localFunctionNumber2)
	{
		Square* square1 = static_cast<Square*>(element1);
		Square* square2 = static_cast<Square*>(element2);
		return this->CouplingTerm(interface, square1, localFunctionNumber1, square2, localFunctionNumber2);
	}*/

	/*double CouplingTerm(ElementInterface* interface, Square* element1, int localFunctionNumber1, Square* element2, int localFunctionNumber2)
	{
		if (!interface->IsBetween(element1, element2))
			return 0;

		Monomial2D* func1 = (Monomial2D*)this->_localFunctions[localFunctionNumber1];
		Monomial2D* func2 = (Monomial2D*)this->_localFunctions[localFunctionNumber2];

		int i1 = func1->DegreeX;
		int j1 = func1->DegreeY;
		int i2 = func2->DegreeX;
		int j2 = func2->DegreeY;

		double* jumpGradFunc1 = new double[2]{ 0, 0 };
		double* jumpGradFunc2 = new double[2]{ 0, 0 };
		double* meanFunc1 = new double[2]{ 0, 0 };
		double* meanFunc2 = new double[2]{ 0, 0 };

		auto n1 = element1->OuterNormalVector(interface);
		auto n2 = element2->OuterNormalVector(interface);

		double integral = 0;
		if (interface == element1->NorthInterface)
		{
			double y = element1->Y + element1->Width;
			if (i1 != 0)
				jumpGradFunc1
			integral = pow(y, j1 + j2) / (i1 + i2 + 1)*(pow(element1->X + element1->Width, i1 + i2 + 1) - pow(element1->X, i1 + i2 + 1)) * (n1[0] * n2[0] + n1[1] * n2[1]);
		}
		else if (interface == element1->SouthInterface)
		{
			double y = element1->Y;
		}
		else if (interface == element1->EastInterface)
		{
			double x = element1->X + element1->Width;
		}
		else if (interface == element1->WestInterface)
		{
			double x = element1->X;
		}
		else
			return 0;
		//return MeanGrad(element1, func1, interface) * Jump(element2, func2, interface) + MeanGrad(element2, func2, interface) * Jump(element1, func1, interface);
	}*/

	/*double PenalizationTerm(ElementInterface* interface, Element* element1, int localFunctionNumber1, Element* element2, int localFunctionNumber2)
	{
		Square* square1 = static_cast<Square*>(element1);
		Square* square2 = static_cast<Square*>(element2);
		return this->PenalizationTerm(interface, square1, localFunctionNumber1, square2, localFunctionNumber2);
	}

	double PenalizationTerm(ElementInterface* interface, Square* element1, int localFunctionNumber1, Square* element2, int localFunctionNumber2)
	{
		//if (!interface->IsBetween(element1, element2))
		//	return 0;

		Monomial2D* func1 = (Monomial2D*)this->_localFunctions[localFunctionNumber1];
		Monomial2D* func2 = (Monomial2D*)this->_localFunctions[localFunctionNumber2];

		int i1 = func1->DegreeX;
		int j1 = func1->DegreeY;
		int i2 = func2->DegreeX;
		int j2 = func2->DegreeY;

		auto n1 = element1->OuterNormalVector(interface);
		auto n2 = element2->OuterNormalVector(interface);

		double integralJump1ScalarJump2 = 0;
		if (interface == element1->NorthInterface)
		{
			double y = element1->Y + element1->Width;
			integralJump1ScalarJump2 = pow(y, j1 + j2) / (i1 + i2 + 1)*(pow(element1->X + element1->Width, i1 + i2 + 1) - pow(element1->X, i1 + i2 + 1)) * (n1[0] * n2[0] + n1[1] * n2[1]);
		}
		else if (interface == element1->SouthInterface)
		{
			double y = element1->Y;
			integralJump1ScalarJump2 = pow(y, j1 + j2) / (i1 + i2 + 1)*(pow(element1->X + element1->Width, i1 + i2 + 1) - pow(element1->X, i1 + i2 + 1)) * (n1[0] * n2[0] + n1[1] * n2[1]);
		}
		else if (interface == element1->EastInterface)
		{
			double x = element1->X + element1->Width;
			integralJump1ScalarJump2 = pow(x, i1 + i2) / (j1 + j2 + 1)*(pow(element1->Y + element1->Width, j1 + j2 + 1) - pow(element1->Y, j1 + j2 + 1)) * (n1[0] * n2[0] + n1[1] * n2[1]);
		}
		else if (interface == element1->WestInterface)
		{
			double x = element1->X;
			integralJump1ScalarJump2 = pow(x, i1 + i2) / (j1 + j2 + 1)*(pow(element1->Y + element1->Width, j1 + j2 + 1) - pow(element1->Y, j1 + j2 + 1)) * (n1[0] * n2[0] + n1[1] * n2[1]);
		}
		else
			return 0;
		return this->_penalizationCoefficient * integralJump1ScalarJump2;
	}*/

	/*double EvalPrimitiveOfTheProductWithConstantY(double y)
	{
		integralJump1ScalarJump2 = pow(y, j1 + j2) / (i1 + i2 + 1)*(pow(element1->X + element1->Width, i1 + i2 + 1) - pow(element1->X, i1 + i2 + 1))
	}*/

	/*function<double(double, double)> PrimitiveProductConstantY(Monomial2D* func1, Monomial2D* func2, double y)
	{
		int i1 = func1->DegreeX;
		int j1 = func1->DegreeY;
		int i2 = func2->DegreeX;
		int j2 = func2->DegreeY;

		function<double(double)> primitive = [y](double x) { 
			return ; 
		};

		return primitive;
	}*/
};