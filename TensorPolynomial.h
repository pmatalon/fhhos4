#pragma once
#include "IBasisFunction.h"
#include <math.h>
using namespace std;

//----------//
//    2D    //
//----------//

class TensorPolynomial2D : public IBasisFunction2D
{
private:
	IBasisFunction1D* _funcX;
	IBasisFunction1D* _funcY;
public:

	TensorPolynomial2D(int localNumber, IBasisFunction1D* funcX, IBasisFunction1D* funcY)
	{
		this->LocalNumber = localNumber;
		this->_funcX = funcX;
		this->_funcY = funcY;
	}

	DefInterval DefinitionInterval() { return this->_funcX->DefinitionInterval(); }

	int GetDegree()
	{
		return this->_funcX->GetDegree() + this->_funcY->GetDegree();
	}

	double Eval(double x, double y)
	{
		return this->_funcX->Eval(x) * this->_funcY->Eval(y);
	}

	double EvalGradX(double x, double y)
	{
		return this->_funcX->EvalDerivative(x) * this->_funcY->Eval(y);
	}

	double EvalGradY(double x, double y)
	{
		return this->_funcX->Eval(x) * this->_funcY->EvalDerivative(y);
	}

	string ToString()
	{
		string polyX = this->_funcX->ToString("X");
		string polyY = this->_funcY->ToString("Y");
		if (polyX.compare("1") == 0)
			return polyY;
		if (polyY.compare("1") == 0)
			return polyX;
		return polyX + " * " + polyY;
	}
};

//----------//
//    3D    //
//----------//

class TensorPolynomial3D : public IBasisFunction3D
{
private:
	IBasisFunction1D* _funcX;
	IBasisFunction1D* _funcY;
	IBasisFunction1D* _funcZ;
public:

	TensorPolynomial3D(int localNumber, IBasisFunction1D* funcX, IBasisFunction1D* funcY, IBasisFunction1D* funcZ)
	{
		this->LocalNumber = localNumber;
		this->_funcX = funcX;
		this->_funcY = funcY;
		this->_funcZ = funcZ;
	}

	DefInterval DefinitionInterval() { return this->_funcX->DefinitionInterval(); }

	int GetDegree()
	{
		return this->_funcX->GetDegree() + this->_funcY->GetDegree() + this->_funcZ->GetDegree();
	}

	double Eval(double x, double y, double z)
	{
		return this->_funcX->Eval(x) * this->_funcY->Eval(y) * this->_funcZ->Eval(z);
	}

	double EvalGradX(double x, double y, double z)
	{
		return this->_funcX->EvalDerivative(x) * this->_funcY->Eval(y) * this->_funcZ->Eval(z);
	}

	double EvalGradY(double x, double y, double z)
	{
		return this->_funcX->Eval(x) * this->_funcY->EvalDerivative(y) * this->_funcZ->Eval(z);
	}

	double EvalGradZ(double x, double y, double z)
	{
		return this->_funcX->Eval(x) * this->_funcY->Eval(y) * this->_funcZ->EvalDerivative(z);
	}

	string ToString()
	{
		string polyX = this->_funcX->ToString("X");
		string polyY = this->_funcY->ToString("Y");
		string polyZ = this->_funcZ->ToString("Z");
		if (polyX.compare("1") == 0 && polyY.compare("1") == 0)
			return polyZ;
		if (polyY.compare("1") == 0 && polyZ.compare("1") == 0)
			return polyX;
		if (polyX.compare("1") == 0 && polyZ.compare("1") == 0)
			return polyY;
		if (polyX.compare("1") == 0)
			return polyY + " * " + polyZ;
		if (polyY.compare("1") == 0)
			return polyX + " * " + polyZ; 
		if (polyZ.compare("1") == 0)
			return polyX + " * " + polyY;
		return polyX + " * " + polyY + " * " + polyZ;
	}
};
