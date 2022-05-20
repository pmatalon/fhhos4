#pragma once
#include "BasisFunction.h"
#include <math.h>
using namespace std;

//----------//
//    2D    //
//----------//

class TensorPolynomial2D : public IBasisFunction2D
{
public:
	IBasisFunction1D* FuncX;
	IBasisFunction1D* FuncY;

	TensorPolynomial2D(int localNumber, IBasisFunction1D* funcX, IBasisFunction1D* funcY)
	{
		this->LocalNumber = localNumber;
		this->FuncX = funcX;
		this->FuncY = funcY;
	}

	int GetDegree() const
	{
		return this->FuncX->GetDegree() + this->FuncY->GetDegree();
	}

	double Eval(double x, double y) const
	{
		return this->FuncX->Eval(x) * this->FuncY->Eval(y);
	}

	double EvalGradX(double x, double y) const
	{
		return this->FuncX->EvalDerivative(x) * this->FuncY->Eval(y);
	}

	double EvalGradY(double x, double y) const
	{
		return this->FuncX->Eval(x) * this->FuncY->EvalDerivative(y);
	}

	string ToString()
	{
		string polyX = this->FuncX->ToString("X");
		string polyY = this->FuncY->ToString("Y");
		if (polyX.compare("1") == 0)
			return to_string(this->LocalNumber) + "\tdegree " + to_string(this->GetDegree()) + "\t" + polyY;
		if (polyY.compare("1") == 0)
			return to_string(this->LocalNumber) + "\tdegree " + to_string(this->GetDegree()) + "\t" + polyX;
		return to_string(this->LocalNumber) + "\tdegree " + to_string(this->GetDegree()) + "\t" + polyX + " * " + polyY;
	}
};

//----------//
//    3D    //
//----------//

#ifdef ENABLE_3D

class TensorPolynomial3D : public IBasisFunction3D
{
private:
	IBasisFunction1D* FuncX;
	IBasisFunction1D* FuncY;
	IBasisFunction1D* _funcZ;
public:

	TensorPolynomial3D(int localNumber, IBasisFunction1D* funcX, IBasisFunction1D* funcY, IBasisFunction1D* funcZ)
	{
		this->LocalNumber = localNumber;
		this->FuncX = funcX;
		this->FuncY = funcY;
		this->_funcZ = funcZ;
	}

	int GetDegree() const
	{
		return this->FuncX->GetDegree() + this->FuncY->GetDegree() + this->_funcZ->GetDegree();
	}

	double Eval(double x, double y, double z) const
	{
		return this->FuncX->Eval(x) * this->FuncY->Eval(y) * this->_funcZ->Eval(z);
	}

	double EvalGradX(double x, double y, double z) const
	{
		return this->FuncX->EvalDerivative(x) * this->FuncY->Eval(y) * this->_funcZ->Eval(z);
	}

	double EvalGradY(double x, double y, double z) const
	{
		return this->FuncX->Eval(x) * this->FuncY->EvalDerivative(y) * this->_funcZ->Eval(z);
	}

	double EvalGradZ(double x, double y, double z) const
	{
		return this->FuncX->Eval(x) * this->FuncY->Eval(y) * this->_funcZ->EvalDerivative(z);
	}

	string ToString()
	{
		string polyX = this->FuncX->ToString("X");
		string polyY = this->FuncY->ToString("Y");
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

#endif // ENABLE_3D