#pragma once
#include "BasisFunction.h"
#include <math.h>
#include <algorithm>
#include <array>
#include <vector>
using namespace std;

// The bases of tensor polynomials (monomials, Legendre) span Q^qDegree + P^pDegree: the products of 1D polynomials whose
// largest partial degree is at most qDegree, or whose degree (the sum of the partial degrees) is at most pDegree. P^p:
// qDegree 0; Q^q: pDegree 0; their sum: the reconstruction of HHO in Q (Program_Diffusion_HHO). The functions are ordered
// by largest partial degree up to qDegree (the first (d + 1)^Dim span Q^d), then by degree: hierarchical in Q, and the
// bases of P keep their order (by degree, then by decreasing partial degree in x).

//----------//
//    2D    //
//----------//

// The partial degrees (in x, y) of the functions of the basis of Q^qDegree + P^pDegree, in their order
inline vector<array<int, 2>> TensorPartialDegrees2D(int qDegree, int pDegree)
{
	vector<array<int, 2>> degrees;
	for (int d = 0; d <= qDegree; d++)
	{
		for (int j = 0; j <= d; j++)
		{
			for (int i = 0; i <= d; i++)
			{
				if (max(i, j) == d)
					degrees.push_back({ i, j });
			}
		}
	}
	for (int degree = qDegree + 1; degree <= pDegree; degree++)
	{
		for (int j = 0; j <= degree; j++)
		{
			int i = degree - j;
			if (max(i, j) > qDegree)
				degrees.push_back({ i, j });
		}
	}
	return degrees;
}

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

	// The largest partial degree: the degree in the sense of Q
	int GetQDegree() const
	{
		return max(this->FuncX->GetDegree(), this->FuncY->GetDegree());
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

// The partial degrees (in x, y, z) of the functions of the basis of Q^qDegree + P^pDegree, in their order
inline vector<array<int, 3>> TensorPartialDegrees3D(int qDegree, int pDegree)
{
	vector<array<int, 3>> degrees;
	for (int d = 0; d <= qDegree; d++)
	{
		for (int k = 0; k <= d; k++)
		{
			for (int j = 0; j <= d; j++)
			{
				for (int i = 0; i <= d; i++)
				{
					if (max({ i, j, k }) == d)
						degrees.push_back({ i, j, k });
				}
			}
		}
	}
	for (int degree = qDegree + 1; degree <= pDegree; degree++)
	{
		for (int degZ = 0; degZ <= degree; degZ++)
		{
			for (int degY = 0; degY <= degree - degZ; degY++)
			{
				int degX = degree - degZ - degY;
				if (max({ degX, degY, degZ }) > qDegree)
					degrees.push_back({ degX, degY, degZ });
			}
		}
	}
	return degrees;
}

class TensorPolynomial3D : public IBasisFunction3D
{
public:
	IBasisFunction1D* FuncX;
	IBasisFunction1D* FuncY;
	IBasisFunction1D* FuncZ;

	TensorPolynomial3D(int localNumber, IBasisFunction1D* funcX, IBasisFunction1D* funcY, IBasisFunction1D* funcZ)
	{
		this->LocalNumber = localNumber;
		this->FuncX = funcX;
		this->FuncY = funcY;
		this->FuncZ = funcZ;
	}

	int GetDegree() const
	{
		return this->FuncX->GetDegree() + this->FuncY->GetDegree() + this->FuncZ->GetDegree();
	}

	// The largest partial degree: the degree in the sense of Q
	int GetQDegree() const
	{
		return max({ this->FuncX->GetDegree(), this->FuncY->GetDegree(), this->FuncZ->GetDegree() });
	}

	double Eval(double x, double y, double z) const
	{
		return this->FuncX->Eval(x) * this->FuncY->Eval(y) * this->FuncZ->Eval(z);
	}

	double EvalGradX(double x, double y, double z) const
	{
		return this->FuncX->EvalDerivative(x) * this->FuncY->Eval(y) * this->FuncZ->Eval(z);
	}

	double EvalGradY(double x, double y, double z) const
	{
		return this->FuncX->Eval(x) * this->FuncY->EvalDerivative(y) * this->FuncZ->Eval(z);
	}

	double EvalGradZ(double x, double y, double z) const
	{
		return this->FuncX->Eval(x) * this->FuncY->Eval(y) * this->FuncZ->EvalDerivative(z);
	}

	string ToString()
	{
		string polyX = this->FuncX->ToString("X");
		string polyY = this->FuncY->ToString("Y");
		string polyZ = this->FuncZ->ToString("Z");
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