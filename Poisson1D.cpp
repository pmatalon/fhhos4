#include "Poisson1D.h"
#include <iostream>
#include "ModalBasis.h"
#include "FileMatrix.h"
#include "FileVector.h"
#include <math.h>

using namespace std;

Poisson1D::Poisson1D(int n, function<double(double)> sourceFunction)
{
	cout << "----------------------------------------" << endl;
	cout << "Poisson1D(n=" << n << ")" << endl;
	cout << "----------------------------------------" << endl;
	this->_n = n;
	this->_sourceFunction = sourceFunction;

	// Grid initialization:
	// [0,1] descretized in 0, 1/n, 2/n, n/n (=> n+1 points)

	this->_x = new double[n+1];
	for (int k = 0; k < n + 1; k++)
		this->_x[k] = (double)k / n;

	cout << "Grid: [0, 1] --> " << (n+1) << " points (" << (n-1) << " interior points + 2 boundary points)" << endl;
}

void Poisson1D::DiscretizeDG(int maxPolynomialDegree, int penalizationCoefficient)
{
	cout << "Discretization: Discontinuous Galerkin SIPG" << endl;
	cout << "\tPolynomial degree: " << maxPolynomialDegree << endl;
	cout << "\tPenalization coefficient: " << penalizationCoefficient << endl;
	cout << "\tBasis of polynomials: monomials" << endl;

	// _n subintervals of [0, 1], maxPolynomialDegree + 1 unknowns per subinterval ==> _n * (maxPolynomialDegree + 1) unknowns
	int nUnknowns = this->_n * (maxPolynomialDegree + 1);
	cout << "Unknowns: " << nUnknowns << endl;

	// Modal basis: 1, X, X^2, ..., X^p
	ModalBasis* basis = new ModalBasis(maxPolynomialDegree, this->_n, penalizationCoefficient, this->_sourceFunction);

	string path = "/mnt/c/Users/pierr/Desktop";
	string fileName = "Poisson1D_n" + to_string(this->_n) + "_DG_SIPG_p" + to_string(maxPolynomialDegree) + "_pen" + to_string(penalizationCoefficient);
	string matrixFilePath = path + "/" + fileName + "_A.dat";
	FileMatrix* fileMatrix = new FileMatrix(nUnknowns, nUnknowns, matrixFilePath);

	string rhsFilePath = path + "/" + fileName + "_b.dat";
	FileVector* fileRHS = new FileVector(rhsFilePath);

	for (int k = 0; k < this->_n; k++) // interval [k/n, (k+1)/n]
	{
		for (int degree = 0; degree < maxPolynomialDegree + 1; degree++)
		{
			int basisFunction1 = k * (maxPolynomialDegree + 1) + degree;

			if (k > 0)
			{
				// Left neighbour (k-1)
				for (int degreeLeft = 0; degreeLeft < maxPolynomialDegree + 1; degreeLeft++)
				{
					int basisFunction2 = (k - 1) * (maxPolynomialDegree + 1) + degreeLeft;
					double couplingTerm = basis->CouplingTerm(this->_x[k], basisFunction1, basisFunction2);
					fileMatrix->Add(basisFunction1 + 1, basisFunction2 + 1, couplingTerm);
				}
			}

			// Current element (block diagonal)
			for (int degreeCurrent = 0; degreeCurrent < maxPolynomialDegree + 1; degreeCurrent++)
			{
				int basisFunction2 = k * (maxPolynomialDegree + 1) + degreeCurrent;
				double volumicTerm = basis->VolumicTerm(basisFunction1, basisFunction2, this->_x[k], this->_x[k + 1]);
				fileMatrix->Add(basisFunction1 + 1, basisFunction2 + 1, volumicTerm);
			}

			if (k < this->_n - 1)
			{
				// Right neighbour (k+1)
				for (int degreeRight = 0; degreeRight < maxPolynomialDegree + 1; degreeRight++)
				{
					int basisFunction2 = (k + 1) * (maxPolynomialDegree + 1) + degreeRight;
					double couplingTerm = basis->CouplingTerm(this->_x[k+1], basisFunction1, basisFunction2);
					fileMatrix->Add(basisFunction1 + 1, basisFunction2 + 1, couplingTerm);
				}
			}

			double rhs = basis->RightHandSide(this->_x[k], this->_x[k + 1], basisFunction1);
			fileRHS->Add(rhs);
		}
	}

	delete fileMatrix;
	delete fileRHS;

	cout << "Matrix exported to " << matrixFilePath << endl;
	cout << "RHS exported to " << rhsFilePath << endl;
}


Poisson1D::~Poisson1D()
{
	delete[] this->_x;
}
