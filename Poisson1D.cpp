#include "Poisson1D.h"
#include <iostream>
#include "ModalBasis.h"

using namespace std;


Poisson1D::Poisson1D(int n, int maxPolynomialDegree)
{
	cout << "new Poisson1D(" << n << ")" << endl;
	this->_n = n;
	this->_polyDegree = maxPolynomialDegree;

	// Grid initialization:
	// [0,1] descretized in 0, 1/n, 2/n, n/n (=> n+1 points)

	this->_x = new double[n+1];
	for (int i = 0; i < n + 1; i++)
		this->_x[i] = (double)i / n;
}

void Poisson1D::DiscretizeDG()
{
	// _n subintervals of [0, 1], _polyDegree + 1 unknowns per subinterval ==> _n * (_polyDegree + 1) unknowns
	int nUnknowns = this->_n * (this->_polyDegree + 1);
	// Modal basis: 1, X, X^2, ..., X^p
	ModalBasis* basis = new ModalBasis(this->_polyDegree, this->_n);



	for (int k = 0; k < this->_n; k++) // interval [k/n, (k+1)/n]
	{
		for (int p = 0; p < this->_polyDegree + 1; p++)
		{
			
		}
	}
	//int nnz = 
}


Poisson1D::~Poisson1D()
{
	delete[] this->_x;
}
