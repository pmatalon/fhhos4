#pragma once
#include <Eigen/Sparse>
#include "Solver.h"
#include "enable_AGMG.h"
using namespace std;

#ifdef TPL_ENABLE_AGMG
extern "C" void dagmg_(int*, double*, int*, int*, double*, double*, int*, int*, int*, int*, double*);
#endif

class AGMG : public Solver
{
private:
	SparseMatrix A;

public:
	double Tolerance;

	AGMG(double tolerance) : Solver() 
	{
#ifdef TPL_ENABLE_AGMG
		this->Tolerance = tolerance;
#else
		assert(false && "AGMG is disabled. The source must be recompiled with the appropriate option to use AGMG.");
#endif
	}

	void Serialize(ostream& os) const override
	{
		os << "AGMG";
	}

	void Setup(const SparseMatrix& A) override
	{
		this->A = A;
	}

	Vector Solve(const Vector& b) override
	{
#ifdef TPL_ENABLE_AGMG
		int ijob = 0, nrest = 1;
		int nun = this->A.cols();
		int nnz = this->A.nonZeros();

		/*      maximal number of iterations */
		int iter = 100;
		/*      tolerance on relative residual norm */
		double tol = this->Tolerance;
		/*      unit number for output messages: 6 => standard output */
		int iprint = 6;

		/*      generate the matrix in required format (CSR) */
		/*        first allocate the vectors with correct size */
		double *a, *f, *x;
		int *ja, *ia;
		ia = new int[nun + 1];
		ja = new int[nnz];
		a = new double[nnz];
		f = new double[nun];
		x = new double[nun];

		/*        next set entries */
		for (int i = 0; i < nnz; i++)
		{
			a[i] = this->A.valuePtr()[i];
			ja[i] = this->A.innerIndexPtr()[i] + 1;
		}
		for (int i = 0; i < nun; i++)
		{
			ia[i] = this->A.outerIndexPtr()[i] + 1;
			f[i] = b[i];
		}
		ia[nun] = this->A.outerIndexPtr()[nun] + 1; // last value

		/*      call agmg
		argument 5 (ijob)  is 0 because we want a complete solve
		argument 7 (nrest) is 1 because we want to use flexible CG
						 (the matrix is symmetric positive definite) */
		dagmg_(&nun, a, ja, ia, f, x, &ijob, &iprint, &nrest, &iter, &tol);

		// copy into solution
		Vector solution = Vector(nun);
		for (int i = 0; i < nun; i++)
			solution[i] = x[i];

		delete ia, ja, a, f, x;
		return solution;
#endif
	}
};