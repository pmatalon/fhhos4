#pragma once
#include "Solver.h"
using namespace std;

#ifdef AGMG_ENABLED
extern "C" void dagmg_(int*, double*, int*, int*, double*, double*, int*, int*, int*, int*, double*);
#endif

enum class AGMGJob : int 
{
	Setup_SolveWithoutInitialGuess_MemoryRelease = 0,
	Setup_SolveWithInitialGuess_MemoryRelease = 10,
	Setup = 1,
	SolveWithoutInitialGuess = 202,
	SolveWithInitialGuess = 212,
	MemoryRelease = -1
};

class AGMG : public Solver
{
private:
	const SparseMatrix* A;

public:
	double Tolerance;
	int MaxIterations;

	AGMG(double tolerance, int maxIterations) : Solver() 
	{
#ifdef AGMG_ENABLED
		this->Tolerance = tolerance;
		this->MaxIterations = maxIterations;
#else
		Utils::FatalError("AGMG is disabled. The source must be recompiled with the appropriate option to use AGMG.");
#endif
	}

	void Serialize(ostream& os) const override
	{
		os << "AGMG";
	}

	void Setup(const SparseMatrix& A) override
	{
		this->A = &A;

#ifdef AGMG_ENABLED
		int n = this->A->cols();
		int nnz = this->A->nonZeros();

		/*      generate the matrix in required format (CSR) */
		/*        first allocate the vectors with correct size */
		double *a;
		int *ja, *ia;
		ia = new int[n + 1];
		ja = new int[nnz];
		a = new double[nnz];

		/*        next set entries */
		for (int i = 0; i < nnz; i++)
		{
			a[i] = this->A->valuePtr()[i];
			ja[i] = this->A->innerIndexPtr()[i] + 1;
		}
		for (int i = 0; i < n; i++)
			ia[i] = this->A->outerIndexPtr()[i] + 1;
		ia[n] = this->A->outerIndexPtr()[n] + 1; // last value

		int job = (int)AGMGJob::Setup;
		int iprint = 6; // unit number for output messages: 6 => standard output
		int nrest = 1; // 1 because we want to use flexible CG (the matrix is symmetric positive definite)
		int iterations;
		double tolerance;
		cout << "before" << endl;
		dagmg_(&n, a, ja, ia, nullptr, nullptr, &job, &iprint, &nrest, &iterations, &tolerance);
		cout << "after" << endl;

		delete ia, ja, a;
#endif
	}

	Vector Solve(const Vector& b) override
	{
#ifdef AGMG_ENABLED
		int n = this->A->cols();

		/*      generate the matrix in required format (CSR) */
		/*        first allocate the vectors with correct size */
		double *f, *x;
		f = new double[n];
		x = new double[n];

		for (int i = 0; i < n; i++)
			f[i] = b[i];

		int job = (int)AGMGJob::SolveWithoutInitialGuess;
		int iprint = 6; // unit number for output messages: 6 => standard output
		int nrest = 1; // 1 because we want to use flexible CG (the matrix is symmetric positive definite)
		int iterations = this->MaxIterations;
		double tolerance = this->Tolerance;
		dagmg_(&n, nullptr, nullptr, nullptr, f, x, &job, &iprint, &nrest, &iterations, &tolerance);

		// copy into solution
		Vector solution = Vector(n);
		for (int i = 0; i < n; i++)
			solution[i] = x[i];

		delete f, x;
		return solution;
#endif
	}
};