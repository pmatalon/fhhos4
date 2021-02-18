#pragma once
#include "IterativeSolver.h"
using namespace std;

#ifdef AGMG_ENABLED
extern "C" void dagmg_(int*, double*, int*, int*, double*, double*, int*, int*, int*, int*, double*);
#endif

enum class AGMGJob : int 
{
	Setup_SolveWithoutInitialGuess_MemoryRelease = 0,
	Setup_SolveWithInitialGuess_MemoryRelease = 10,
	Setup = 1,
	SolveWithoutInitialGuess = 202, // should be  2 (problem with RowMajor/ColorMajor)
	SolveWithInitialGuess = 212,    // should be 12 (problem with RowMajor/ColorMajor)
	MemoryRelease = -1
};

class AGMG : public IterativeSolver
{
public:
	AGMG() : IterativeSolver()
	{
#ifndef AGMG_ENABLED
		Utils::Error("AGMG is disabled. The source must be recompiled with the appropriate option to use AGMG.");
#endif
	}

	void Serialize(ostream& os) const override
	{
		os << "AGMG";
	}

	void Setup(const SparseMatrix& A) override
	{
		IterativeSolver::Setup(A);

		int n = A.cols();
		int nnz = A.nonZeros();

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
			a[i] = A.valuePtr()[i];
			ja[i] = A.innerIndexPtr()[i] + 1;
		}
		for (int i = 0; i < n; i++)
			ia[i] = A.outerIndexPtr()[i] + 1;
		ia[n] = A.outerIndexPtr()[n] + 1; // last value

		int job = (int)AGMGJob::Setup;
		int iprint = 6; // unit number for output messages: 6 => standard output
		int nrest = 1; // 1 because we want to use flexible CG (the matrix is symmetric positive definite)
		int iterations;
		double tolerance;
#ifdef AGMG_ENABLED
		dagmg_(&n, a, ja, ia, nullptr, nullptr, &job, &iprint, &nrest, &iterations, &tolerance);
#endif
		delete ia, ja, a;
	}

	void Solve(const Vector& b, Vector& initialGuess, bool zeroInitialGuess) override
	{
		// TODO: Manage initialGuess

		int n = this->Matrix->cols();

		/*        first allocate the vectors with correct size */
		double *f, *x;
		f = new double[n];
		x = new double[n];

		for (int i = 0; i < n; i++)
			f[i] = b[i];

		if (!zeroInitialGuess)
		{
			for (int i = 0; i < n; i++)
				x[i] = initialGuess[i];
		}

		int job = zeroInitialGuess ? (int)AGMGJob::SolveWithoutInitialGuess : (int)AGMGJob::SolveWithInitialGuess;
		int iprint = this->PrintIterationResults ? 6 : -1; // unit number for output messages: 6 => standard output
		int nrest = 1; // 1 because we want to use flexible CG (the matrix is symmetric positive definite)
		int iterations = this->MaxIterations;
		double tolerance = this->Tolerance;
#ifdef AGMG_ENABLED
		dagmg_(&n, nullptr, nullptr, nullptr, f, x, &job, &iprint, &nrest, &iterations, &tolerance);
#endif
		// copy into solution
		for (int i = 0; i < n; i++)
			initialGuess[i] = x[i];

		this->IterationCount = iterations;

		delete f, x;
	}
};