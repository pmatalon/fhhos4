#pragma once
#include <cstdint>

// The number of functions of P^degree (the polynomials of degree at most degree), or of Q^degree (at most degree in each
// variable) if polynomialSpaceQ, in dim variables. Standalone: also compiled into the preconditioner for PETSc
// (Solver/DiscretizationInference.h).
inline std::int64_t PolynomialSpaceDimension(int dim, int degree, bool polynomialSpaceQ)
{
	std::int64_t n = 1;
	for (int i = 1; i <= dim; i++)
		n = polynomialSpaceQ ? n * ((std::int64_t)degree + 1) : n * ((std::int64_t)degree + i) / i; // C(degree + i, i): exact at each step
	return n;
}
