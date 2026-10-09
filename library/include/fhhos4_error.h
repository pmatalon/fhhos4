// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Pierre Matalon
#pragma once
#include <stdexcept>

// Part of the public header <fhhos4>, also included by the internals of fhhos4 (Utils::FatalError()) without the rest.

// Symbols exported by the shared library fhhos4_AMG, which is compiled with -fvisibility=hidden: everything else
// (fhhos4's internals, its copy of Eigen's functions) stays private to the library.
#ifndef FHHOS4_API
#define FHHOS4_API __attribute__((visibility("default")))
#endif

namespace fhhos4
{
	// Thrown by fhhos4 on any error: invalid input, divergence of a solver, unmanaged option... The program fhhos4
	// prints its message and exits with EXIT_FAILURE (see main.cpp).
	class FHHOS4_API Error : public std::runtime_error
	{
	public:
		using std::runtime_error::runtime_error;
	};
}
