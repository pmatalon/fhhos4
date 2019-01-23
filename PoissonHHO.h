#pragma once
#include "Problem.h"
#include "IMesh.h"
#include "FunctionalBasisWithObjects.h"
#include "ElementInterface.h"

template <class IBasisFunction>
class PoissonHHO : Problem
{
	PoissonHHO(string solutionName) : Problem(solutionName)
	{	}

	void Discretize(IMesh* mesh, FunctionalBasis<IBasisFunction>* basis, IPoisson_HHOTerms<IBasisFunction>* hho)
	{

	}
};