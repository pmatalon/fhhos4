#pragma once
#include "../FunctionalBasis.h"
#include "LagrangeP1.h"
#include "../../Utils/Utils.h"

template <int Dim>
class LagrangeBasis : public FunctionalBasis<Dim>
{
public:
	bool   IsHierarchical()                const override { return false; }
	string BasisCode()                     const override { return "lagrange"; }
	static string Code()                                  { return "lagrange"; };
};


class LagrangeP1Basis : public LagrangeBasis<2>
{
private:
	LagrangeP1_Node1 _lagrangeNode1;
	LagrangeP1_Node2 _lagrangeNode2;
	LagrangeP1_Node3 _lagrangeNode3;
public:
	LagrangeP1Basis(int maxPolynomialDegree)
	{
		if (maxPolynomialDegree != 1)
			Utils::FatalError("Only piecewise linear approximation is implemented with the Lagrange basis.");
	}

	vector<BasisFunction<2>*> LocalFunctions() override
	{
		return { &_lagrangeNode1 , &_lagrangeNode2, &_lagrangeNode3 };
	}

	int GetDegree() const override { return 1; }
	int Size()      const override { return 3; }

	FunctionalBasis<2>* CreateSameBasisForDegree(int degree) override
	{
		assert(false);
		return nullptr;
	}

	FunctionalBasis<2>* CreateLowerDegreeBasis(int degree) override
	{
		assert(false);
		return nullptr;
	}
};