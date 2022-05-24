#pragma once
#include "OrthogBasisFunctionOnCstJacShape.h"
#include "OrthogonalBasis.h"

template <int Dim>
class OrthogonalBasisOnCstJacShape : public FunctionalBasis<Dim>
{
private:
	vector<OrthogBasisFunctionOnCstJacShape<Dim>> _localFunctions;
public:
	OrthogonalBasis<Dim>* RefShapeBasis;

private:
	OrthogonalBasisOnCstJacShape() {}
public:
	OrthogonalBasisOnCstJacShape(OrthogonalBasis<Dim>* refShapeBasis, double detJacobian, bool normalize)
	{
		RefShapeBasis = refShapeBasis;
		_localFunctions.reserve(refShapeBasis->Size());
		for (const auto& phi : refShapeBasis->_localFunctions)
			_localFunctions.emplace_back(&phi, detJacobian, normalize);
	}

	vector<BasisFunction<Dim>*> LocalFunctions() override
	{
		vector<BasisFunction<Dim>*> list;
		list.reserve(_localFunctions.size());
		for (auto& phi : _localFunctions)
			list.push_back(&phi);
		return list;
	}

	bool IsNormalized() const
	{
		return _localFunctions[0].IsNormalized();
	}

	FunctionalBasis<Dim>* CreateSameBasisForDegree(int degree) override
	{
		assert(false && "TO IMPLEMENT");
		return nullptr;
	}

	FunctionalBasis<Dim>* CreateLowerDegreeBasis(int degree) override
	{
		OrthogonalBasisOnCstJacShape* lowerBasis = new OrthogonalBasisOnCstJacShape();
		for (auto& phi : _localFunctions)
		{
			if (phi.GetDegree() <= degree)
				lowerBasis->_localFunctions.push_back(phi);
		}
		return lowerBasis;
	}

	int  Size()                const override { return (int)_localFunctions.size(); }
	string BasisCode()         const override { assert(false); return ""; }
	bool IsHierarchical()      const override { assert(false); return false; }
	int  GetDegree()           const override { assert(false); return 0; }
	bool UsePolynomialSpaceQ() const override { assert(false); return false; }
};