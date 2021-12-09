#pragma once
#include "DiffusionTestCase.h"
using namespace std;

template <int Dim>
class VirtualDiffusionTestCase : public DiffusionTestCase<Dim>
{
public:
	VirtualDiffusionTestCase() {}

	VirtualDiffusionTestCase(DomFunction sourceFunction, DiffusionField<Dim>& diffField)
	{
		this->SourceFunction = sourceFunction;
		this->DiffField = diffField;
	}

	string Code() override
	{
		return "virtual";
	}
	string Description() override
	{
		return "Virtual diffusion test case";
	}
};
