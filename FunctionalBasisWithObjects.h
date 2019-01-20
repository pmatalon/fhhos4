#pragma once
#include <string>
#include <vector>
#include <map>
#include "Element.h"
#include "Face.h"
#include "IBasisFunction.h"

class FunctionalBasisWithObjects
{
public:
	vector<BasisFunction*> LocalFunctions;
	virtual std::string Name() = 0;

	virtual int GetDegree() = 0;

	int NumberOfLocalFunctionsInElement(Element* element)
	{
		return static_cast<int>(this->LocalFunctions.size());
	}

	BigNumber GlobalFunctionNumber(Element* element, BasisFunction* phi)
	{
		return element->Number * static_cast<int>(this->LocalFunctions.size()) + phi->LocalNumber; // the numbers start at 0
	}

	virtual ~FunctionalBasisWithObjects() {}
};