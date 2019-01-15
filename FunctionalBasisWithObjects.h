#pragma once
#include <string>
#include <vector>
#include <map>
#include "Element.h"
#include "ElementInterface.h"

template <class IBasisFunction>
class FunctionalBasisWithObjects
{
protected:
	map<int, IBasisFunction*> _localFunctions;
	//vector<IBasisFunction*> _localFunctions;

public:
	virtual std::string Name() = 0;

	virtual int GetDegree() = 0;

	int NumberOfLocalFunctionsInElement(Element* element)
	{
		return static_cast<int>(this->_localFunctions.size());
	}

	IBasisFunction* GetLocalBasisFunction(Element* element, int localFunctionNumber)
	{
		return this->_localFunctions[localFunctionNumber];
	}

	BigNumber GlobalFunctionNumber(Element* element, int localFunctionNumber)
	{
		return element->Number * static_cast<int>(this->_localFunctions.size()) + localFunctionNumber; // the numbers start at 0
	}

	virtual ~FunctionalBasisWithObjects() {}
};