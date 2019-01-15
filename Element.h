#pragma once
#include <vector>
#include "Utils.h"
class Face;

enum class StandardElementCode
{
	None,
	Interval,
	Square,
	Cube
};

class Element 
{
public:
	BigNumber Number;
	std::vector<Face*> Faces;

	Element(BigNumber number)
	{
		this->Number = number;
	}

	virtual StandardElementCode StdElementCode() = 0;

	virtual double* OuterNormalVector(Face* interface) = 0;

	virtual ~Element() {}
};