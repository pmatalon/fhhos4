#pragma once
#include "Point.h"
#include "../Utils/Types.h"

class Vertex : public DomPoint
{
public:
	BigNumber Number;

	Vertex(BigNumber number, double x) : DomPoint(x)
	{
		this->Number = number;
	}
	Vertex(BigNumber number, double x, double y) : DomPoint(x, y)
	{
		this->Number = number;
	}
	Vertex(BigNumber number, double x, double y, double z) : DomPoint(x, y, z)
	{
		this->Number = number;
	}
};