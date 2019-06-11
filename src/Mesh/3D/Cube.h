#pragma once
#include "Parallelepiped.h"

class Cube : public Parallelepiped
{
public:
	Cube(int number, double x, double y, double z, double width) :
		Element(number),
		Parallelepiped(number, x, y, z, width, width, width)
	{}
};