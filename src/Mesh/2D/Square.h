#pragma once
#include "Rectangle.h"

class Square : public Rectangle
{
public:
	Square(int number, double x, double y, double width) :
		Element(number),
		Rectangle(number, x, y, width, width)
	{}
};