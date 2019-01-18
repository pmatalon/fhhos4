/*#pragma once
class DefInterval
{
public:
	static DefInterval &MinusOne_One() { static DefInterval interval(-1, 1); return interval; }
	static DefInterval &Zero_One() { static DefInterval interval(0, 1); return interval; }
	double Left;
	double Right;
	double Length;
private:
	DefInterval(double left, double right)
	{
		this->Left = left;
		this->Right = right;
		this->Length = right - left;
	}
};*/