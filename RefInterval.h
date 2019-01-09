#pragma once
class RefInterval
{
public:
	static RefInterval &MinusOne_One() { static RefInterval interval(-1, 1); return interval; }
	static RefInterval &Zero_One() { static RefInterval interval(0, 1); return interval; }
	double Left;
	double Right;
	double Length;
private:
	RefInterval(double left, double right)
	{
		this->Left = left;
		this->Right = right;
		this->Length = right - left;
	}
};