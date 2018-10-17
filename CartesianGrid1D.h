#pragma once
class CartesianGrid1D
{
private:
	int _n;
	double* _x;
public:
	CartesianGrid1D(int n);

	inline ~CartesianGrid1D()
	{
		delete[] this->_x;
	}

	inline int NElements()
	{
		return this->_n;
	}

	inline double X(int point)
	{
		return this->_x[point];
	}

	inline double XRight(int element)
	{
		return this->_x[element + 1];
	}

	inline double XLeft(int element)
	{
		return this->_x[element];
	}

	inline int GetInterface(int element1, int element2)
	{
		if (element1 == element2 + 1)
			return element1;
		if (element1 == element2 - 1)
			return element2;
		return -1;
	}

	inline int LeftInterface(int element)
	{
		return element;
	}

	inline int RightInterface(int element)
	{
		return element + 1;
	}

	inline bool IsLeftInterface(int element, int point)
	{
		return point == element;
	}

	inline bool IsRightInterface(int element, int point)
	{
		return point == element + 1;
	}

	inline bool IsBoundaryLeft(int point)
	{
		return point == 0;
	}
	inline bool IsBoundaryRight(int point)
	{
		return point == this->_n;
	}

	inline bool IsFirstElement(int element)
	{
		return element == 0;
	}
	inline bool IsLastElement(int element)
	{
		return element == this->_n - 1;
	}

};