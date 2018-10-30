#pragma once
class CartesianGrid1D
{
private:
	BigNumber _n;
	double* _x;
public:
	CartesianGrid1D(BigNumber n)
	{
		this->_n = n;

		// [0,1] descretized in 0, 1/n, 2/n, n/n (=> n+1 points)
		this->_x = new double[n + 1];
		for (BigNumber k = 0; k < n + 1; k++)
			this->_x[k] = (double)k / n;
	}

	inline ~CartesianGrid1D()
	{
		delete[] this->_x;
	}

	inline int NElements()
	{
		return this->_n;
	}

	inline double X(BigNumber point)
	{
		return this->_x[point];
	}

	inline double XRight(BigNumber element)
	{
		return this->_x[element + 1];
	}

	inline double XLeft(BigNumber element)
	{
		return this->_x[element];
	}

	inline int GetInterface(BigNumber element1, BigNumber element2)
	{
		if (element1 == element2 + 1)
			return element1;
		if (element1 == element2 - 1)
			return element2;
		return -1;
	}

	inline int LeftInterface(BigNumber element)
	{
		return element;
	}

	inline int RightInterface(BigNumber element)
	{
		return element + 1;
	}

	inline bool IsLeftInterface(BigNumber element, BigNumber point)
	{
		return point == element;
	}

	inline bool IsRightInterface(BigNumber element, BigNumber point)
	{
		return point == element + 1;
	}

	inline bool IsBoundaryLeft(BigNumber point)
	{
		return point == 0;
	}
	inline bool IsBoundaryRight(BigNumber point)
	{
		return point == this->_n;
	}

	inline bool IsFirstElement(BigNumber element)
	{
		return element == 0;
	}
	inline bool IsLastElement(BigNumber element)
	{
		return element == this->_n - 1;
	}

};