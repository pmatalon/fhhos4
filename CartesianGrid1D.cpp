#include "CartesianGrid1D.h"

CartesianGrid1D::CartesianGrid1D(int n)
{
	this->_n = n;

	// [0,1] descretized in 0, 1/n, 2/n, n/n (=> n+1 points)
	this->_x = new double[n + 1];
	for (int k = 0; k < n + 1; k++)
		this->_x[k] = (double)k / n;
}

//inline CartesianGrid1D::~CartesianGrid1D()
//{
//	delete[] this->_x;
//}
//
//inline int CartesianGrid1D::NElements()
//{
//	return this->_n;
//}
//
//inline double CartesianGrid1D::X(int point)
//{
//	return this->_x[point];
//}
//
//inline double CartesianGrid1D::XRight(int element)
//{
//	return this->_x[element + 1];
//}
//
//inline double CartesianGrid1D::XLeft(int element)
//{
//	return this->_x[element];
//}
//
//inline int CartesianGrid1D::GetInterface(int element1, int element2)
//{
//	if (element1 == element2 + 1)
//		return element1;
//	if (element1 == element2 - 1)
//		return element2;
//	return -1;
//}
//
//inline int CartesianGrid1D::LeftInterface(int element)
//{
//	return element;
//}
//
//inline int CartesianGrid1D::RightInterface(int element)
//{
//	return element + 1;
//}
//
//inline bool CartesianGrid1D::IsLeftInterface(int element, int point)
//{
//	return point == element;
//}
//
//inline bool CartesianGrid1D::IsRightInterface(int element, int point)
//{
//	return point == element + 1;
//}
//
//inline bool CartesianGrid1D::IsBoundaryLeft(int point)
//{
//	return point == 0;
//}
//inline bool CartesianGrid1D::IsBoundaryRight(int point)
//{
//	return point == this->_n;
//}
//
//inline bool CartesianGrid1D::IsFirstElement(int element)
//{
//	return element == 0;
//}
//inline bool CartesianGrid1D::IsLastElement(int element)
//{
//	return element == this->_n - 1;
//}
