#pragma once
#include <exception>
#include <string>
using namespace std;

class AgglomerationException : public exception
{
private:
	string _error;
public:
	AgglomerationException(string error) : _error(error)
	{}

	virtual const char* what() const throw()
	{
		string s = "AgglomerationException: " + _error;
		return s.c_str();
	}
};