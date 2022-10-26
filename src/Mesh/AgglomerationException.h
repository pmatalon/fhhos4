#pragma once
#include <exception>
#include <string>
using namespace std;

class AgglomerationException : public exception
{
private:
	string _error;
public:
	AgglomerationException(string error)
	{
		_error = "AgglomerationException: " + error;
	}

	virtual const char* what() const throw()
	{
		return _error.c_str();
	}
};