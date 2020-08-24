#pragma once
#include <string>
#include <vector>
using namespace std;

class PhysicalGroup
{
public:
	int Id;
	string Name;

	PhysicalGroup(int id)
	{
		this->Id = id;
	}
};