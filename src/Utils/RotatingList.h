#pragma once
using namespace std;

template <typename T>
class RotatingList
{
private:
	vector<T> _list;
	int _i = 0;
public:
	RotatingList(vector<T> v) : _list(v) {}

	size_t Size() { return _list.size(); }

	T Get()
	{ 
		return _list[_i];
	}

	T GetNext()
	{
		int j = _i + 1;
		if (j == _list.size())
			j = 0;
		return _list[j];
	}

	T GetPrevious()
	{
		int j = _i - 1;
		if (j < 0)
			j = _list.size() - 1;
		return _list[j];
	}

	size_t Index() { return _i; }

	void MoveNext()
	{
		_i++;
		if (_i == _list.size())
			_i = 0;
	}

	void MoveBack()
	{
		_i--;
		if (_i < 0)
			_i = _list.size() - 1;
	}

	T GetAndMoveNext()
	{
		T obj = Get();
		MoveNext();
		return obj;
	}

	void Reset() { _i = 0; }

	void GoTo(size_t index) { _i = index; }

	void GoTo(T obj)
	{
		for (int j = 0; j < _list.size(); j++)
		{
			if (_list[j] == obj)
			{
				_i = j;
				return;
			}
		}
		assert(false && "Object not found in the list");
	}
};