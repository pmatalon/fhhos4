#pragma once
#include <time.h>
#include <chrono>
using namespace std;

class Duration
{
private:
	int _hh;
	int _mm;
	int _ss;
	int _ms;
public:
	int InMilliseconds;

	Duration(int durationInMilliseconds)
	{
		InMilliseconds = durationInMilliseconds;

		int rest = durationInMilliseconds;
		_hh = rest / (3600 * 1000);
		rest = rest - _hh * 3600 * 1000;
		_mm = rest / (60 * 1000);
		rest = rest - _mm * 60 * 1000;
		_ss = rest / 1000;
		_ms = rest - _ss * 1000;
	}

	friend ostream& operator<<(ostream& os, const Duration& d)
	{
		stringstream ss;
		if (d._hh < 10)
			ss << "0";
		ss << d._hh << ":";
		if (d._mm < 10)
			ss << "0";
		ss << d._mm << ":";
		if (d._ss < 10)
			ss << "0";
		ss << d._ss << ".";
		if (d._ms < 100)
			ss << "0";
		if (d._ms < 10)
			ss << "0";
		ss << d._ms;
		os << ss.str();
		return os;
	}
};

class Timer
{
private:
	clock_t _cpu_start;
	clock_t _cpu_stop;

	chrono::time_point<chrono::high_resolution_clock> _elapsed_start;
	chrono::time_point<chrono::high_resolution_clock> _elapsed_stop;
public:
	Timer() {}

	void Start()
	{
		_cpu_start = clock();
		_elapsed_start = chrono::high_resolution_clock::now();
	}

	void Stop()
	{
		_cpu_stop = clock();
		_elapsed_stop = chrono::high_resolution_clock::now();
	}

	Duration CPU() const
	{
		double span = (double)(_cpu_stop - _cpu_start);
		Duration d((int)((double)span / CLOCKS_PER_SEC * 1000));
		return d;
	}

	Duration Elapsed() const
	{
		double durationInMilliseconds = chrono::duration_cast<chrono::duration<double, std::milli>>(_elapsed_stop - _elapsed_start).count();
		Duration d(durationInMilliseconds);
		return d;
	}
};