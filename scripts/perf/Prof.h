#pragma once
// Temporary wall-clock profiling of code regions, per multigrid level (no sampling profiler in WSL2).
// Not part of the build: include it in the files to instrument, then revert them, e.g.
//
//     #include "../../../../scripts/perf/Prof.h"   // path relative to the instrumented file
//     ProfLevel = this->Number;                     // level to which the next timings are attributed
//     PROF_START(build);
//     mesh.Build();
//     PROF_STOP(build, "mesh.Build");
//
// The times are accumulated per region and level, and printed at exit, sorted by total (wall-clock seconds).
// Only time regions executed by one thread at a time (not inside a parallel loop).
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

inline int ProfLevel = 0;

inline std::map<std::string, std::map<int, double>>& ProfRegistry()
{
	// Never destroyed: read by the atexit handler, which runs after the destruction of the function-local statics
	static auto* registry = new std::map<std::string, std::map<int, double>>();
	return *registry;
}

inline void ProfAdd(const std::string& name, std::chrono::steady_clock::time_point start)
{
	ProfRegistry()[name][ProfLevel] += std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

#define PROF_START(id) auto _prof_start_##id = std::chrono::steady_clock::now()
#define PROF_STOP(id, name) ProfAdd(name, _prof_start_##id)

inline void ProfPrint()
{
	int nLevels = 0;
	std::vector<std::pair<std::string, double>> totals;
	for (auto& [name, perLevel] : ProfRegistry())
	{
		double total = 0;
		for (auto& [level, seconds] : perLevel)
		{
			total += seconds;
			nLevels = std::max(nLevels, level + 1);
		}
		totals.push_back({ name, total });
	}
	std::sort(totals.begin(), totals.end(), [](auto& a, auto& b) { return a.second > b.second; });

	std::cout << "---- PROFILE (wall-clock seconds) ----" << std::endl;
	std::cout << std::left << std::setw(48) << "region" << std::right << std::setw(9) << "total";
	for (int l = 0; l < nLevels; l++)
		std::cout << std::setw(8) << ("L" + std::to_string(l));
	std::cout << std::endl;
	for (auto& [name, total] : totals)
	{
		std::cout << std::left << std::setw(48) << name << std::right << std::fixed << std::setprecision(3) << std::setw(9) << total;
		auto& perLevel = ProfRegistry()[name];
		for (int l = 0; l < nLevels; l++)
			std::cout << std::setw(8) << (perLevel.count(l) ? perLevel[l] : 0.0);
		std::cout << std::endl;
	}
}

inline int _profPrintAtExit = (std::atexit([] { ProfPrint(); }), 0);
