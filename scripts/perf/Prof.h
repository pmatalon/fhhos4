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
// The times are accumulated per region and level, and printed at exit, sorted by total (wall-clock seconds), with the
// minor page faults of the region (first writes of fresh memory): those of the process ("faults") and of the calling
// thread ("master"). In a parallel loop, the calling thread takes its share of the faults (1/threads); beyond it, the
// region writes fresh memory sequentially (e.g. a vector<T>(n) zero-filled by the calling thread).
// PROF_STOP_AT(id, name, level): PROF_STOP with the level set first; ProfScope{"name"}: a region until the end of the
// scope.
// Only time regions executed by one thread at a time (not inside a parallel loop).
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <sys/resource.h>

inline int ProfLevel = 0;

inline std::map<std::string, std::map<int, double>>& ProfRegistry()
{
	// Never destroyed: read by the atexit handler, which runs after the destruction of the function-local statics
	static auto* registry = new std::map<std::string, std::map<int, double>>();
	return *registry;
}

// Minor page faults of the process (first touches of fresh memory), per region
inline long ProfFaults() { rusage u; getrusage(RUSAGE_SELF, &u); return u.ru_minflt; }
inline long ProfThreadFaults() { rusage u; getrusage(RUSAGE_THREAD, &u); return u.ru_minflt; }
inline std::map<std::string, long>& ProfThreadFaultRegistry()
{
	static auto* registry = new std::map<std::string, long>();
	return *registry;
}
inline std::map<std::string, long>& ProfFaultRegistry()
{
	static auto* registry = new std::map<std::string, long>();
	return *registry;
}
struct ProfStartPoint { std::chrono::steady_clock::time_point T; long F; long FT = ProfThreadFaults(); };

inline void ProfAdd(const std::string& name, ProfStartPoint start)
{
	ProfRegistry()[name][ProfLevel] += std::chrono::duration<double>(std::chrono::steady_clock::now() - start.T).count();
	ProfFaultRegistry()[name] += ProfFaults() - start.F;
	ProfThreadFaultRegistry()[name] += ProfThreadFaults() - start.FT;
}

#define PROF_START(id) ProfStartPoint _prof_start_##id{ std::chrono::steady_clock::now(), ProfFaults() }
#define PROF_STOP(id, name) ProfAdd(name, _prof_start_##id)
#define PROF_STOP_AT(id, name, lvl) do { ProfLevel = (lvl); PROF_STOP(id, name); } while (0)

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
	std::cout << std::left << std::setw(48) << "region" << std::right << std::setw(9) << "total" << std::setw(10) << "faults" << std::setw(10) << "master";
	for (int l = 0; l < nLevels; l++)
		std::cout << std::setw(8) << ("L" + std::to_string(l));
	std::cout << std::endl;
	for (auto& [name, total] : totals)
	{
		std::cout << std::left << std::setw(48) << name << std::right << std::fixed << std::setprecision(3) << std::setw(9) << total << std::setw(10) << ProfFaultRegistry()[name] << std::setw(10) << ProfThreadFaultRegistry()[name];
		auto& perLevel = ProfRegistry()[name];
		for (int l = 0; l < nLevels; l++)
			std::cout << std::setw(8) << (perLevel.count(l) ? perLevel[l] : 0.0);
		std::cout << std::endl;
	}
}

inline int _profPrintAtExit = (std::atexit([] { ProfPrint(); }), 0);

// RAII: adds the time from its construction to its destruction (after the destructors of the later locals)
struct ProfScope
{
	std::string Name;
	ProfStartPoint Start{ std::chrono::steady_clock::now(), ProfFaults() };
	~ProfScope() { ProfAdd(Name, Start); }
};
