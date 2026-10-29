#pragma once
// Temporary measurement of the memory of each process (resident set) between points of the code. Not part of the build:
// include it in the files to instrument, then revert them, e.g.
//
//     #include "../../scripts/perf/Mem.h"   // path relative to the instrumented file
//     MemPoint("assembly");
//
// Prints on stderr (all the processes of mpiexec; the rank from MPICH's PMI_RANK) the resident memory now and its peak
// since the previous point (VmRSS and VmHWM of /proc/self/status, then VmHWM reset by /proc/self/clear_refs), in MB.
// The peak of the process is the largest of its points: put one at the end. Call it outside the parallel regions.
#include <cstdio>
#include <cstdlib>

inline void MemPoint(const char* label)
{
	long rss = 0, hwm = 0;
	char line[256];
	if (FILE* f = std::fopen("/proc/self/status", "r"))
	{
		while (std::fgets(line, sizeof line, f))
		{
			std::sscanf(line, "VmRSS: %ld", &rss);
			std::sscanf(line, "VmHWM: %ld", &hwm);
		}
		std::fclose(f);
	}
	if (FILE* f = std::fopen("/proc/self/clear_refs", "w"))
	{
		std::fputs("5", f);
		std::fclose(f);
	}
	const char* rank = std::getenv("PMI_RANK");
	std::fprintf(stderr, "MEM rank %s | %-26s | rss %8.1f MB | peak since previous %8.1f MB\n", rank ? rank : "0", label, rss / 1024.0, hwm / 1024.0);
}
