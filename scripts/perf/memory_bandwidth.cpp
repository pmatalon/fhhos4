// The memory limits of a node for U-AMG's access patterns, on 1, 2, 4... threads up to OMP_NUM_THREADS: the first write
// of fresh memory (its page faults), a streaming read (sum of a 512 MB array) and a random gather (sum of a[idx[i]],
// one cache miss per element), best of 5. Argument "serial": the read arrays first written by one thread (all their
// pages on its NUMA node, as an array zero-filled by the calling thread); default: by all the threads (each page on the
// NUMA node of the thread that wrote it, if the threads are bound to their cores: OMP_PLACES=cores OMP_PROC_BIND=close).
// Build: $CXX -O3 -fopenmp memory_bandwidth.cpp (scripts/cluster/run.pbs, step memory).
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <omp.h>

static double SecondsSince(std::chrono::steady_clock::time_point start)
{
	return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

int main(int argc, char** argv)
{
	bool serialPlacement = argc > 1 && std::strcmp(argv[1], "serial") == 0;
	size_t n = (size_t)64 << 20; // 64M doubles = 512 MB
	size_t nGathers = n / 4;
	size_t bytes = n * sizeof(double);
	// malloc: untouched pages, placed by their first write below (a std::vector would zero-fill them on this thread)
	double* a = (double*)std::malloc(bytes);
	uint32_t* idx = (uint32_t*)std::malloc(nGathers * sizeof(uint32_t));
	#pragma omp parallel for schedule(static) if(!serialPlacement)
	for (size_t i = 0; i < n; i++)
		a[i] = (double)(i % 1000);
	#pragma omp parallel for schedule(static) if(!serialPlacement)
	for (size_t i = 0; i < nGathers; i++)
	{
		uint64_t z = i * 0x9E3779B97F4A7C15ull;
		z ^= z >> 29;
		idx[i] = (uint32_t)(z % n);
	}

	int maxThreads = omp_get_max_threads();
	volatile double sink = 0;
	printf("placement %s, %d threads at most, OMP_PROC_BIND=%s\n", serialPlacement ? "serial" : "parallel", maxThreads,
	       std::getenv("OMP_PROC_BIND") ? std::getenv("OMP_PROC_BIND") : "(unset)");
	for (int t = 1;; t = std::min(2 * t, maxThreads))
	{
		double firstWrite = 1e9, stream = 1e9, gather = 1e9;
		for (int rep = 0; rep < 5; rep++)
		{
			// Above glibc's mmap threshold: fresh pages at each malloc, unmapped by free
			double* fresh = (double*)std::malloc(bytes);
			auto start = std::chrono::steady_clock::now();
			#pragma omp parallel for num_threads(t) schedule(static)
			for (size_t i = 0; i < n; i++)
				fresh[i] = (double)i;
			firstWrite = std::min(firstWrite, SecondsSince(start));
			sink = sink + fresh[n / 3];
			std::free(fresh);

			double sum = 0;
			start = std::chrono::steady_clock::now();
			#pragma omp parallel for num_threads(t) reduction(+:sum) schedule(static)
			for (size_t i = 0; i < n; i++)
				sum += a[i];
			stream = std::min(stream, SecondsSince(start));

			double gathered = 0;
			start = std::chrono::steady_clock::now();
			#pragma omp parallel for num_threads(t) reduction(+:gathered) schedule(static)
			for (size_t i = 0; i < nGathers; i++)
				gathered += a[idx[i]];
			gather = std::min(gather, SecondsSince(start));
			sink = sink + sum + gathered;
		}
		printf("threads %d: first write %.1f GB/s, stream %.1f GB/s, random gather %.0f M/s\n", t, bytes / firstWrite * 1e-9,
		       bytes / stream * 1e-9, nGathers / gather * 1e-6);
		if (t == maxThreads)
			break;
	}
	std::free(a);
	std::free(idx);
}
