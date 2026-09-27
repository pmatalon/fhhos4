#pragma once

// Small out-param populated by Program_*::Execute() so callers (tests, in particular) can
// read the outcome of a run without scraping stdout. Left untouched (nullptr) by the normal
// CLI path.
struct ProgramResults
{
	int IterationCount = -1;
	double L2Error = -1;
};
