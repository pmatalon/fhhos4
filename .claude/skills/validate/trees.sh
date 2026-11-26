#!/bin/bash
# trees.sh <tree>...: the final round of a validation, in each build tree given, one after the other (one heavy job at
# a time): the build skill's build.sh, then the tests of the tree, then one summary line per tree. The tests: in build,
# all but the 8 slowest (-E HPConfig: run them too before a release); in build-mpi, the program on several processes
# (-R ProgramMPI: the rest is build's and build-mpi-lib's); in the other trees, all theirs.
# Usage (from anywhere in the repository): .claude/skills/validate/trees.sh build build-mpi build-mpi-lib ...
# Exit status: 0 if every tree built and passed its tests.
[ $# -gt 0 ] || { echo "Usage: trees.sh <build tree>..." >&2; exit 2; }
cd "$(git rev-parse --show-toplevel)" || exit 2
summary=()
failed=0
for tree in "$@"; do
	case "$tree" in
		build) tests=(-E HPConfig) ;;
		build-mpi) tests=(-R ProgramMPI) ;;
		*) tests=() ;;
	esac
	output=$(.claude/skills/build/build.sh "$tree" -- "${tests[@]}" 2>&1)
	status=$?
	echo "$output"
	echo
	result=$(echo "$output" | grep -E -m1 "tests passed|error:|WATCHDOG|Another heavy job|not configured")
	[ $status -eq 0 ] || failed=1
	summary+=("$tree (${tests[*]:-all tests}): $([ $status -eq 0 ] && echo OK || echo FAILED), ${result:-no test summary}")
done
printf '%s\n' "${summary[@]}"
exit $failed
