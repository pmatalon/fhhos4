#!/bin/bash
# Puts what to bring back of jobs of scripts/cluster/run.pbs into one archive: per job, its results directory, its PBS
# output (src/<name>/<job name>.o<job>, written by PBS when the job ends) and qstat -fx (resources used, exit status).
# A running job gives what it has written so far.
# Usage, on the login node: cd <base> && bash src/<name>/scripts/cluster/collect.sh <job>...
# (<job>: the id qsub printed, with or without its server suffix). Writes <base>/fhhos4-results-<job>[-<job>...].tar.gz.
set -eu
shopt -s nullglob
[ $# -gt 0 ] || { echo "Usage: cd <base> && bash src/<name>/scripts/cluster/collect.sh <job>..."; exit 1; }
[ -d results ] || { echo "$(pwd): no results/ here: run from <base>"; exit 1; }
jobs=()
files=()
for job in "$@"; do
	job=${job%%.*}
	jobs+=("$job")
	found=(results/*-"$job" src/*/*.[oe]"$job")
	[ ${#found[@]} -gt 0 ] || echo "$job: neither results/*-$job nor src/*/*.o$job in $(pwd)"
	qstat -fx "$job" > "results/qstat-$job.txt" 2>&1 || echo "$job: qstat -fx failed (results/qstat-$job.txt)"
	files+=("${found[@]}" "results/qstat-$job.txt")
done
archive=$(pwd -P)/fhhos4-results-$(IFS=-; echo "${jobs[*]}").tar.gz
tar czf "$archive" "${files[@]}"
ls -lh "$archive"
