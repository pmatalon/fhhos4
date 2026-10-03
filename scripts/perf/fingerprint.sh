#!/bin/bash
# Iteration tables of the logs of a run directory (e.g. written by amg_papers.sh), reduced to
# (iteration, inner iterations, residual): the "fingerprint" of the results, without the times.
#
# Usage: fingerprint.sh <run dir>               prints the fingerprint
#        fingerprint.sh <run dir A> <run dir B> compares two runs
fingerprint() {
	for f in "$1"/*.log; do
		echo "## $(basename "$f")"
		grep -E '^ *[0-9]+ +(/ *[0-9]+ +)?[0-9.]+e[+-][0-9]+' "$f" | awk '{ if ($2 == "/") print $1, $3, $4; else print $1, "-", $2 }'
	done
}

if [ $# -eq 1 ]; then
	fingerprint "$1"
elif [ $# -eq 2 ]; then
	if diff <(fingerprint "$1") <(fingerprint "$2"); then
		echo "Identical iteration tables ($(ls "$1"/*.log | wc -l) runs)"
	else
		echo "The iteration tables differ (see above)"
		exit 1
	fi
else
	echo "Usage: fingerprint.sh <run dir> [<run dir to compare with>]"
	exit 2
fi
