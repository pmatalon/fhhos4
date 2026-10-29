#!/bin/bash
# Scaling of the distributed U-AMG (MPI): processes x threads per process, the processes bound to their cores (MPICH's
# -bind-to core:<threads>), elapsed setup and solve times (best of 2) and iterations. Run nothing else meanwhile.
# Mind the memory: PERFORMANCE.md "Distributed memory (MPI)" gives the peak per process of the program.
#
# Usage (from build-mpi/, conda env activated): scaling_mpi.sh <output dir> <arguments of fhhos4...>
# CONFIGS: the configurations "processes,threads" (default: "1,1 1,8 2,1 2,4 4,1 4,2 8,1").
# Writes one log per run and summary.txt.
OUT=${1:?Usage: scaling_mpi.sh <output dir> <arguments of fhhos4...>}
shift
ARGS=("$@")
mkdir -p "$OUT"
: > "$OUT/summary.txt"
CONFIGS=${CONFIGS:-"1,1 1,8 2,1 2,4 4,1 4,2 8,1"}
for config in $CONFIGS; do
	np=${config%,*}; th=${config#*,}
	best_setup=999999; best_solve=999999
	for rep in 1 2; do
		log="$OUT/np${np}_t${th}_r${rep}.log"
		OMP_NUM_THREADS=$th timeout 1800 mpiexec -n $np -bind-to core:$th ./bin/fhhos4 "${ARGS[@]}" > "$log" 2>&1
		# The rows of the table after the solve: "Setup | CPU | Elapsed | ...", "Solving | ..."
		setup=$(awk -F'|' '/^Setup /{gsub(/ /,"",$3); split($3,t,":"); print t[1]*3600+t[2]*60+t[3]}' "$log")
		solve=$(awk -F'|' '/^Solving /{gsub(/ /,"",$3); split($3,t,":"); print t[1]*3600+t[2]*60+t[3]}' "$log")
		its=$(grep -m1 -oE "^[0-9]+ iterations" "$log")
		[ -n "$setup" ] && best_setup=$(echo "$setup $best_setup" | awk '{print ($1<$2)?$1:$2}')
		[ -n "$solve" ] && best_solve=$(echo "$solve $best_solve" | awk '{print ($1<$2)?$1:$2}')
	done
	echo "$np processes x $th threads: setup $best_setup s, solve $best_solve s, $its" | tee -a "$OUT/summary.txt"
done
