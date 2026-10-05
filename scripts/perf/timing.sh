#!/bin/bash
# Setup and solve times of U-AMG on the reference cases of PERFORMANCE.md, sequential (-threads 1) and parallel
# (-threads 0: all logical cores). Run nothing else meanwhile (no build): the noise is already about +-10%.
#
# Usage (from build/, conda env activated): timing.sh <output dir> [fhhos4 binary, default: ./bin/fhhos4] [repetitions, default: 2]
# Writes one log per run, timing.txt (one line per run) and best.txt (best elapsed setup and solve times of each case).
OUT=${1:?Usage: timing.sh <output dir> [fhhos4 binary] [repetitions]}
BIN=${2:-./bin/fhhos4}
REPS=${3:-2}
mkdir -p "$OUT"
: > "$OUT/timing.txt"

seconds() { # hh:mm:ss.mmm -> seconds
	echo "$1" | awk -F: '{ print $1 * 3600 + $2 * 60 + $3 }'
}

run() {
	local label="$1"; shift
	for threads in 1 0; do
		for r in $(seq 1 $REPS); do
			local log="$OUT/$(echo "$label" | tr -c 'A-Za-z0-9._\n-' '_')_t${threads}_r${r}.log"
			"$BIN" "$@" -threads $threads -no-cache > "$log" 2>&1
			local it=$(grep -E '^[0-9]+ iterations' "$log" | tail -1 | awk '{print $1}')
			local setup=$(grep '^Setup ' "$log" | tail -1 | awk -F'|' '{gsub(/ /,"",$3); print $3}')
			local solve=$(grep '^Solving |' "$log" | tail -1 | awk -F'|' '{gsub(/ /,"",$3); print $3}')
			printf "%-26s | %-10s | %3s it | setup %6.2f s | solve %6.2f s\n" "$label" "$([ $threads = 1 ] && echo sequential || echo parallel)" "${it:--}" \
				"$(seconds "${setup:-0:0:0}")" "$(seconds "${solve:-0:0:0}")" | tee -a "$OUT/timing.txt"
		done
	done
}

run "Cube-tet k0 n32 cp6"      -geo cube -mesh tetra -k 0 -n 32 -s fcguamg -coarsening-prolong 6
run "Cube-tet k0 n32 cp4"      -geo cube -mesh tetra -k 0 -n 32 -s fcguamg -coarsening-prolong 4
run "Cube-tet k0 n32 cp5"      -geo cube -mesh tetra -k 0 -n 32 -s fcguamg -coarsening-prolong 5
run "Cube-tet k0 n32 cp3"      -geo cube -mesh tetra -k 0 -n 32 -s fcguamg -coarsening-prolong 3
run "Cube-cart-aniso100 n64"   -geo cube -mesh cart -mesher inhouse -k 0 -n 64 -aniso 100 -s fcguamg
# -hp-cs h: the hp-coarsening of the measurements recorded in PERFORMANCE.md (the default is p_h since 2026-10-05)
run "Cube-tet k2 n16"          -geo cube -mesh tetra -k 2 -n 16 -s fcguamg -hp-cs h

# Best of the repetitions, per case and number of threads
awk -F'|' '{ key = $1 "|" $2; split($4, s, " "); split($5, t, " ");
             if (!(key in setup) || s[2] < setup[key]) setup[key] = s[2];
             if (!(key in solve) || t[2] < solve[key]) solve[key] = t[2];
             if (!(key in seen)) { seen[key] = 1; order[n++] = key } }
     END { for (i = 0; i < n; i++) printf "%s| best setup %6.2f s | best solve %6.2f s\n", order[i], setup[order[i]], solve[order[i]] }' \
	"$OUT/timing.txt" | tee "$OUT/best.txt"
