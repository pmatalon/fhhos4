#!/bin/bash
# Bitwise comparison of the results of two binaries: the solution vectors (exported with 17 significant digits, i.e.
# exactly; -0 and 0 are considered equal) and the iteration tables (without the times), sequential (-threads 1) and
# parallel (-threads 0). Stronger than fingerprint.sh (3-digit residuals): checks that a change of the solve phase is
# bit-identical (status OK). For a change that reorders arithmetic operations, DIFF and RESIDUALS-DIFF are expected,
# but not ITERATIONS-DIFF (different numbers of iterations or inner iterations).
# The cases cover the code paths of the smoothers: scalar and block (sizes 3, 4, 6, 7, 8, 10) Gauss-Seidel, all the
# directions, from x = 0 or not, with exact zeros in the right-hand side (-source one), standalone and in multigrids;
# block SOR (also used by the standalone solvers bgs and rbgs).
# 35 cases x 2 thread counts x 2 binaries, ~10 minutes.
#
# Known differences: the hybrid smoothers (hbgs, hrbgs) in parallel vary from one run to the next (data race),
# and 'fcgaggregamg' fails at k=1 (exit code 1, no vectors).
#
# Usage (from build/, conda env activated): solutions.sh <output dir> <binary before> [binary after, default: ./bin/fhhos4]
# Writes one directory per run and result.txt: case | threads | status | iteration lines | exit code | arguments.
OUT=${1:?Usage: solutions.sh <output dir> <binary before> [binary after]}
BEFORE=${2:?Usage: solutions.sh <output dir> <binary before> [binary after]}
AFTER=${3:-./bin/fhhos4}
mkdir -p "$OUT"
: > "$OUT/result.txt"
i=0

# Whether two exported vectors have the same values (-0 and 0 being equal)
same_values() {
	local normalize='{ if (NF == 1 && $1 == 0) print 0; else print }'
	cmp -s <(awk "$normalize" "$1") <(awk "$normalize" "$2")
}

# Iteration numbers and inner iteration counts of an iteration table, without the residuals
iterations() {
	awk '{ if ($2 == "/") print $1, $3; else print $1 }' "$1"
}

check() {
	i=$((i+1))
	for threads in 1 0; do
		for tag in before after; do
			local bin=$([ $tag = before ] && echo "$BEFORE" || echo "$AFTER")
			local d="$OUT/c${i}_t${threads}_$tag"
			rm -rf "$d"; mkdir -p "$d"
			timeout 600 "$bin" "$@" -threads $threads -no-cache -export solvect -o "$d" > "$d/log" 2>&1
			echo $? > "$d/exit"
			# iteration tables without the times: the last column (remaining time) and the CPU time (2 columns before)
			grep -E '^ *[0-9]+ +(/ *[0-9-]+ +)?[0-9.]+e[+-][0-9]+' "$d/log" \
				| awk '{ if ($NF ~ /^([0-9]+:[0-9]+:[0-9]+|-)$/) { $NF = ""; $(NF-2) = "" } print }' > "$d/iters"
		done
		local b="$OUT/c${i}_t${threads}_before" a="$OUT/c${i}_t${threads}_after"
		local status=OK
		[ "$(ls "$b"/*.dat 2>/dev/null | wc -l)" -eq 0 ] && status="NO-VECTORS"
		for f in "$b"/*.dat; do
			[ -e "$f" ] || continue
			same_values "$f" "$a/$(basename "$f")" || status="DIFF($(basename "$f"))"
		done
		if ! cmp -s "$b/iters" "$a/iters"; then
			if cmp -s <(iterations "$b/iters") <(iterations "$a/iters"); then
				status="$status RESIDUALS-DIFF"
			else
				status="$status ITERATIONS-DIFF"
			fi
		fi
		cmp -s "$b/exit" "$a/exit" || status="$status EXIT-DIFF"
		printf "%-3s | t%s | %-28s | %4s it | exit %s | %s\n" "$i" "$threads" "$status" "$(wc -l < "$b/iters")" "$(cat "$b/exit")" "$*" | tee -a "$OUT/result.txt"
	done
}

# U-AMG, default smoothers (bgs/rbgs), per face block size
check -geo cube -mesh tetra -k 0 -n 16 -s fcguamg
check -geo cube -mesh tetra -k 1 -n 8 -s fcguamg
check -geo cube -mesh tetra -k 2 -n 8 -s fcguamg
check -geo cube -mesh tetra -k 3 -n 4 -s fcguamg
check -geo square -mesh tri -k 2 -n 32 -s fcguamg
check -geo square -mesh tri -k 6 -n 8 -s fcguamg
check -geo cube -mesh cart -mesher inhouse -k 1 -n 8 -aniso 100 -s fcguamg
# Other smoothers and cycles
check -geo cube -mesh tetra -k 2 -n 8 -s fcguamg -smoothers hbgs,hrbgs
check -geo cube -mesh tetra -k 0 -n 16 -s fcguamg -smoothers sgs
check -geo cube -mesh tetra -k 1 -n 8 -s fcguamg -smoothers sbgs
check -geo cube -mesh tetra -k 0 -n 16 -s fcguamg -smoothers ags
check -geo cube -mesh tetra -k 1 -n 8 -s fcguamg -smoothers abgs
check -geo cube -mesh tetra -k 1 -n 8 -s fcguamg -smoothers rbgs,bgs
check -geo cube -mesh tetra -k 2 -n 8 -s fcguamg -smoothers rbgs,bgs
check -geo cube -mesh tetra -k 0 -n 16 -s fcguamg -cycle V,2,2
check -geo cube -mesh tetra -k 1 -n 8 -s fcguamg -cycle W,2,1
check -geo cube -mesh tetra -k 0 -n 16 -s uamg
check -geo cube -mesh tetra -k 2 -n 8 -s fcguamg -smoothers bsor,rbsor -relax 1.2
check -geo cube -mesh tetra -k 1 -n 8 -s fcguamg -smoothers sbsor,sbsor -relax 0.9
# Exact zeros in the right-hand side (-source one: 2D only)
check -geo square -mesh tri -k 0 -n 32 -s fcguamg -source one
check -geo square -mesh tri -k 0 -n 32 -s fcguamg -source one -smoothers rgs,gs
check -geo square -mesh tri -k 2 -n 16 -s fcguamg -source one
check -geo square -mesh tri -k 3 -n 16 -s fcguamg -source one -smoothers rbgs,bgs
check -geo square -mesh tri -k 5 -n 8 -s fcguamg -source one -smoothers rbgs,bgs
check -geo square -mesh tri -k 7 -n 4 -s fcguamg -source one
# Other multigrids
check -geo cube -mesh tetra -k 0 -n 16 -s fcgaggregamg
check -geo cube -mesh tetra -k 1 -n 8 -s fcgaggregamg
check -geo square -mesh cart -mesher inhouse -k 1 -n 32 -s mg
check -geo square -mesh cart -mesher inhouse -k 0 -n 32 -s mg -cycle V,1,1
# Standalone smoothers
check -geo square -mesh cart -mesher inhouse -k 0 -n 16 -s gs
check -geo square -mesh cart -mesher inhouse -k 0 -n 16 -s rgs
check -geo square -mesh cart -mesher inhouse -k 0 -n 16 -s sgs
check -geo square -mesh cart -mesher inhouse -k 0 -n 16 -s sgs -source one
check -geo square -mesh cart -mesher inhouse -k 1 -n 16 -s bgs
check -geo square -mesh cart -mesher inhouse -k 1 -n 16 -s bsor -relax 1.3
