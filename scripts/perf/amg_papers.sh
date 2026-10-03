#!/bin/bash
# Runs the AMG configurations of the papers at reduced sizes (the paper sizes don't fit in 13 GB of RAM):
# reproducibility/2023_AMG_for_hybrid_methods.md (U-AMG and C-AMG) and the U-AMG runs of
# reproducibility/2023_biharmonic_problem.md. 30 runs, ~6 minutes.
# Compare the iteration tables of two runs with fingerprint.sh: a change that keeps the results must give identical ones.
#
# Usage (from build/, conda env activated): amg_papers.sh <output dir> [fhhos4 binary, default: ./bin/fhhos4]
# Writes one log per run, and summary.txt: label | exit code | iterations | setup time | solve time (elapsed).
OUT=${1:?Usage: amg_papers.sh <output dir> [fhhos4 binary]}
BIN=${2:-./bin/fhhos4}
mkdir -p "$OUT"
: > "$OUT/summary.txt"

run() {
	local label="$1"; shift
	local log="$OUT/$(echo "$label" | tr -c 'A-Za-z0-9._\n-' '_').log"
	timeout 600 "$BIN" "$@" -no-cache > "$log" 2>&1
	local code=$?
	local it=$(grep -E '^[0-9]+ iterations' "$log" | tail -1 | awk '{print $1}')
	local setup=$(grep '^Setup ' "$log" | tail -1 | awk -F'|' '{gsub(/ /,"",$3); print $3}')
	local solve=$(grep '^Solving |' "$log" | tail -1 | awk -F'|' '{gsub(/ /,"",$3); print $3}')
	printf "%-58s | %3s | %4s | %s | %s\n" "$label" "$code" "${it:--}" "${setup:--}" "${solve:--}" | tee -a "$OUT/summary.txt"
}

# Figure 4.2 / Table 4.2: test cases x solvers (paper size in parentheses), and Figure 4.4
for s in fcguamg fcgaggregamg; do
	run "Cube-cart n32 (128) $s"          -geo cube -mesh cart -mesher inhouse -k 0 -n 32 -s $s
	run "Cube-tet n16 (64) $s"            -geo cube -mesh tetra -k 0 -n 16 -s $s
	run "Complex-tet n4 (16) $s"          -geo platewith4holes -k 0 -n 4 -s $s
	run "Heterog1e8 n256 (2048) $s"       -geo square4quadrants -k 0 -n 256 -heterog 1e8 -s $s
	run "Cube-cart-aniso100 n32 (128) $s" -geo cube -mesh cart -mesher inhouse -k 0 -n 32 -aniso 100 -s $s
	run "Cube-tet-aniso20 n16 (64) $s"    -geo cube -mesh tetra -k 0 -n 16 -aniso 20 -s $s
	run "Fig4.4 Cube-tet n32 $s"          -geo cube -mesh tetra -k 0 -n 32 -s $s
done

# Section 4.3.3 / Table 4.3: coarsening strategies
for cs in mpa dpa; do
	run "Cube-tet n16 fcguamg -cs $cs" -geo cube -mesh tetra -k 0 -n 16 -s fcguamg -cs $cs
done

# Table 4.4: coarsening prolongations
for cp in 3 4 5 6; do
	run "Cube-tet n16 fcguamg -coarsening-prolong $cp"           -geo cube -mesh tetra -k 0 -n 16 -s fcguamg -coarsening-prolong $cp
	run "Cube-cart-aniso100 n32 fcguamg -coarsening-prolong $cp" -geo cube -mesh cart -mesher inhouse -k 0 -n 32 -aniso 100 -s fcguamg -coarsening-prolong $cp
done

# Biharmonic paper, Tables 1-3 (U-AMG for the Laplacian problems)
for n in 8 16; do
	for p in s no; do
		run "Bihar cube-tet n$n -bihar-prec $p fcguamg" -pb bihar -geo cube -source exp -mesh tetra -not-compute-errors -s fcguamg -hp-cs p_h -nbh-depth 2 -bihar-prec-solver bicgstab -tol 1e-8 -k 0 -bihar-prec $p -n $n
	done
done
# Biharmonic paper, commented-out "Heuristics" (k=3)
for o in 0 2; do
	run "Bihar square-tri k3 n32 -opt2 $o fcguamg" -pb bihar -geo square -source poly -s fcguamg -mesh tri -cs r -k 3 -n 32 -tol 1e-8 -bihar-prec p -nbh-depth 2 -opt2 $o
done
