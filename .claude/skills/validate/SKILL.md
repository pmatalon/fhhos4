---
name: validate
description: Validates a change of fhhos4's code against the version before it - saves that version, compares the solutions bit for bit (solutions.sh), runs ctest and the other build trees the change reaches. Use before editing solver or library code (to save the version before) and after, before saying a change keeps the results.
---

# Validate

A change keeps the results when it keeps the papers' iteration counts (`PERFORMANCE.md`, "Rules"): an exact change
keeps the solutions to the last bit, a change that reorders operations keeps the iteration counts. Two checks, each
blind to what the other sees: `scripts/perf/solutions.sh` runs ~40 cases with two versions and compares their
solutions bit for bit (relative: no reference results); `ctest` checks the iteration counts against the papers'
(absolute: not the bits). The comparison needs the **before version**, saved before the change's first rebuild.

## Steps

1. **Save the before version**, before rebuilding the change (from the repository root):
   ```bash
   .claude/skills/validate/save_before.sh build/before-<change>
   ```
   It copies the program and its library and writes `build/before-<change>/fhhos4.sh`, which runs them (why a
   wrapper: the script's header). Done when the version it prints is the commit before the change (`-dirty`: built
   from a modified tree, which is then the before version). Already rebuilt: build the commit before in a directory
   of its own (`git archive <commit> | tar -x -C build/<name>-src`, configured with `-S build/<name>-src`, the `build`
   skill), whose program runs as it is.
2. **Pick the checks** by the kind of change:

   | The change | `solutions.sh` must give |
   |---|---|
   | exact: an optimization, a refactoring, the memory | `OK` everywhere |
   | reorders arithmetic operations | `DIFF`, `RESIDUALS-DIFF` allowed; `ITERATIONS-DIFF` never |
   | a new option (a change of the method is opt-in) | `OK` everywhere: the defaults unchanged |

   Always `ctest` (step 5). Code compiled on several processes too (`DistributedMatrix`, the distributed passes,
   anything under `ENABLE_MPI`): `build-mpi-lib` (all its tests) and `build-mpi` (`-R ProgramMPI`). The library's
   interface or behaviour: the `clients` skill.
3. **Build** the change with the `build` skill (`build`).
4. **Compare**, from `build/` with the conda env activated, in the background (~12 minutes, a heavy job: alone, as
   in the `build` skill):
   ```bash
   ../scripts/perf/solutions.sh before-<change>/solutions before-<change>/fhhos4.sh ./bin/fhhos4
   ```
   Done when every run of `result.txt` has a status allowed by step 2's row or a known difference (the script's
   header): count them (`awk -F'|' '{print $3}' result.txt | sort | uniq -c`). An `ITERATIONS-DIFF` or an
   `EXIT-DIFF` means the change alters the behaviour: find why in the two runs' `log` before going on.
5. **`ctest`** with the `build` skill: `.claude/skills/build/build.sh build -- -E HPConfig` (the 8 slowest skipped:
   run them too before a release). Done at 0 failures. A test comparing two code paths bit for bit
   (`CoarseningPassTest`: the flat passes against the mesh passes) fails on a reordering made in one path only: make
   it in both.
6. **The other trees** of step 2, with the `build` skill.
7. **Report**: the statuses counted, the `ctest` result, the trees run, and those the change reaches left stale.
