---
name: validate
description: Validates a change of fhhos4's code against the version before it - saves that version, iterates on the fastest build tree, then runs the final round once (solutions.sh bit for bit, ctest in every build tree the change reaches). Use before editing solver or library code (to save the version before), and once the change is settled, before saying it keeps the results.
---

# Validate

A change keeps the results when it keeps the papers' iteration counts (`PERFORMANCE.md`, "Rules"): an exact change
keeps the solutions to the last bit, a change that reorders operations keeps the iteration counts. Two checks, each
blind to what the other sees: `scripts/perf/solutions.sh` runs ~40 cases with two versions and compares their
solutions bit for bit (relative: no reference results); `ctest` checks the iteration counts against the papers'
(absolute: not the bits). The comparison needs the **before version**, saved before the change's first rebuild.

The **final round** (`solutions.sh` and every tree the change reaches, one after the other) takes 25-45 minutes: run it
once, on the settled change, and iterate on one fast tree before it.

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
2. **Iterate** in the fastest tree that compiles the change, with the tests it concerns (the `build` skill):
   `build-mpi-lib` (the library and its MPI tests, ~3 minutes) for library or distributed code, `build` with
   `-- -R <tests>` for the program. Done when the change is **settled**: its code, names and interface final. The new
   public names it brings (a parameter, an enum and its values, an option) that the user didn't fix: show them to the
   user now, since a rename after the final round runs it again.
3. **Pick the final checks** by the kind of change:

   | The change | `solutions.sh` must give |
   |---|---|
   | exact: an optimization, a refactoring, the memory | `OK` everywhere |
   | reorders arithmetic operations | `DIFF`, `RESIDUALS-DIFF` allowed; `ITERATIONS-DIFF` never |
   | a new option (a change of the method is opt-in) | `OK` everywhere: the defaults unchanged |
   | the same computations: renames, comments, documentation, tests only | not run: the trees' `ctest` suffices |

   The trees: `build` always; the program on several processes (its distributed code, anything under `ENABLE_MPI`)
   also `build-mpi`; code the library compiles (`library/`, the solver code it shares, `CMakeLists.txt`) all five,
   `build`, `build-mpi`, `build-mpi-lib`, `build-python`, `build-python-petsc`. The library's interface or behaviour:
   also the `clients` skill.
4. **Compare**, if step 3 asks: build `build` (the `build` skill), then from `build/` with the conda env activated, in
   the background (~12 minutes, a heavy job: alone, as in the `build` skill):
   ```bash
   ../scripts/perf/solutions.sh before-<change>/solutions before-<change>/fhhos4.sh ./bin/fhhos4
   ```
   Done when every run of `result.txt` has a status allowed by step 3's row or a known difference (the script's
   header): count them (`awk -F'|' '{print $3}' result.txt | sort | uniq -c`). An `ITERATIONS-DIFF` or an
   `EXIT-DIFF` means the change alters the behaviour: find why in the two runs' `log` before going on.
5. **The trees** of step 3, in one command, in the background:
   ```bash
   .claude/skills/validate/trees.sh build build-mpi build-mpi-lib build-python build-python-petsc
   ```
   It builds and tests each tree with the `build` skill's script, one after the other (the tests of each tree: its
   header), then prints one line per tree. Done when every line says `OK`. A test comparing two code paths bit for
   bit (`CoarseningPassTest`: the flat passes against the mesh passes) fails on a reordering made in one path only:
   make it in both.
6. **Report**: the statuses counted, the line of each tree, and the trees the change reaches left stale.
