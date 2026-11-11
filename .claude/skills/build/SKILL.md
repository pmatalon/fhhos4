---
name: build
description: Builds fhhos4 in a build directory (an existing one or a new one) one job at a time under a memory watchdog, then optionally runs its tests. Use before any cmake or ninja run: a build, a rebuild, a compile check, a build and test.
---

# Build

`Program.cpp` alone takes most of a build's memory and time (`CLAUDE.md`, "Build"; more than 5.5 GB in Debug), and a
machine that runs out of memory may crash (WSL does). Every build therefore goes through `build.sh`: it refuses to
start beside another of your compilations or fhhos4 runs, activates the directory's conda env (a compiler outside
conda builds in the current environment: load the modules or the Spack env first), compiles with `-j1` under 8 GB
available, and kills the compilation below 1 GB.

## Steps

1. **Pick the directory** from the table below. The script prints the build type and options that are ON: the
   `CMakeCache.txt` is the source of truth.
2. **Configure** only a new directory, or an existing one whose option changes (seconds, in the foreground), with the
   env and options of the directory's document:
   ```bash
   eval "$(conda shell.bash hook)" && conda activate <env>
   cmake -S . -B <dir> -G Ninja -DCMAKE_BUILD_TYPE=Release <options>   # an existing directory: cmake -D<OPTION>=ON <dir>
   ```
   `conda activate` needs the `eval` line: the Bash tool's shell has no `conda init`. Without conda, load the modules
   or the Spack env instead of these two lines. An env that doesn't exist yet: `conda env create -f
   conda/environment.yml` (`fhhos4`); the others, and the envs' quirks: the directory's document.
3. **Build**, with `run_in_background: true` (a build with `Program.cpp` outlasts the Bash tool's 2-minute default
   timeout), then wait for the notification:
   ```bash
   .claude/skills/build/build.sh <dir> [targets]            # build only
   .claude/skills/build/build.sh <dir> [targets] -- [ctest arguments]   # then the tests, if the build succeeded
   ```
   Arguments after `--` go to `ctest` (e.g. `-R UAMGLibrary`, `-E HPConfig`), run in the directory's env:
   `ctest` exists only inside the envs. Several directories: chain the calls with `&&` in one command.
4. **Read the outcome**: the script prints the compiler errors, the last line of ninja, the exit status and the log
   (`<dir>/build.log`, `<dir>/test.log`); with tests, the ctest summary and the failed tests.
   - `Another heavy job is running`: someone else is compiling or running fhhos4 (another session, the user). Wait
     for it to end (Monitor on `pgrep -x 'cc1plus|ninja|fhhos4.*'`); leave it running.
   - `WATCHDOG`: memory ran short. Report `free -h` to the user and wait for memory to be freed: a retry with the
     same memory dies the same way.
   - A link error `file truncated`: an interrupted build left that object half-written; delete it and rebuild (the
     script already deletes the empty ones).

## Which directory

| Directory | For | Env and options |
|---|---|---|
| `build` | the program and `fhhos4_tests`; the final validation | `CLAUDE.md`, "Build" |
| `build-mpi-lib` | the library on several processes, `fhhos4_mpi_tests`: the fast loop of MPI work | `docs/mpi.md` |
| `build-mpi` | the program on several processes (`-- -R ProgramMPI`) | `docs/mpi.md` |
| `build-mpi-bigindex` | a change of index or state handling (`CLAUDE.md`, "Test"): `-Wtype-limits` warnings in `build.log` | `build-mpi`'s, `-DENABLE_64BIT_INDEX=ON`, configured with `CXXFLAGS="$CXXFLAGS -Wtype-limits"` |
| `build-python` | the Python module | `docs/python.md` |
| `build-python-mpi` | the Python module on several processes | `docs/python.md` |
| `build-python-openmpi` | NGSolve on several processes | `docs/python.md` |
| `build-python-petsc` | the preconditioner for PETSc | `docs/library.md` |

- **Debug while iterating** on light translation units: `build-debug` (env `fhhos4`, `-DCMAKE_BUILD_TYPE=Debug
  -DENABLE_PROGRAM=ON`, create it when needed), one object at a time (`ninja -C build-debug -t targets all | grep
  <File>.cpp.o` names the target). A change to a header that `Program.cpp` includes compiles in the Release `build`
  directly: its Debug compilation exceeds the memory.
- **Release** for timings, the iteration-count tests and the final validation.
- **Every tree compiles the library**: after a change to `library/`, to the solver code it shares, or to
  `CMakeLists.txt`, the other directories are stale. Rebuild those the change concerns before saying it works there,
  and name the ones left stale to the user.
