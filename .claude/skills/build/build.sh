#!/usr/bin/env bash
# Builds one configured build directory with ninja, alone, under a memory watchdog; then, with --, runs its tests.
# Usage: build.sh <build directory> [ninja arguments: targets, -j N, -n] [-- [ctest arguments]]
set -u
dir=$(realpath "$1")
shift
ninja_arguments=()
while [ $# -gt 0 ] && [ "$1" != -- ]; do
    ninja_arguments+=("$1")
    shift
done
run_tests=$#
[ $# -gt 0 ] && shift
cache="$dir/CMakeCache.txt"
[ -f "$cache" ] || { echo "$dir is not configured: run cmake first (the build skill)"; exit 2; }

# Only this user's processes: on a shared machine (a cluster's login node), the others' are not ours to wait for or kill.
user=$(id -u)
if pgrep -u "$user" -a -x 'cc1plus|ninja|fhhos4.*'; then
    echo "Another heavy job is running (above): wait for it, one heavy job at a time."
    exit 2
fi

# The conda env is the one of the compiler that configured the directory.
compiler=$(sed -n 's/^CMAKE_CXX_COMPILER:[A-Z]*=//p' "$cache")
env=$(echo "$compiler" | sed -n 's|.*/envs/\([^/]*\)/bin/.*|\1|p')
if [ -n "$env" ]; then
    eval "$(conda shell.bash hook)" && conda activate "$env" || exit 2
else
    env="current (the compiler $compiler is in no conda env: modules or Spack must be loaded)"
fi
echo "$dir: env $env, $(sed -n 's/^CMAKE_BUILD_TYPE:[A-Z]*=//p' "$cache"), ON: $(sed -n 's/^ENABLE_\([A-Z0-9_]*\):BOOL=ON$/\1/p' "$cache" | tr '\n' ' ')"

# An interrupted build leaves empty objects newer than their sources, which ninja then takes as up to date.
find "$dir" -name '*.o' -empty -print -delete

available=$(awk '/MemAvailable/{print $2}' /proc/meminfo)
jobs=()
# The pool of 2 (Program.cpp + a test file) went from 6.9 GB available to below 1.5 GB.
[ "$available" -lt 8000000 ] && jobs=(-j1)

log="$dir/build.log"
ninja -C "$dir" "${jobs[@]}" "${ninja_arguments[@]}" > "$log" 2>&1 &
build=$!
fired=
while kill -0 "$build" 2>/dev/null; do
    available=$(awk '/MemAvailable/{print $2}' /proc/meminfo)
    if [ "$available" -lt 1500000 ]; then
        pkill -u "$user" -x cc1plus  # -x: the process name; -f would also match the Bash tool's shell running this script
        fired=$available
    fi
    sleep 2
done
wait "$build"
status=$?

grep -E 'error:|FAILED:' "$log" | head -40
tail -1 "$log"
[ -n "$fired" ] && echo "WATCHDOG: compilation killed at $fired kB available"
echo "build exit $status (${jobs[*]:-default pool}), full log: $log"
[ $status -ne 0 ] || [ "$run_tests" -eq 0 ] && exit $status

log="$dir/test.log"
ctest --test-dir "$dir" --output-on-failure "$@" > "$log" 2>&1
status=$?
sed -n '/tests passed/,$p' "$log"
echo "test exit $status, full log: $log"
exit $status
