#!/bin/bash
# Compiles uamg_harness.cpp with the compiler, flags and definitions of the build in build/ (~30 s, ~1.5 GB of RAM).
# Usage (from build/, conda env activated): build_uamg_harness.sh <output binary> [<source root, default: ..>]
#   The source root lets compile the same harness against another version of src/ (e.g. a `git archive` of HEAD,
#   to compare before/after).
OUT=${1:?Usage: build_uamg_harness.sh <output binary> [source root]}
ROOT=$(realpath "${2:-..}")
HERE=$(dirname "$(realpath "$0")")
CXX=$(grep '^CMAKE_CXX_COMPILER:' CMakeCache.txt | cut -d= -f2)
BLOCK=$(grep -A6 'build CMakeFiles/fhhos4_core.dir/src/Program.cpp.o' build.ninja)
FLAGS=$(echo "$BLOCK" | grep '^  FLAGS = ' | sed 's/^  FLAGS = //')
DEFINES=$(echo "$BLOCK" | grep '^  DEFINES = ' | sed 's/^  DEFINES = //')
INCLUDES=$(echo "$BLOCK" | grep '^  INCLUDES = ' | sed 's/^  INCLUDES = //' | sed "s|-I[^ ]*/src |-I$ROOT/src |")
eval "$CXX $FLAGS $DEFINES $INCLUDES -I$HERE $HERE/uamg_harness.cpp -o $OUT -fopenmp"
