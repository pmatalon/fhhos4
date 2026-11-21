#!/bin/bash
# save_before.sh <dir> [build dir, default: build]: saves the program and the library of a build directory (the
# version before a change) into <dir>, and writes <dir>/fhhos4.sh, which runs them: the library preloaded (the
# program's rpath names <build dir>/lib, where the next build puts the new library: a copied program alone would run
# with it) and the meshes of the sources (FHHOS4_DATA_DIR: the copy is no longer next to them). Prints the version of
# the saved program (its commit, "-dirty" if built from a modified tree).
set -e
DIR=${1:?Usage: save_before.sh <dir> [build dir, default: build]}
BUILD=$(readlink -f "${2:-build}")
ROOT=$(git -C "$BUILD" rev-parse --show-toplevel)
LIB=$(readelf -d "$BUILD/bin/fhhos4" | sed -n 's/.*NEEDED.*\[\(libfhhos4[^]]*\)\].*/\1/p')
[ -n "$LIB" ] || { echo "$BUILD/bin/fhhos4 does not link libfhhos4" >&2; exit 1; }
mkdir -p "$DIR"
DIR=$(readlink -f "$DIR")
cp -p "$BUILD/bin/fhhos4" "$DIR/fhhos4"
cp -pL "$BUILD/lib/$LIB" "$DIR/$LIB"
cat > "$DIR/fhhos4.sh" <<EOF
#!/bin/bash
# The version saved by save_before.sh from $BUILD: its program with its library
FHHOS4_DATA_DIR=$ROOT/data LD_PRELOAD=$DIR/$LIB exec $DIR/fhhos4 "\$@"
EOF
chmod +x "$DIR/fhhos4.sh"
"$DIR/fhhos4.sh" -h | grep -m1 '^Version'
