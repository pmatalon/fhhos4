#!/bin/bash
# Makes the two files to copy to a cluster without internet, where scripts/cluster/run.pbs uses them:
#   fhhos4-env-<hash>.tar.gz  the conda env (conda-pack): copied once, then again only when <hash> changes (the packages
#                             of conda/environment.yml, the versions of the local env fhhos4, GLIBC)
#   fhhos4-src-<name>.tar.gz  the working tree (<name>: git describe), ~1 MB
# Usage (from the repository): scripts/cluster/pack.sh [output directory, default ~/code/fhhos4-cluster]
# GLIBC: the cluster's (ldd --version, default 2.34): no package of the env, nor the compiler's sysroot, needs a newer
# one. The env is also kept here (fhhos4-cluster, the cluster's packages): reproduce a failure of the cluster in it.
set -eu
out=$(realpath -m "${1:-$HOME/code/fhhos4-cluster}")
glibc=${GLIBC:-2.34}
env=fhhos4-cluster
cd "$(git rev-parse --show-toplevel)"
eval "$(conda shell.bash hook)"

# The packages of conda/environment.yml at the versions of the local env fhhos4 (those the release workflow pins), then
# conda-pack, and GNU time for the memory of the runs.
mapfile -t packages < <(sed -n '/^dependencies:/,$s/^  - //p' conda/environment.yml)
mapfile -t pins < <(conda list -n fhhos4 --export | sed -n 's/^\(gmsh\|eigen\|cgal-cpp\|gcc\|gxx\|mpich\|mumps-mpi\)=\([^=]*\)=.*/\1=\2/p')
spec=$(printf '%s\n' "${packages[@]}" "${pins[@]}" conda-pack time "glibc $glibc" | sort)
hash=$(echo "$spec" | sha1sum | cut -c1-8)
prefix=$(conda info --base)/envs/$env
if [ "$(cat "$prefix/fhhos4-cluster.spec" 2>/dev/null)" != "$spec" ]; then
	conda env remove -y -n $env > /dev/null 2>&1 || true
	CONDA_OVERRIDE_GLIBC=$glibc conda create -y -n $env -c conda-forge --override-channels "${packages[@]}" "${pins[@]}" \
		conda-pack time
	echo "$spec" > "$prefix/fhhos4-cluster.spec"
fi
mkdir -p "$out"
[ -f "$out/fhhos4-env-$hash.tar.gz" ] ||
	"$prefix/bin/conda-pack" -p "$prefix" -o "$out/fhhos4-env-$hash.tar.gz" --n-threads -1

# The tracked and the new files of the working tree, under src/<name>/, with the commit and the env's hash.
name=$(git describe --always --dirty)
tmp=$(mktemp -d)
mkdir -p "$tmp/scripts/cluster"
echo "$hash" > "$tmp/scripts/cluster/env-hash"
{ git rev-parse HEAD; echo "$name"; git status --short; } > "$tmp/COMMIT"
git ls-files -z --cached --others --exclude-standard |
	tar czf "$out/fhhos4-src-$name.tar.gz" --ignore-failed-read --transform "s|^|src/$name/|" \
		--null -T - -C "$tmp" COMMIT scripts/cluster/env-hash 2> >(grep -v 'Cannot stat' >&2)
rm -r "$tmp"

ls -lh "$out/fhhos4-env-$hash.tar.gz" "$out/fhhos4-src-$name.tar.gz"
cat << EOF

Copy to <base> on the cluster (a directory visible from the compute nodes):
  fhhos4-env-$hash.tar.gz   (unless <base> already has it)
  fhhos4-src-$name.tar.gz
Then, on the cluster:
  cd <base> && tar xzf fhhos4-src-$name.tar.gz && cd src/$name && qsub scripts/cluster/run.pbs
Bring back <base>/results/$name-<job>.tar.gz.
EOF
