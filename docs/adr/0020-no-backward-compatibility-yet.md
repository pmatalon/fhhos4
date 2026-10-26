---
Status: accepted, to revisit once the library has external users
---

# No backward compatibility in the library's API yet

The library has no external users yet (2026-10-14): a renamed header, CMake target, function or parameter is renamed
outright, without an alias or deprecated shim, and the client codes are updated instead, despite the SameMajorVersion
promise of the CMake package. Breaking changes still follow the version rules (a new public member of `Solver`: a new
minor version, since it breaks the ABI). Revisit once the conda-forge and Spack packages have users.
