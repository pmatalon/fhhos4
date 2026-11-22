// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Pierre Matalon
#pragma once

// fhhos4's preconditioner for PETSc, -pc_type fhhos4 (libfhhos4_petsc, built with -DENABLE_PETSC_PC=ON): see
// library/doc/petsc.md. A PETSc program loads it without being recompiled with the option
// -dll_append <dir>/libfhhos4_petsc.so (PETSc then calls PetscDLLibraryRegister_fhhos4_petsc()); a program linked to
// it calls PetscDLLibraryRegister_fhhos4_petsc() once, after PetscInitialize().
#include <petscpc.h>

#define PCFHHOS4 "fhhos4"

// Registers PCFHHOS4 (PCRegister())
PETSC_EXTERN PetscErrorCode PetscDLLibraryRegister_fhhos4_petsc(void);

// The operator as one AIJ matrix in your numbering (instead of a MATNEST of the blocks): the rows of its cells and of
// its faces, before the setup. Each process gives rows it owns, of whole cells and faces, their DoFs in the library's
// order (the blocks of each cell and face contiguous, the components interleaved); the setup extracts the blocks
// (MatCreateSubMatrix()). The rows in neither index set (e.g. Dirichlet DoFs eliminated by MatZeroRowsColumns()) get
// the identity. Composes the index sets on the PC as "fhhos4_cells" and "fhhos4_faces" (petsc4py:
// pc.compose("fhhos4_cells", cells)).
PETSC_EXTERN PetscErrorCode PCFHHOS4SetIS(PC pc, IS cells, IS faces);
