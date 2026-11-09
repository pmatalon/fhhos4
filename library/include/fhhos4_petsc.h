// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Pierre Matalon
#pragma once

// fhhos4's preconditioner for PETSc, -pc_type fhhos4 (libfhhos4_petsc, built with -DENABLE_PETSC_PC=ON): see
// README_library.md. A PETSc program loads it without being recompiled with the option
// -dll_append <dir>/libfhhos4_petsc.so (PETSc then calls PetscDLLibraryRegister_fhhos4_petsc()); a program linked to
// it calls PetscDLLibraryRegister_fhhos4_petsc() once, after PetscInitialize().
#include <petscpc.h>

#define PCFHHOS4 "fhhos4"

// Registers PCFHHOS4 (PCRegister())
PETSC_EXTERN PetscErrorCode PetscDLLibraryRegister_fhhos4_petsc(void);
