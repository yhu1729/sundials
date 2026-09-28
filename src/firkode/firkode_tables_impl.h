/* -----------------------------------------------------------------
 * Programmer(s): Yifan Hu @ UMBC
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025-2026, Lawrence Livermore National Security,
 * University of Maryland Baltimore County, and the SUNDIALS contributors.
 * Copyright (c) 2013-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * Copyright (c) 2002-2013, Lawrence Livermore National Security.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * Private (package-internal) FIRKODE table routines.
 * -----------------------------------------------------------------*/

#ifndef _FIRKODE_TABLES_IMPL_H
#define _FIRKODE_TABLES_IMPL_H

#include <firkode/firkode_tables.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Construct the Radau IIA table with s stages by computing the nodes
   (roots of P_s(2x-1) - P_{s-1}(2x-1)) and the collocation
   coefficients numerically, regardless of whether a closed-form
   table exists for s. Used to cross-check the closed-form tables. */
FIRKodeTable firkTable_RadauIIA_Computed(int s);

/* Fill the derived quantities (Ainv, d, e, P) from c, A, b, and
   gamma0. Returns 0 on success, -1 on failure. */
int firkTable_Derive(FIRKodeTable T);

#ifdef __cplusplus
}
#endif

#endif
