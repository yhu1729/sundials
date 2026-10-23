/* -----------------------------------------------------------------------------
 * Programmer(s): Yifan Hu @ UMBC
 * -----------------------------------------------------------------------------
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
 * -----------------------------------------------------------------------------
 * ACSL predicates for N_Vector.
 * ---------------------------------------------------------------------------*/

#ifndef _NVECTOR_ACSL_H
#define _NVECTOR_ACSL_H

#include <nvector/nvector_serial.h>

/*@
  predicate vector_layout_serial{L}(N_Vector x) =
    \valid_read(x) &&
    \valid_read(NV_CONTENT_S(x)) &&
    NV_LENGTH_S(x) >= 0 &&
    \separated(
      x,
      NV_CONTENT_S(x),
      NV_DATA_S(x) + (0 .. NV_LENGTH_S(x) - 1));

  predicate valid_read_vector_serial{L}(N_Vector x) =
    vector_layout_serial{L}(x) &&
    \valid_read(NV_DATA_S(x) + (0 .. NV_LENGTH_S(x) - 1));

  predicate valid_vector_serial{L}(N_Vector x) =
    vector_layout_serial{L}(x) &&
    \valid(NV_DATA_S(x) + (0 .. NV_LENGTH_S(x) - 1));

  predicate separated_vector_serial{L}(N_Vector x, N_Vector y) =
    \separated(
      x,
      NV_CONTENT_S(x),
      NV_DATA_S(x) + (0 .. NV_LENGTH_S(x) - 1),
      y,
      NV_CONTENT_S(y),
      NV_DATA_S(y) + (0 .. NV_LENGTH_S(y) - 1));
*/

#endif
