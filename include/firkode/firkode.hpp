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
 * C++ specific FIRKODE definitions.
 * ---------------------------------------------------------------------------*/

#ifndef _FIRKODE_HPP
#define _FIRKODE_HPP

#include <sundials/sundials_classview.hpp>

#include <firkode/firkode.h>

namespace sundials {
namespace experimental {

struct FIRKodeDeleter
{
  void operator()(void* v) { FIRKodeFree(&v); }
};

using FIRKodeView = ClassView<void*, FIRKodeDeleter>;

} // namespace experimental
} // namespace sundials

#endif
