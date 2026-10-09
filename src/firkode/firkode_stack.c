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
 * The stage stack (design decision D4, Section 8.3): the s stage
 * vectors of length N in one vector of length sN, for the stacked
 * Krylov solver and the stacked norm.
 *
 * The stack is an NVECTOR_MANYVECTOR over the stage vectors or, when
 * SUNDIALS is built with MPI and y0 has a communicator, an
 * NVECTOR_MPIMANYVECTOR over that communicator. Both compute the dot
 * product and the WRMS norm correctly, but the plain kind reduces
 * once per stage while the MPI kind reduces once; the plain kind's
 * max norm, min and L1 norm would be rank-local over distributed
 * stages (FIRKODE does not use them on a stack).
 *
 * The kind is decided once, in FIRKodeInit, and every access to a
 * stage dispatches on it rather than on the vector ID, which is that
 * of a ManyVector for either kind and also for a ManyVector y0. The
 * one-stage method has no stack.
 * -----------------------------------------------------------------*/

#include <stdlib.h>

#include <nvector/nvector_manyvector.h>
#include <sundials/sundials_config.h>
#if SUNDIALS_MPI_ENABLED
#include <nvector/nvector_mpimanyvector.h>
#endif

#include "firkode_impl.h"

/*---------------------------------------------------------------
  firkStackSetKind

  Decides the kind of the stages' stack from y0. The MPI kind needs
  the local reductions of the stages to reduce once over all of them;
  a y0 with a communicator but without nvdotprodlocal or
  nvwsqrsumlocal is rejected with FIRK_ILL_INPUT.
  ---------------------------------------------------------------*/
int firkStackSetKind(FIRKodeMem firk_mem, N_Vector y0)
{
  firk_mem->stack_kind = FIRK_STACK_PLAIN;

#if SUNDIALS_MPI_ENABLED
  if (N_VGetCommunicator(y0) != SUN_COMM_NULL)
  {
    if (y0->ops->nvdotprodlocal == NULL || y0->ops->nvwsqrsumlocal == NULL)
    {
      firkProcessError(firk_mem, FIRK_ILL_INPUT, __LINE__, __func__, __FILE__,
                       MSGFIRK_STACK_LOCAL_OPS);
      return FIRK_ILL_INPUT;
    }
    firk_mem->stack_kind = FIRK_STACK_MPI;
  }
#else
  (void)y0;
#endif

  return FIRK_SUCCESS;
}

/*---------------------------------------------------------------
  firkStackNew

  Creates a stack of the decided kind over the s vectors vecs, which
  it does not own; returns NULL on failure. The MPI kind duplicates
  the communicator of vecs[0], a collective operation.
  ---------------------------------------------------------------*/
N_Vector firkStackNew(FIRKodeMem firk_mem, int s, N_Vector* vecs)
{
#if SUNDIALS_MPI_ENABLED
  if (firk_mem->stack_kind == FIRK_STACK_MPI)
  {
    return N_VMake_MPIManyVector(N_VGetCommunicator(vecs[0]), (sunindextype)s,
                                 vecs, firk_mem->sunctx);
  }
#endif
  return N_VNew_ManyVector((sunindextype)s, vecs, firk_mem->sunctx);
}

/*---------------------------------------------------------------
  firkStackWeightView

  Creates the stacked weight vector: a stack whose s stages are all
  the vector w. It owns no data, is never cloned, and its
  destruction leaves w intact (H16). Returns NULL on failure.
  ---------------------------------------------------------------*/
N_Vector firkStackWeightView(FIRKodeMem firk_mem, int s, N_Vector w)
{
  N_Vector* vecs;
  N_Vector view;
  int i;

  vecs = (N_Vector*)malloc((size_t)s * sizeof(N_Vector));
  if (vecs == NULL) { return NULL; }
  for (i = 0; i < s; i++) { vecs[i] = w; }
  view = firkStackNew(firk_mem, s, vecs);
  free(vecs);
  return view;
}

/*---------------------------------------------------------------
  firkStage

  Returns the i-th stage of a stack; with one stage there is no
  stack and the "stack" is the stage itself.
  ---------------------------------------------------------------*/
N_Vector firkStage(FIRKodeMem firk_mem, N_Vector stk, int i)
{
  if (firk_mem->s_alloc < 2) { return stk; }
#if SUNDIALS_MPI_ENABLED
  if (firk_mem->stack_kind == FIRK_STACK_MPI)
  {
    return N_VGetSubvector_MPIManyVector(stk, (sunindextype)i);
  }
#else
  (void)firk_mem;
#endif
  return N_VGetSubvector_ManyVector(stk, (sunindextype)i);
}

/*---------------------------------------------------------------
  firkCloneVectorArray

  Clones w count times. Returns NULL, with nothing left allocated, if
  any clone fails, also one that returns NULL without reporting an
  error to the context, which N_VCloneVectorArray would miss (H14).
  The result is released by N_VDestroyVectorArray.
  ---------------------------------------------------------------*/
N_Vector* firkCloneVectorArray(int count, N_Vector w)
{
  N_Vector* vs;
  int i, j;

  vs = (N_Vector*)malloc((size_t)count * sizeof(N_Vector));
  if (vs == NULL) { return NULL; }
  for (j = 0; j < count; j++)
  {
    vs[j] = N_VClone(w);
    if (vs[j] == NULL)
    {
      for (i = 0; i < j; i++) { N_VDestroy(vs[i]); }
      free(vs);
      return NULL;
    }
  }
  return vs;
}
