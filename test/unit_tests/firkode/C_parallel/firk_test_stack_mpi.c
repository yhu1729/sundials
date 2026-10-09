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
 * Unit test of the stage stack with a distributed vector (design Section 8.3,
 * decision D4) on y' = diag(lambda) y, NLOC components per rank, lambda_g =
 * -10^(g mod 4) (1 + g / N) for the global index g, y(0) = 1, the three-stage
 * method, and SPGMR with a user jtimes and an exact block preconditioner:
 *
 * - a y0 with a communicator selects the MPIManyVector kind (white box);
 * - FIRKodeInit rejects a y0 with a communicator but without nvdotprodlocal
 *   or nvwsqrsumlocal (H9);
 * - no forbidden operation of Table 8.8 on the MPIManyVector stacks, the
 *   multi-vector operations only inside the stacked FGMRES (H9, H15), and
 *   no clone during the steps (H14);
 * - fixed steps on all ranks agree with a serial run of the whole problem
 *   on rank 0 to 1e-12, with the same Newton iterations.
 *
 * Pass criterion: every check holds on every rank.
 * ---------------------------------------------------------------------------*/

#include <mpi.h>
#include <stdio.h>

#include "firkode/firkode.h"
#include "firkode/firkode_impl.h" /* the stacks */
#include "firkode/firkode_ls.h"
#include "firkode/firkode_ls_impl.h" /* the stacked FGMRES */
#include "nvector/nvector_mpimanyvector.h"
#include "nvector/nvector_parallel.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_math.h"
#include "sunlinsol/sunlinsol_spgmr.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)
#define NLOC 4

static int rank = 0, nprocs = 1;
static sunindextype NGLOB = NLOC;

/* the components [offset, offset + length) of the global problem */
struct block
{
  sunindextype offset;
};

static sunrealtype lambda(sunindextype g)
{
  return -SUNRpowerI(SUN_RCONST(10.0), (int)(g % 4)) *
         (ONE + (sunrealtype)g / (sunrealtype)NGLOB);
}

static sunindextype local_length(N_Vector v)
{
  return (N_VGetVectorID(v) == SUNDIALS_NVEC_PARALLEL) ? NV_LOCLENGTH_P(v)
                                                       : N_VGetLength(v);
}

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  struct block* b = (struct block*)user_data;
  sunrealtype* yd = N_VGetArrayPointer(y);
  sunrealtype* fd = N_VGetArrayPointer(ydot);
  sunindextype i, n = local_length(y);
  (void)t;
  for (i = 0; i < n; i++) { fd[i] = lambda(b->offset + i) * yd[i]; }
  return 0;
}

static int jtimes(N_Vector v, N_Vector Jv, sunrealtype t, N_Vector y,
                  N_Vector fy, void* user_data, N_Vector tmp)
{
  struct block* b = (struct block*)user_data;
  sunrealtype* vd = N_VGetArrayPointer(v);
  sunrealtype* jd = N_VGetArrayPointer(Jv);
  sunindextype i, n = local_length(v);
  (void)t;
  (void)y;
  (void)fy;
  (void)tmp;
  for (i = 0; i < n; i++) { jd[i] = lambda(b->offset + i) * vd[i]; }
  return 0;
}

static int psolve(sunrealtype t, N_Vector y, N_Vector fy, N_Vector r, N_Vector z,
                  sunrealtype gamma, sunrealtype delta, int lr, void* user_data)
{
  struct block* b = (struct block*)user_data;
  sunrealtype* rd = N_VGetArrayPointer(r);
  sunrealtype* zd = N_VGetArrayPointer(z);
  sunindextype i, n = local_length(r);
  (void)t;
  (void)y;
  (void)fy;
  (void)delta;
  (void)lr;
  for (i = 0; i < n; i++)
  {
    zd[i] = rd[i] / (ONE - gamma * lambda(b->offset + i));
  }
  return 0;
}

/*---------------------------------------------------------------
  Clone counting (H14) and the operations of the stacks (H9, H15)
  ---------------------------------------------------------------*/

static long int nclones               = 0;
static long int ntrapped              = 0;
static long int nmulti_outside        = 0;
static int in_fgmr                    = 0;
static int (*fgmr_solve)(SUNLinearSolver, SUNMatrix, N_Vector, N_Vector,
                         sunrealtype) = NULL;

static N_Vector clone_counting(N_Vector w)
{
  nclones++;
  return N_VClone_Parallel(w); /* copies the operations of w */
}

static sunrealtype trap_norm(N_Vector x)
{
  (void)x;
  ntrapped++;
  return ZERO;
}

static sunrealtype trap_wnorm(N_Vector x, N_Vector w)
{
  (void)x;
  (void)w;
  ntrapped++;
  return ZERO;
}

static sunbooleantype trap_test2(N_Vector x, N_Vector z)
{
  (void)x;
  (void)z;
  ntrapped++;
  return SUNTRUE;
}

static sunbooleantype trap_test3(N_Vector c, N_Vector x, N_Vector m)
{
  (void)c;
  (void)x;
  (void)m;
  ntrapped++;
  return SUNTRUE;
}

static sunindextype trap_length(N_Vector x)
{
  (void)x;
  ntrapped++;
  return 0;
}

static sunrealtype* trap_data(N_Vector x)
{
  (void)x;
  ntrapped++;
  return NULL;
}

static SUNErrCode linearcombination_counting(int nvec, sunrealtype* c,
                                             N_Vector* X, N_Vector z)
{
  if (!in_fgmr) { nmulti_outside++; }
  return N_VLinearCombination_MPIManyVector(nvec, c, X, z);
}

static SUNErrCode scaleaddmulti_counting(int nvec, sunrealtype* a, N_Vector x,
                                         N_Vector* Y, N_Vector* Z)
{
  if (!in_fgmr) { nmulti_outside++; }
  return N_VScaleAddMulti_MPIManyVector(nvec, a, x, Y, Z);
}

static SUNErrCode dotprodmulti_counting(int nvec, N_Vector x, N_Vector* Y,
                                        sunrealtype* dotprods)
{
  if (!in_fgmr) { nmulti_outside++; }
  return N_VDotProdMulti_MPIManyVector(nvec, x, Y, dotprods);
}

static void replace_stack_ops(N_Vector stk)
{
  stk->ops->nvmaxnorm           = trap_norm;
  stk->ops->nvmin               = trap_norm;
  stk->ops->nvl1norm            = trap_norm;
  stk->ops->nvwl2norm           = trap_wnorm;
  stk->ops->nvminquotient       = trap_wnorm;
  stk->ops->nvinvtest           = trap_test2;
  stk->ops->nvconstrmask        = trap_test3;
  stk->ops->nvgetlocallength    = trap_length;
  stk->ops->nvgetarraypointer   = trap_data;
  stk->ops->nvlinearcombination = linearcombination_counting;
  stk->ops->nvscaleaddmulti     = scaleaddmulti_counting;
  stk->ops->nvdotprodmulti      = dotprodmulti_counting;
}

static int fgmr_solve_marked(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                             N_Vector b, sunrealtype tol)
{
  int retval;
  in_fgmr = 1;
  retval  = fgmr_solve(S, A, x, b, tol);
  in_fgmr = 0;
  return retval;
}

static int check(int ok, const char* what)
{
  if (!ok)
  {
    printf("  FAIL (rank %d): %s\n", rank, what);
    return 1;
  }
  return 0;
}

/*---------------------------------------------------------------
  The runs
  ---------------------------------------------------------------*/

/* fixed steps to tout with the vector y (serial or parallel); returns the
   number of Newton iterations in *nni */
static int run_fixed(SUNContext sunctx, N_Vector y, struct block* b,
                     long int* nni, int instrument)
{
  int fails = 0;
  sunrealtype tret;
  SUNLinearSolver LS;
  FIRKodeMem fm;
  FIRKLsMem ls;
  long int nclones_step1;
  char what[160];
  void* mem;

  N_VConst(ONE, y);
  LS  = SUNLinSol_SPGMR(y, SUN_PREC_LEFT, 0, sunctx);
  mem = FIRKodeCreate(sunctx);
  if (FIRKodeInit(mem, f, ZERO, y) || FIRKodeSetUserData(mem, b) ||
      FIRKodeSStolerances(mem, SUN_RCONST(1.0e-8), SUN_RCONST(1.0e-8)) ||
      FIRKodeSetLinearSolver(mem, LS, NULL) ||
      FIRKodeSetJacTimes(mem, NULL, jtimes) ||
      FIRKodeSetPreconditioner(mem, NULL, psolve) ||
      FIRKodeSetFixedStep(mem, SUN_RCONST(0.01)))
  {
    return 1;
  }

  fm = (FIRKodeMem)mem;
  if (instrument)
  {
    /* the MPI kind, and its stacks' operations replaced before the stacked
       FGMRES clones them */
    fails += check(fm->stack_kind == FIRK_STACK_MPI &&
                     N_VGetVectorID(fm->g_stk) == SUNDIALS_NVEC_MPIMANYVECTOR,
                   "a y0 with a communicator gives MPIManyVector stacks");
    ntrapped = nmulti_outside = 0;
    replace_stack_ops(fm->g_stk);
    replace_stack_ops(fm->delta_stk);
    replace_stack_ops(fm->ewt_stk);
    if (fm->linit(fm) != 0) { return fails + 1; }
    ls                     = (FIRKLsMem)fm->lmem;
    fgmr_solve             = ls->SPFGMR->ops->solve;
    ls->SPFGMR->ops->solve = fgmr_solve_marked;
  }

  fails += check(FIRKodeEvolve(mem, SUN_RCONST(0.1), y, &tret, FIRK_ONE_STEP) ==
                   FIRK_SUCCESS,
                 "the first step succeeds");
  nclones_step1 = nclones;
  fails += check(FIRKodeEvolve(mem, SUN_RCONST(0.1), y, &tret, FIRK_NORMAL) ==
                   FIRK_SUCCESS,
                 "the integration succeeds");
  FIRKodeGetNumNonlinSolvIters(mem, nni);

  if (instrument)
  {
    snprintf(what, sizeof(what),
             "H9: %ld forbidden operations on the MPIManyVector stacks",
             ntrapped);
    fails += check(ntrapped == 0, what);
    snprintf(what, sizeof(what),
             "H15: %ld multi-vector operations on stacks outside FGMRES",
             nmulti_outside);
    fails += check(nmulti_outside == 0, what);
    snprintf(what, sizeof(what), "H14: %ld clones during the steps",
             nclones - nclones_step1);
    fails += check(nclones == nclones_step1, what);
  }

  FIRKodeFree(&mem);
  SUNLinSolFree(LS);
  return fails;
}

int main(int argc, char* argv[])
{
  SUNContext sunctx = NULL;
  int fails = 0, allfails = 0, flag;
  long int nni_par = 0, nni_ser = 0;
  sunrealtype diff = ZERO;
  sunrealtype yglob[64];
  sunindextype i;
  N_Vector y, ybad;
  struct block bpar, bser;
  void* mem;
  char what[160];

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  NGLOB = (sunindextype)nprocs * NLOC;
  if (NGLOB > 64 || SUNContext_Create(MPI_COMM_WORLD, &sunctx))
  {
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  /* the distributed y0, counting its clones */
  y               = N_VNew_Parallel(MPI_COMM_WORLD, NLOC, NGLOB, sunctx);
  y->ops->nvclone = clone_counting;
  bpar.offset     = (sunindextype)rank * NLOC;

  /* H9: Init rejects a y0 without the local reductions */
  ybad                      = N_VClone_Parallel(y);
  ybad->ops->nvdotprodlocal = NULL;
  mem                       = FIRKodeCreate(sunctx);
  flag                      = FIRKodeInit(mem, f, ZERO, ybad);
  fails += check(flag == FIRK_ILL_INPUT,
                 "H9: Init rejects a y0 without nvdotprodlocal");
  FIRKodeFree(&mem);
  ybad->ops->nvdotprodlocal = y->ops->nvdotprodlocal;
  ybad->ops->nvwsqrsumlocal = NULL;
  mem                       = FIRKodeCreate(sunctx);
  flag                      = FIRKodeInit(mem, f, ZERO, ybad);
  fails += check(flag == FIRK_ILL_INPUT,
                 "H9: Init rejects a y0 without nvwsqrsumlocal");
  FIRKodeFree(&mem);
  N_VDestroy(ybad);

  /* the distributed run, instrumented */
  fails += run_fixed(sunctx, y, &bpar, &nni_par, 1);

  /* gather it to rank 0 and compare with a serial run of the whole
     problem there */
  MPI_Gather(N_VGetArrayPointer(y), NLOC, MPI_SUNREALTYPE, yglob, NLOC,
             MPI_SUNREALTYPE, 0, MPI_COMM_WORLD);
  if (rank == 0)
  {
    SUNContext serctx = NULL;
    N_Vector yser;
    if (SUNContext_Create(SUN_COMM_NULL, &serctx))
    {
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    yser        = N_VNew_Serial(NGLOB, serctx);
    bser.offset = 0;
    fails += run_fixed(serctx, yser, &bser, &nni_ser, 0);
    for (i = 0; i < NGLOB; i++)
    {
      diff = SUNMAX(diff, SUNRabs(yglob[i] - NV_Ith_S(yser, i)));
    }
    snprintf(what, sizeof(what),
             "%d ranks agree with the serial run: difference %g, Newton "
             "iterations %ld and %ld",
             nprocs, (double)diff, nni_par, nni_ser);
    fails += check(diff <= SUN_RCONST(1.0e-12) && nni_par == nni_ser, what);
    N_VDestroy(yser);
    SUNContext_Free(&serctx);
  }

  N_VDestroy(y);
  SUNContext_Free(&sunctx);

  MPI_Allreduce(&fails, &allfails, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (rank == 0)
  {
    if (allfails) { printf("\n%d test(s) FAILED\n", allfails); }
    else { printf("All tests passed\n"); }
  }
  MPI_Finalize();
  return allfails ? 1 : 0;
}
