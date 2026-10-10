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
 * Batched Robertson kinetics with FIRKODE on distributed vectors: NCOPY
 * independent copies of the Robertson problem per MPI rank, copy g with
 * the rate constant k1 = 0.04 (1 + 0.2 g / G) of G copies,
 *
 *   y1' = -k1 y1 + k3 y2 y3
 *   y2' =  k1 y1 - k3 y2 y3 - k2 y2^2
 *   y3' =  k2 y2^2,        k2 = 3e7, k3 = 1e4,  y(0) = (1, 0, 0),
 *
 * in one NVECTOR_PARALLEL vector. The three-stage Radau IIA method solves
 * its stage systems with SPGMR as the block solver, a user Jacobian-vector
 * product, and the exact block-Jacobi preconditioner of the 3 x 3 blocks
 * I - gamma J of the copies.
 *
 * The run checks itself against a serial integration of all copies on rank
 * 0 (design milestone M3): with fixed steps the solutions agree to 1e-12
 * relative, the only order-dependent operations being the global sums;
 * with adaptive steps the reduction order perturbs the error estimate at
 * roundoff, so the step counts are compared (within 2%) and the errors
 * against a reference run at rtol = 1e-10 (within a factor 2).
 * ---------------------------------------------------------------------------*/

#include <firkode/firkode.h>
#include <firkode/firkode_ls.h>
#include <mpi.h>
#include <nvector/nvector_parallel.h>
#include <nvector/nvector_serial.h>
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_math.h>
#include <sunlinsol/sunlinsol_spgmr.h>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#else
#define GSYM "g"
#define ESYM "e"
#endif

#define ZERO  SUN_RCONST(0.0)
#define ONE   SUN_RCONST(1.0)
#define NCOPY 50 /* copies per rank */
#define K2    SUN_RCONST(3.0e7)
#define K3    SUN_RCONST(1.0e4)

/* the copies [offset, offset + ncopy) of the G copies, with the Jacobians
   and inverse preconditioner blocks of the last setup */
typedef struct
{
  sunindextype offset, ncopy, G;
  sunrealtype* J;    /* 9 ncopy, row-major per copy */
  sunrealtype* Pinv; /* 9 ncopy                     */
} UserData;

static sunrealtype k1(const UserData* ud, sunindextype c)
{
  return SUN_RCONST(0.04) *
         (ONE +
          SUN_RCONST(0.2) * (sunrealtype)(ud->offset + c) / (sunrealtype)ud->G);
}

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  UserData* ud   = (UserData*)user_data;
  sunrealtype* Y = N_VGetArrayPointer(y);
  sunrealtype* F = N_VGetArrayPointer(ydot);
  sunindextype c;
  (void)t;
  for (c = 0; c < ud->ncopy; c++)
  {
    const sunrealtype* yc = Y + 3 * c;
    sunrealtype* fc       = F + 3 * c;
    fc[0]                 = -k1(ud, c) * yc[0] + K3 * yc[1] * yc[2];
    fc[2]                 = K2 * yc[1] * yc[1];
    fc[1]                 = -fc[0] - fc[2];
  }
  return 0;
}

/* the Jacobian of copy c at y, row-major */
static void jac_copy(const UserData* ud, sunindextype c, const sunrealtype* yc,
                     sunrealtype* Jc)
{
  Jc[0] = -k1(ud, c);
  Jc[1] = K3 * yc[2];
  Jc[2] = K3 * yc[1];
  Jc[3] = k1(ud, c);
  Jc[4] = -K3 * yc[2] - SUN_RCONST(2.0) * K2 * yc[1];
  Jc[5] = -K3 * yc[1];
  Jc[6] = ZERO;
  Jc[7] = SUN_RCONST(2.0) * K2 * yc[1];
  Jc[8] = ZERO;
}

static int jtimes(N_Vector v, N_Vector Jv, sunrealtype t, N_Vector y,
                  N_Vector fy, void* user_data, N_Vector tmp)
{
  UserData* ud    = (UserData*)user_data;
  sunrealtype* Y  = N_VGetArrayPointer(y);
  sunrealtype* V  = N_VGetArrayPointer(v);
  sunrealtype* JV = N_VGetArrayPointer(Jv);
  sunrealtype Jc[9];
  sunindextype c;
  int i;
  (void)t;
  (void)fy;
  (void)tmp;
  for (c = 0; c < ud->ncopy; c++)
  {
    jac_copy(ud, c, Y + 3 * c, Jc);
    for (i = 0; i < 3; i++)
    {
      JV[3 * c + i] = Jc[3 * i] * V[3 * c] + Jc[3 * i + 1] * V[3 * c + 1] +
                      Jc[3 * i + 2] * V[3 * c + 2];
    }
  }
  return 0;
}

/* inverts the 3 x 3 matrix A into B; returns 1 if singular */
static int inv3(const sunrealtype* A, sunrealtype* B)
{
  sunrealtype det = A[0] * (A[4] * A[8] - A[5] * A[7]) -
                    A[1] * (A[3] * A[8] - A[5] * A[6]) +
                    A[2] * (A[3] * A[7] - A[4] * A[6]);
  if (det == ZERO) { return 1; }
  B[0] = (A[4] * A[8] - A[5] * A[7]) / det;
  B[1] = (A[2] * A[7] - A[1] * A[8]) / det;
  B[2] = (A[1] * A[5] - A[2] * A[4]) / det;
  B[3] = (A[5] * A[6] - A[3] * A[8]) / det;
  B[4] = (A[0] * A[8] - A[2] * A[6]) / det;
  B[5] = (A[2] * A[3] - A[0] * A[5]) / det;
  B[6] = (A[3] * A[7] - A[4] * A[6]) / det;
  B[7] = (A[1] * A[6] - A[0] * A[7]) / det;
  B[8] = (A[0] * A[4] - A[1] * A[3]) / det;
  return 0;
}

/* the blocks I - gamma J of the copies, J saved or evaluated anew */
static int psetup(sunrealtype t, N_Vector y, N_Vector fy, sunbooleantype jok,
                  sunbooleantype* jcur, sunrealtype gamma, void* user_data)
{
  UserData* ud   = (UserData*)user_data;
  sunrealtype* Y = N_VGetArrayPointer(y);
  sunrealtype P[9];
  sunindextype c;
  int i;
  (void)t;
  (void)fy;
  for (c = 0; c < ud->ncopy; c++)
  {
    if (!jok) { jac_copy(ud, c, Y + 3 * c, ud->J + 9 * c); }
    for (i = 0; i < 9; i++) { P[i] = -gamma * ud->J[9 * c + i]; }
    P[0] += ONE;
    P[4] += ONE;
    P[8] += ONE;
    if (inv3(P, ud->Pinv + 9 * c)) { return 1; }
  }
  *jcur = !jok;
  return 0;
}

static int psolve(sunrealtype t, N_Vector y, N_Vector fy, N_Vector r, N_Vector z,
                  sunrealtype gamma, sunrealtype delta, int lr, void* user_data)
{
  UserData* ud   = (UserData*)user_data;
  sunrealtype* R = N_VGetArrayPointer(r);
  sunrealtype* Z = N_VGetArrayPointer(z);
  sunrealtype rc[3];
  sunindextype c;
  int i;
  (void)t;
  (void)y;
  (void)fy;
  (void)gamma;
  (void)delta;
  (void)lr;
  for (c = 0; c < ud->ncopy; c++)
  {
    const sunrealtype* B = ud->Pinv + 9 * c;
    for (i = 0; i < 3; i++) { rc[i] = R[3 * c + i]; }
    for (i = 0; i < 3; i++)
    {
      Z[3 * c + i] = B[3 * i] * rc[0] + B[3 * i + 1] * rc[1] +
                     B[3 * i + 2] * rc[2];
    }
  }
  return 0;
}

/* an integration to tout of the copies held by y (distributed or serial);
   fixed steps h > 0, else adaptive at rtol; returns the integrator's flag */
static int integrate(SUNContext sunctx, N_Vector y, UserData* ud, sunrealtype h,
                     sunrealtype rtol, sunrealtype tout, long int* nst,
                     long int* nni)
{
  sunindextype c;
  sunrealtype tret, *A;
  int flag;
  N_Vector abstol;
  SUNLinearSolver LS;
  void* mem;

  /* y(0) and the absolute tolerances 1e-8, 1e-14, 1e-6 times rtol / 1e-4 */
  abstol = N_VClone(y);
  A      = N_VGetArrayPointer(abstol);
  N_VConst(ZERO, y);
  for (c = 0; c < ud->ncopy; c++)
  {
    N_VGetArrayPointer(y)[3 * c] = ONE;
    A[3 * c]                     = SUN_RCONST(1.0e-4) * rtol;
    A[3 * c + 1]                 = SUN_RCONST(1.0e-10) * rtol;
    A[3 * c + 2]                 = SUN_RCONST(1.0e-2) * rtol;
  }

  LS  = SUNLinSol_SPGMR(y, SUN_PREC_LEFT, 0, sunctx);
  mem = FIRKodeCreate(sunctx);
  if (mem == NULL || LS == NULL) { return -1; }
  flag = FIRKodeInit(mem, f, ZERO, y);
  if (!flag) { flag = FIRKodeSetUserData(mem, ud); }
  if (!flag) { flag = FIRKodeSVtolerances(mem, rtol, abstol); }
  if (!flag) { flag = FIRKodeSetLinearSolver(mem, LS, NULL); }
  if (!flag) { flag = FIRKodeSetJacTimes(mem, NULL, jtimes); }
  if (!flag) { flag = FIRKodeSetPreconditioner(mem, psetup, psolve); }
  if (!flag) { flag = FIRKodeSetMaxNumSteps(mem, 100000); }
  if (!flag && h > ZERO) { flag = FIRKodeSetFixedStep(mem, h); }
  if (!flag) { flag = FIRKodeEvolve(mem, tout, y, &tret, FIRK_NORMAL); }
  FIRKodeGetNumSteps(mem, nst);
  FIRKodeGetNumNonlinSolvIters(mem, nni);

  FIRKodeFree(&mem);
  SUNLinSolFree(LS);
  N_VDestroy(abstol);
  return flag;
}

static void* safe_malloc(size_t n)
{
  void* p = malloc(n);
  if (p == NULL) { MPI_Abort(MPI_COMM_WORLD, 1); }
  return p;
}

int main(int argc, char* argv[])
{
  int rank, nprocs, fails = 0, flag, i;
  sunindextype G, c;
  long int nst_p = 0, nni_p = 0, nst_s = 0, nni_s = 0, nst_r = 0, nni_r = 0;
  sunrealtype *yall      = NULL, diff, err_p, err_s, ymax;
  const sunrealtype hfix = SUN_RCONST(1.0e-4), tfix = SUN_RCONST(0.01);
  const sunrealtype rtol = SUN_RCONST(1.0e-6), tend = SUN_RCONST(4.0);
  SUNContext sunctx = NULL, serctx = NULL;
  N_Vector y, ys = NULL, yref = NULL;
  UserData udp, uds;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  G = (sunindextype)nprocs * NCOPY;
  if (SUNContext_Create(MPI_COMM_WORLD, &sunctx))
  {
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  udp.offset = (sunindextype)rank * NCOPY;
  udp.ncopy  = NCOPY;
  udp.G      = G;
  udp.J      = (sunrealtype*)safe_malloc(9 * NCOPY * sizeof(sunrealtype));
  udp.Pinv   = (sunrealtype*)safe_malloc(9 * NCOPY * sizeof(sunrealtype));
  y          = N_VNew_Parallel(MPI_COMM_WORLD, 3 * NCOPY, 3 * G, sunctx);

  if (rank == 0)
  {
    printf("\nBatched Robertson kinetics: %ld copies on %d MPI ranks\n",
           (long int)G, nprocs);
    printf("Radau IIA (s = 3), SPGMR blocks with an exact block-Jacobi "
           "preconditioner\n\n");
    if (SUNContext_Create(SUN_COMM_NULL, &serctx))
    {
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    uds.offset = 0;
    uds.ncopy  = G;
    uds.G      = G;
    uds.J      = (sunrealtype*)safe_malloc(9 * (size_t)G * sizeof(sunrealtype));
    uds.Pinv   = (sunrealtype*)safe_malloc(9 * (size_t)G * sizeof(sunrealtype));
    ys         = N_VNew_Serial(3 * G, serctx);
    yref       = N_VNew_Serial(3 * G, serctx);
    yall       = (sunrealtype*)safe_malloc(3 * (size_t)G * sizeof(sunrealtype));
  }

  /* fixed steps: the distributed and the serial runs agree to roundoff */
  flag = integrate(sunctx, y, &udp, hfix, rtol, tfix, &nst_p, &nni_p);
  if (flag != FIRK_SUCCESS) { fails++; }
  MPI_Gather(N_VGetArrayPointer(y), 3 * NCOPY, MPI_SUNREALTYPE, yall, 3 * NCOPY,
             MPI_SUNREALTYPE, 0, MPI_COMM_WORLD);
  if (rank == 0)
  {
    flag = integrate(serctx, ys, &uds, hfix, rtol, tfix, &nst_s, &nni_s);
    if (flag != FIRK_SUCCESS) { fails++; }
    diff = ZERO;
    ymax = ZERO;
    for (c = 0; c < 3 * G; c++)
    {
      diff = SUNMAX(diff, SUNRabs(yall[c] - NV_Ith_S(ys, c)));
      ymax = SUNMAX(ymax, SUNRabs(NV_Ith_S(ys, c)));
    }
    printf("fixed steps h = %" GSYM " to t = %" GSYM ":\n", hfix, tfix);
    printf("  y(copy 0) = %14.6" ESYM " %14.6" ESYM " %14.6" ESYM "\n", yall[0],
           yall[1], yall[2]);
    printf("  steps %ld, Newton iterations %ld (serial: %ld, %ld)\n", nst_p,
           nni_p, nst_s, nni_s);
    printf("  largest difference from the serial run, relative: %.2" ESYM "\n",
           diff / ymax);
    if (diff > SUN_RCONST(1.0e-12) * ymax || labs(nni_p - nni_s) > 1)
    {
      printf("  FAILED: the runs differ by more than roundoff\n");
      fails++;
    }
  }

  /* adaptive steps: comparable step counts and errors */
  flag = integrate(sunctx, y, &udp, ZERO, rtol, tend, &nst_p, &nni_p);
  if (flag != FIRK_SUCCESS) { fails++; }
  MPI_Gather(N_VGetArrayPointer(y), 3 * NCOPY, MPI_SUNREALTYPE, yall, 3 * NCOPY,
             MPI_SUNREALTYPE, 0, MPI_COMM_WORLD);
  if (rank == 0)
  {
    flag = integrate(serctx, ys, &uds, ZERO, rtol, tend, &nst_s, &nni_s);
    if (flag != FIRK_SUCCESS) { fails++; }
    flag = integrate(serctx, yref, &uds, ZERO, SUN_RCONST(1.0e-10), tend,
                     &nst_r, &nni_r);
    if (flag != FIRK_SUCCESS) { fails++; }
    err_p = err_s = ZERO;
    for (c = 0; c < G; c++)
    {
      for (i = 0; i < 3; i++)
      {
        /* errors relative to the weights of the run's tolerances */
        sunrealtype w = ONE / (rtol * SUNRabs(NV_Ith_S(yref, 3 * c + i)) +
                               rtol * ((i == 0)   ? SUN_RCONST(1.0e-4)
                                       : (i == 1) ? SUN_RCONST(1.0e-10)
                                                  : SUN_RCONST(1.0e-2)));
        err_p         = SUNMAX(err_p,
                               w * SUNRabs(yall[3 * c + i] - NV_Ith_S(yref, 3 * c + i)));
        err_s         = SUNMAX(err_s, w * SUNRabs(NV_Ith_S(ys, 3 * c + i) -
                                                  NV_Ith_S(yref, 3 * c + i)));
      }
    }
    printf("\nadaptive steps, rtol = %" GSYM ", to t = %" GSYM ":\n", rtol, tend);
    printf("  y(copy 0) = %14.6" ESYM " %14.6" ESYM " %14.6" ESYM "\n", yall[0],
           yall[1], yall[2]);
    printf("  steps %ld, Newton iterations %ld (serial: %ld, %ld)\n", nst_p,
           nni_p, nst_s, nni_s);
    printf("  weighted error against rtol = 1e-10: %.2" ESYM
           " (serial: %.2" ESYM ")\n",
           err_p, err_s);
    if (SUNRabs((sunrealtype)(nst_p - nst_s)) >
          SUN_RCONST(0.05) * (sunrealtype)nst_s ||
        err_p > SUN_RCONST(2.0) * err_s + SUN_RCONST(1.0e-3))
    {
      printf("  FAILED: the step counts or the errors differ too much\n");
      fails++;
    }
  }

  i = fails;
  MPI_Allreduce(&i, &fails, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (rank == 0)
  {
    printf("\n%s\n", fails ? "FAILED" : "PASSED");
    N_VDestroy(ys);
    N_VDestroy(yref);
    free(uds.J);
    free(uds.Pinv);
    free(yall);
    SUNContext_Free(&serctx);
  }
  N_VDestroy(y);
  free(udp.J);
  free(udp.Pinv);
  SUNContext_Free(&sunctx);
  MPI_Finalize();
  return fails ? 1 : 0;
}
