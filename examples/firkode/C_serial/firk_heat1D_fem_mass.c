/*---------------------------------------------------------------
 * Programmer(s): Yifan Hu @ UMBC
 *---------------------------------------------------------------
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
 *---------------------------------------------------------------
 * Example problem:
 *
 * The following test simulates the 1D heat equation
 *    u_t = k*u_xx + f(t,x)
 * for t in [0, 1], x in [0, 1], with homogeneous Dirichlet
 * boundary conditions and the source term
 *    f(t,x) = (k*pi^2 - 1) exp(-t) sin(pi x)
 * chosen so that the exact solution is u(t,x) = exp(-t) sin(pi x).
 *
 * The equation is discretized in space with linear (P1) finite
 * elements on a uniform mesh of N nodes. The Galerkin
 * semi-discretization has the form
 *    M u' = -k K u + M f(t),
 * where M and K are the tridiagonal mass and stiffness matrices
 * and the source term is interpolated at the mesh nodes. The
 * boundary rows of M are replaced by the identity and those of K
 * and f by zero so that the boundary values remain zero.
 *
 * This program integrates the semi-discrete system with the
 * 3-stage Radau IIA method of FIRKODE, using band SUNLinearSolver
 * objects for both the Newton matrices M - gamma*J and the mass
 * matrix M, with user-supplied Jacobian and mass matrix routines.
 *
 * The maximum nodal error with respect to the exact solution is
 * printed at 10 output times, and run statistics are printed at
 * the end.
 *---------------------------------------------------------------*/

/* Header files */
#include <firkode/firkode.h>    /* prototypes for FIRKODE fcts., consts */
#include <firkode/firkode_ls.h> /* FIRKODE linear solver interface      */
#include <math.h>
#include <nvector/nvector_serial.h> /* serial N_Vector types, fcts., macros */
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_types.h> /* defs. of sunrealtype, sunindextype, etc */
#include <sunlinsol/sunlinsol_band.h> /* access to band SUNLinearSolver       */
#include <sunmatrix/sunmatrix_band.h> /* access to band SUNMatrix             */

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)
#define TWO  SUN_RCONST(2.0)
#define FOUR SUN_RCONST(4.0)
#define SIX  SUN_RCONST(6.0)
#define PI   SUN_RCONST(3.14159265358979323846)

/* user data structure */
typedef struct
{
  sunindextype N; /* number of nodes       */
  sunrealtype dx; /* mesh spacing          */
  sunrealtype k;  /* diffusion coefficient */
}* UserData;

/* User-supplied Functions Called by the Solver */
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int MassMatrix(sunrealtype t, SUNMatrix M, void* user_data,
                      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private helper functions */
static sunrealtype exact(sunrealtype t, sunrealtype x);
static sunrealtype max_error(sunrealtype t, N_Vector y, UserData udata);
static int check_flag(void* flagvalue, const char* funcname, int opt);

/* Main Program */
int main(void)
{
  /* general problem parameters */
  sunrealtype T0   = SUN_RCONST(0.0);     /* initial time */
  sunrealtype Tf   = SUN_RCONST(1.0);     /* final time */
  int Nt           = 10;                  /* total number of output times */
  sunrealtype rtol = SUN_RCONST(1.0e-6);  /* relative tolerance */
  sunrealtype atol = SUN_RCONST(1.0e-10); /* absolute tolerance */
  UserData udata   = NULL;
  sunrealtype* data;
  sunindextype N = 101;             /* spatial mesh size */
  sunrealtype k  = SUN_RCONST(1.0); /* heat conductivity */
  sunindextype i;

  /* general problem variables */
  int flag;                   /* reusable error-checking flag */
  N_Vector y          = NULL; /* empty vector for storing solution */
  SUNMatrix A         = NULL; /* empty matrix for the Newton systems */
  SUNMatrix M         = NULL; /* empty mass matrix */
  SUNLinearSolver LS  = NULL; /* empty linear solver object */
  SUNLinearSolver MLS = NULL; /* empty mass matrix solver object */
  void* firkode_mem   = NULL; /* empty FIRKODE memory structure */
  sunrealtype t, dTout, tout;
  int iout;
  long int nst, nst_a, nfe, nsetups, nje, nni, nnf, ncfn, netf;
  long int nmset, nms, nMv;

  /* Create the SUNDIALS context object for this simulation */
  SUNContext ctx;
  flag = SUNContext_Create(SUN_COMM_NULL, &ctx);
  if (check_flag(&flag, "SUNContext_Create", 1)) { return 1; }

  /* allocate and fill udata structure */
  udata = (UserData)malloc(sizeof(*udata));
  if (check_flag((void*)udata, "malloc", 2)) { return 1; }
  udata->N  = N;
  udata->k  = k;
  udata->dx = ONE / (sunrealtype)(N - 1); /* mesh spacing */

  /* Initial problem output */
  printf("\n1D Heat PDE test problem (P1 finite elements with mass matrix):\n");
  printf("  N = %li\n", (long int)udata->N);
  printf("  diffusion coefficient:  k = %" GSYM "\n", udata->k);

  /* Initialize data structures */
  y = N_VNew_Serial(N, ctx); /* Create serial vector for solution */
  if (check_flag((void*)y, "N_VNew_Serial", 0)) { return 1; }
  data = N_VGetArrayPointer(y);
  for (i = 0; i < N; i++) { data[i] = exact(T0, udata->dx * (sunrealtype)i); }

  /* Call FIRKodeCreate and FIRKodeInit to initialize the integrator and
     specify the right-hand side function in M y'=f(t,y), the initial time
     T0, and the initial dependent variable vector y. */
  firkode_mem = FIRKodeCreate(ctx);
  if (check_flag((void*)firkode_mem, "FIRKodeCreate", 0)) { return 1; }
  flag = FIRKodeInit(firkode_mem, f, T0, y);
  if (check_flag(&flag, "FIRKodeInit", 1)) { return 1; }

  /* Set routines */
  flag = FIRKodeSetUserData(firkode_mem, (void*)udata); /* Pass udata */
  if (check_flag(&flag, "FIRKodeSetUserData", 1)) { return 1; }
  flag = FIRKodeSetMaxNumSteps(firkode_mem, 10000); /* Increase max num steps */
  if (check_flag(&flag, "FIRKodeSetMaxNumSteps", 1)) { return 1; }
  flag = FIRKodeSStolerances(firkode_mem, rtol, atol); /* Specify tolerances */
  if (check_flag(&flag, "FIRKodeSStolerances", 1)) { return 1; }

  /* Initialize band matrices and linear solvers for the Newton systems and
     the mass matrix (both tridiagonal) */
  A = SUNBandMatrix(N, 1, 1, ctx);
  if (check_flag((void*)A, "SUNBandMatrix", 0)) { return 1; }
  LS = SUNLinSol_Band(y, A, ctx);
  if (check_flag((void*)LS, "SUNLinSol_Band", 0)) { return 1; }
  M = SUNBandMatrix(N, 1, 1, ctx);
  if (check_flag((void*)M, "SUNBandMatrix", 0)) { return 1; }
  MLS = SUNLinSol_Band(y, M, ctx);
  if (check_flag((void*)MLS, "SUNLinSol_Band", 0)) { return 1; }

  /* Linear solver interface */
  flag = FIRKodeSetLinearSolver(firkode_mem, LS, A); /* Attach matrix and LS */
  if (check_flag(&flag, "FIRKodeSetLinearSolver", 1)) { return 1; }
  flag = FIRKodeSetJacFn(firkode_mem, Jac); /* Set the Jacobian routine */
  if (check_flag(&flag, "FIRKodeSetJacFn", 1)) { return 1; }

  /* Mass matrix solver interface: the mass matrix is time-independent */
  flag = FIRKodeSetMassLinearSolver(firkode_mem, MLS, M, SUNFALSE);
  if (check_flag(&flag, "FIRKodeSetMassLinearSolver", 1)) { return 1; }
  flag = FIRKodeSetMassFn(firkode_mem, MassMatrix); /* Set the mass routine */
  if (check_flag(&flag, "FIRKodeSetMassFn", 1)) { return 1; }

  /* Main time-stepping loop: calls FIRKodeEvolve to perform the integration,
     then prints results.  Stops when the final time has been reached */
  t     = T0;
  dTout = (Tf - T0) / Nt;
  tout  = T0 + dTout;
  printf("        t        max error\n");
  printf("   -------------------------\n");
  printf("  %10.6" FSYM "  %10.3" ESYM "\n", t, max_error(t, y, udata));
  for (iout = 0; iout < Nt; iout++)
  {
    flag = FIRKodeEvolve(firkode_mem, tout, y, &t,
                         FIRK_NORMAL); /* call integrator */
    if (check_flag(&flag, "FIRKodeEvolve", 1)) { break; }
    printf("  %10.6" FSYM "  %10.3" ESYM "\n", t,
           max_error(t, y, udata)); /* print solution error */
    if (flag >= 0)
    { /* successful solve: update output time */
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    }
    else
    { /* unsuccessful solve: break */
      fprintf(stderr, "Solver failure, stopping integration\n");
      break;
    }
  }
  printf("   -------------------------\n");

  /* Get/print some final statistics on how the solve progressed */
  flag = FIRKodeGetNumSteps(firkode_mem, &nst);
  check_flag(&flag, "FIRKodeGetNumSteps", 1);
  flag = FIRKodeGetNumStepAttempts(firkode_mem, &nst_a);
  check_flag(&flag, "FIRKodeGetNumStepAttempts", 1);
  flag = FIRKodeGetNumRhsEvals(firkode_mem, &nfe);
  check_flag(&flag, "FIRKodeGetNumRhsEvals", 1);
  flag = FIRKodeGetNumLinSolvSetups(firkode_mem, &nsetups);
  check_flag(&flag, "FIRKodeGetNumLinSolvSetups", 1);
  flag = FIRKodeGetNumErrTestFails(firkode_mem, &netf);
  check_flag(&flag, "FIRKodeGetNumErrTestFails", 1);
  flag = FIRKodeGetNumStepSolveFails(firkode_mem, &ncfn);
  check_flag(&flag, "FIRKodeGetNumStepSolveFails", 1);
  flag = FIRKodeGetNumNonlinSolvIters(firkode_mem, &nni);
  check_flag(&flag, "FIRKodeGetNumNonlinSolvIters", 1);
  flag = FIRKodeGetNumNonlinSolvConvFails(firkode_mem, &nnf);
  check_flag(&flag, "FIRKodeGetNumNonlinSolvConvFails", 1);
  flag = FIRKodeGetNumJacEvals(firkode_mem, &nje);
  check_flag(&flag, "FIRKodeGetNumJacEvals", 1);
  flag = FIRKodeGetNumMassSetups(firkode_mem, &nmset);
  check_flag(&flag, "FIRKodeGetNumMassSetups", 1);
  flag = FIRKodeGetNumMassSolves(firkode_mem, &nms);
  check_flag(&flag, "FIRKodeGetNumMassSolves", 1);
  flag = FIRKodeGetNumMassMult(firkode_mem, &nMv);
  check_flag(&flag, "FIRKodeGetNumMassMult", 1);

  printf("\nFinal Solver Statistics:\n");
  printf("   Internal solver steps = %li (attempted = %li)\n", nst, nst_a);
  printf("   Total RHS evals = %li\n", nfe);
  printf("   Total linear solver setups = %li\n", nsetups);
  printf("   Total number of Jacobian evaluations = %li\n", nje);
  printf("   Total number of Newton iterations = %li\n", nni);
  printf("   Total number of nonlinear solver convergence failures = %li\n", nnf);
  printf("   Total number of error test failures = %li\n", netf);
  printf("   Total number of failed steps from solver failure = %li\n", ncfn);
  printf("   Total mass matrix setups = %li\n", nmset);
  printf("   Total mass matrix solves = %li\n", nms);
  printf("   Total mass times evals = %li\n", nMv);

  /* Clean up and return with successful completion */
  N_VDestroy(y);             /* Free y vector */
  free(udata);               /* Free user data */
  FIRKodeFree(&firkode_mem); /* Free integrator memory */
  SUNLinSolFree(LS);         /* Free linear solvers */
  SUNLinSolFree(MLS);
  SUNMatDestroy(A); /* Free matrices */
  SUNMatDestroy(M);
  SUNContext_Free(&ctx); /* Free context */

  return 0;
}

/*--------------------------------
 * Functions called by the solver
 *--------------------------------*/

/* f routine to compute the right-hand side of M u' = f(t,u):
   f = -k K u + M fsrc(t), with zero boundary rows */
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  UserData udata    = (UserData)user_data; /* access problem data */
  sunindextype N    = udata->N;            /* set variable shortcuts */
  sunrealtype k     = udata->k;
  sunrealtype dx    = udata->dx;
  sunrealtype* Y    = N_VGetArrayPointer(y); /* access data arrays */
  sunrealtype* Ydot = N_VGetArrayPointer(ydot);
  sunrealtype c, fl, fc, fr;
  sunindextype i;

  /* source amplitude: f(t,x) = c sin(pi x) */
  c = (k * PI * PI - ONE) * exp(-t);

  /* boundary rows */
  Ydot[0]     = ZERO;
  Ydot[N - 1] = ZERO;

  /* interior rows: stiffness and mass-weighted source */
  for (i = 1; i < N - 1; i++)
  {
    fl      = c * sin(PI * dx * (sunrealtype)(i - 1));
    fc      = c * sin(PI * dx * (sunrealtype)i);
    fr      = c * sin(PI * dx * (sunrealtype)(i + 1));
    Ydot[i] = -k * (-Y[i - 1] + TWO * Y[i] - Y[i + 1]) / dx +
              dx * (fl + FOUR * fc + fr) / SIX;
  }

  return 0; /* Return with success */
}

/* Jacobian routine to compute J(t,u) = -k K (band, mu = ml = 1) */
static int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData udata = (UserData)user_data; /* access problem data */
  sunindextype N = udata->N;            /* set variable shortcuts */
  sunrealtype k  = udata->k;
  sunrealtype dx = udata->dx;
  sunindextype i;

  SUNMatZero(J); /* initialize Jacobian to zero (boundary rows stay zero) */

  for (i = 1; i < N - 1; i++)
  {
    SM_ELEMENT_B(J, i, i - 1) = k / dx;
    SM_ELEMENT_B(J, i, i)     = -TWO * k / dx;
    SM_ELEMENT_B(J, i, i + 1) = k / dx;
  }

  return 0; /* Return with success */
}

/* Routine to compute the mass matrix M (band, mu = ml = 1) */
static int MassMatrix(sunrealtype t, SUNMatrix M, void* user_data,
                      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData udata = (UserData)user_data; /* access problem data */
  sunindextype N = udata->N;            /* set variable shortcuts */
  sunrealtype dx = udata->dx;
  sunindextype i;

  SUNMatZero(M); /* initialize mass matrix to zero */

  /* boundary rows: identity */
  SM_ELEMENT_B(M, 0, 0)         = ONE;
  SM_ELEMENT_B(M, N - 1, N - 1) = ONE;

  /* interior rows: P1 mass matrix */
  for (i = 1; i < N - 1; i++)
  {
    SM_ELEMENT_B(M, i, i - 1) = dx / SIX;
    SM_ELEMENT_B(M, i, i)     = TWO * dx / SUN_RCONST(3.0);
    SM_ELEMENT_B(M, i, i + 1) = dx / SIX;
  }

  return 0; /* Return with success */
}

/*-------------------------------
 * Private helper functions
 *-------------------------------*/

/* exact solution u(t,x) = exp(-t) sin(pi x) */
static sunrealtype exact(sunrealtype t, sunrealtype x)
{
  return exp(-t) * sin(PI * x);
}

/* maximum nodal error with respect to the exact solution */
static sunrealtype max_error(sunrealtype t, N_Vector y, UserData udata)
{
  sunrealtype* Y  = N_VGetArrayPointer(y);
  sunrealtype err = ZERO, e;
  sunindextype i;

  for (i = 0; i < udata->N; i++)
  {
    e = fabs(Y[i] - exact(t, udata->dx * (sunrealtype)i));
    if (e > err) { err = e; }
  }
  return err;
}

/* Check function return value...
    opt == 0 means SUNDIALS function allocates memory so check if
             returned NULL pointer
    opt == 1 means SUNDIALS function returns a flag so check if
             flag >= 0
    opt == 2 means function allocates memory so check if returned
             NULL pointer
*/
static int check_flag(void* flagvalue, const char* funcname, int opt)
{
  int* errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL)
  {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  /* Check if flag < 0 */
  else if (opt == 1)
  {
    errflag = (int*)flagvalue;
    if (*errflag < 0)
    {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
              funcname, *errflag);
      return 1;
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL)
  {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  return 0;
}

/*---- end of file ----*/
