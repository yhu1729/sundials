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
 * Example problem:
 *
 * The following test simulates the Robertson problem,
 * corresponding to the kinetics of an autocatalytic reaction.
 * This is an ODE system with 3 components, Y = [u,v,w], satisfying
 * the equations,
 *    du/dt = -0.04*u + 1e4*v*w
 *    dv/dt = 0.04*u - 1e4*v*w - 3e7*v^2
 *    dw/dt = 3e7*v^2
 * for t in the interval [0.0, 1e11], with initial conditions
 * Y0 = [1,0,0].
 *
 * While integrating the system, we use the rootfinding feature
 * to find the times at which either u=1e-4 or w=1e-2.
 *
 * This program solves the problem with the 3-stage Radau IIA method
 * of FIRKODE, using the dense SUNLinearSolver and a user-supplied
 * Jacobian routine.
 *
 * 12 outputs are printed at times 0.4*10^k, and run statistics are
 * printed at the end.
 *---------------------------------------------------------------*/

/* Header files */
#include <firkode/firkode.h>    /* prototypes for FIRKODE fcts., consts */
#include <firkode/firkode_ls.h> /* FIRKODE linear solver interface      */
#include <math.h>
#include <nvector/nvector_serial.h> /* serial N_Vector types, fcts., macros */
#include <stdio.h>
#include <sundials/sundials_types.h> /* defs. of 'sunrealtype', 'sunindextype' */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver     */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix           */

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#else
#define GSYM "g"
#define ESYM "e"
#endif

/* User-supplied Functions Called by the Solver */
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int g(sunrealtype t, N_Vector y, sunrealtype* gout, void* user_data);

/* Private function to check function return values */
static int check_flag(void* flagvalue, const char* funcname, int opt);

/* Main Program */
int main(void)
{
  /* general problem parameters */
  sunrealtype T0    = SUN_RCONST(0.0);  /* initial time */
  sunrealtype T1    = SUN_RCONST(0.4);  /* first output time */
  sunrealtype TMult = SUN_RCONST(10.0); /* output time multiplication factor */
  int Nt            = 12;               /* total number of output times */
  sunindextype NEQ  = 3;                /* number of dependent vars. */

  /* rootfinding variables */
  int rootsfound[2];
  int rtflag; /* root info flag */

  /* general problem variables */
  int flag;                  /* reusable error-checking flag */
  N_Vector y         = NULL; /* empty vector for storing solution */
  N_Vector atols     = NULL; /* empty vector for absolute tolerances */
  SUNMatrix A        = NULL; /* empty matrix for linear solver */
  SUNLinearSolver LS = NULL; /* empty linear solver object */
  void* firkode_mem  = NULL; /* empty FIRKODE memory structure */
  FILE* UFID;
  sunrealtype t, tout;
  int iout;
  long int nst, nst_a, nfe, nsetups, nje, nni, nnf, ncfn, netf, nge;

  /* set up the initial conditions */
  sunrealtype u0     = SUN_RCONST(1.0);
  sunrealtype v0     = SUN_RCONST(0.0);
  sunrealtype w0     = SUN_RCONST(0.0);
  sunrealtype reltol = SUN_RCONST(1.0e-4);

  /* Create the SUNDIALS context object for this simulation */
  SUNContext ctx;
  flag = SUNContext_Create(SUN_COMM_NULL, &ctx);
  if (check_flag(&flag, "SUNContext_Create", 1)) { return 1; }

  /* Initial problem output */
  printf("\nRobertson ODE test problem (with rootfinding):\n");
  printf("    initial conditions:  u0 = %" GSYM ",  v0 = %" GSYM
         ",  w0 = %" GSYM "\n",
         u0, v0, w0);

  /* Initialize data structures */
  y = N_VNew_Serial(NEQ, ctx); /* Create serial vector for solution */
  if (check_flag((void*)y, "N_VNew_Serial", 0)) { return 1; }
  NV_Ith_S(y, 0) = u0; /* Set initial conditions into y */
  NV_Ith_S(y, 1) = v0;
  NV_Ith_S(y, 2) = w0;

  atols = N_VClone(y); /* Create serial vector absolute tolerances */
  if (check_flag((void*)atols, "N_VNew_Serial", 0)) { return 1; }

  /* Set absolute tolerances */
  NV_Ith_S(atols, 0) = SUN_RCONST(1.0e-8);
  NV_Ith_S(atols, 1) = SUN_RCONST(1.0e-11);
  NV_Ith_S(atols, 2) = SUN_RCONST(1.0e-8);

  /* Call FIRKodeCreate and FIRKodeInit to initialize the integrator and
     specify the right-hand side function in y'=f(t,y), the initial time
     T0, and the initial dependent variable vector y. */
  firkode_mem = FIRKodeCreate(ctx);
  if (check_flag((void*)firkode_mem, "FIRKodeCreate", 0)) { return 1; }
  flag = FIRKodeInit(firkode_mem, f, T0, y);
  if (check_flag(&flag, "FIRKodeInit", 1)) { return 1; }

  /* Set routines */
  flag = FIRKodeSetMaxNumSteps(firkode_mem, 100000); /* Increase max num steps */
  if (check_flag(&flag, "FIRKodeSetMaxNumSteps", 1)) { return 1; }
  flag = FIRKodeSVtolerances(firkode_mem, reltol, atols); /* Specify tolerances */
  if (check_flag(&flag, "FIRKodeSVtolerances", 1)) { return 1; }

  /* Specify the root-finding function, having 2 equations */
  flag = FIRKodeRootInit(firkode_mem, 2, g);
  if (check_flag(&flag, "FIRKodeRootInit", 1)) { return 1; }

  /* Initialize dense matrix data structure and solver */
  A = SUNDenseMatrix(NEQ, NEQ, ctx);
  if (check_flag((void*)A, "SUNDenseMatrix", 0)) { return 1; }
  LS = SUNLinSol_Dense(y, A, ctx);
  if (check_flag((void*)LS, "SUNLinSol_Dense", 0)) { return 1; }

  /* Linear solver interface */
  flag = FIRKodeSetLinearSolver(firkode_mem, LS, A); /* Attach matrix and LS */
  if (check_flag(&flag, "FIRKodeSetLinearSolver", 1)) { return 1; }
  flag = FIRKodeSetJacFn(firkode_mem, Jac); /* Set the Jacobian routine */
  if (check_flag(&flag, "FIRKodeSetJacFn", 1)) { return 1; }

  /* Open output stream for results, output comment line */
  UFID = fopen("solution.txt", "w");
  fprintf(UFID, "# t u v w\n");

  /* output initial condition to disk */
  fprintf(UFID, " %.16" ESYM " %.16" ESYM " %.16" ESYM " %.16" ESYM "\n", T0,
          NV_Ith_S(y, 0), NV_Ith_S(y, 1), NV_Ith_S(y, 2));

  /* Main time-stepping loop: calls FIRKodeEvolve to perform the integration,
     then prints results.  Stops when the final time has been reached */
  t = T0;
  printf("        t             u             v             w\n");
  printf("   -----------------------------------------------------\n");
  printf("  %12.5" ESYM "  %12.5" ESYM "  %12.5" ESYM "  %12.5" ESYM "\n", t,
         NV_Ith_S(y, 0), NV_Ith_S(y, 1), NV_Ith_S(y, 2));
  tout = T1;
  iout = 0;
  while (1)
  {
    flag = FIRKodeEvolve(firkode_mem, tout, y, &t,
                         FIRK_NORMAL); /* call integrator */
    if (check_flag(&flag, "FIRKodeEvolve", 1)) { break; }
    printf("  %12.5" ESYM "  %12.5" ESYM "  %12.5" ESYM "  %12.5" ESYM
           "\n", /* access/print solution */
           t, NV_Ith_S(y, 0), NV_Ith_S(y, 1), NV_Ith_S(y, 2));
    fprintf(UFID, " %.16" ESYM " %.16" ESYM " %.16" ESYM " %.16" ESYM "\n", t,
            NV_Ith_S(y, 0), NV_Ith_S(y, 1), NV_Ith_S(y, 2));
    if (flag == FIRK_ROOT_RETURN)
    { /* check if a root was found */
      rtflag = FIRKodeGetRootInfo(firkode_mem, rootsfound);
      if (check_flag(&rtflag, "FIRKodeGetRootInfo", 1)) { return 1; }
      printf("      rootsfound[] = %3d %3d\n", rootsfound[0], rootsfound[1]);
    }
    if (flag >= 0)
    { /* successful solve: update output time */
      iout++;
      tout *= TMult;
    }
    else
    { /* unsuccessful solve: break */
      fprintf(stderr, "Solver failure, stopping integration\n");
      break;
    }
    if (iout == Nt) { break; /* stop after enough outputs */ }
  }
  printf("   -----------------------------------------------------\n");
  fclose(UFID);

  /* Print some final statistics */
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
  flag = FIRKodeGetNumGEvals(firkode_mem, &nge);
  check_flag(&flag, "FIRKodeGetNumGEvals", 1);

  printf("\nFinal Solver Statistics:\n");
  printf("   Internal solver steps = %li (attempted = %li)\n", nst, nst_a);
  printf("   Total RHS evals = %li\n", nfe);
  printf("   Total linear solver setups = %li\n", nsetups);
  printf("   Total number of Jacobian evaluations = %li\n", nje);
  printf("   Total number of Newton iterations = %li\n", nni);
  printf("   Total root-function g evals = %li\n", nge);
  printf("   Total number of nonlinear solver convergence failures = %li\n", nnf);
  printf("   Total number of error test failures = %li\n", netf);
  printf("   Total number of failed steps from solver failure = %li\n", ncfn);

  /* Clean up and return with successful completion */
  N_VDestroy(y);             /* Free y vector */
  N_VDestroy(atols);         /* Free atols vector */
  FIRKodeFree(&firkode_mem); /* Free integrator memory */
  SUNLinSolFree(LS);         /* Free linear solver */
  SUNMatDestroy(A);          /* Free A matrix */
  SUNContext_Free(&ctx);     /* Free context */

  return 0;
}

/*-------------------------------
 * Functions called by the solver
 *-------------------------------*/

/* f routine to compute the ODE RHS function f(t,y). */
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype u = NV_Ith_S(y, 0); /* access current solution */
  sunrealtype v = NV_Ith_S(y, 1);
  sunrealtype w = NV_Ith_S(y, 2);

  /* Fill in ODE RHS function */
  NV_Ith_S(ydot, 0) = -0.04 * u + 1.e4 * v * w;
  NV_Ith_S(ydot, 1) = 0.04 * u - 1.e4 * v * w - 3.e7 * v * v;
  NV_Ith_S(ydot, 2) = 3.e7 * v * v;

  return 0; /* Return with success */
}

/* Jacobian routine to compute J(t,y) = df/dy. */
static int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  sunrealtype v = NV_Ith_S(y, 1); /* access current solution */
  sunrealtype w = NV_Ith_S(y, 2);
  SUNMatZero(J); /* initialize Jacobian to zero */

  /* Fill in the Jacobian of the ODE RHS function */
  SM_ELEMENT_D(J, 0, 0) = -0.04;
  SM_ELEMENT_D(J, 0, 1) = 1.e4 * w;
  SM_ELEMENT_D(J, 0, 2) = 1.e4 * v;

  SM_ELEMENT_D(J, 1, 0) = 0.04;
  SM_ELEMENT_D(J, 1, 1) = -1.e4 * w - 6.e7 * v;
  SM_ELEMENT_D(J, 1, 2) = -1.e4 * v;

  SM_ELEMENT_D(J, 2, 1) = 6.e7 * v;

  return 0; /* Return with success */
}

/* g routine to compute the root-finding function g(t,y). */
static int g(sunrealtype t, N_Vector y, sunrealtype* gout, void* user_data)
{
  sunrealtype u = NV_Ith_S(y, 0); /* access current solution */
  sunrealtype w = NV_Ith_S(y, 2);

  gout[0] = u - SUN_RCONST(0.0001); /* check for u == 1e-4 */
  gout[1] = w - SUN_RCONST(0.01);   /* check for w == 1e-2 */

  return 0; /* Return with success */
}

/*-------------------------------
 * Private helper functions
 *-------------------------------*/

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
