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
 * The following is a simple example problem with analytical
 * solution,
 *    dy/dt = lambda*y + 1/(1+t^2) - lambda*atan(t)
 * for t in the interval [0.0, 10.0], with initial condition: y=0.
 *
 * The stiffness of the problem is directly proportional to the
 * value of "lambda".  The value of lambda should be negative to
 * result in a well-posed ODE; for values with magnitude larger
 * than 100 the problem becomes quite stiff.
 *
 * This program solves the problem with the s-stage Radau IIA
 * method of FIRKODE (s = 3 by default, or given as the first
 * command-line argument), using the dense linear solver and a
 * user-supplied Jacobian. Output is printed every 1.0 units of
 * time (10 total) together with the error against the analytical
 * solution y(t) = atan(t). Run statistics (optional outputs) are
 * printed at the end.
 * -----------------------------------------------------------------*/

/* Header files */
#include <firkode/firkode.h>    /* prototypes for FIRKODE fcts., consts */
#include <firkode/firkode_ls.h> /* FIRKODE linear solver interface      */
#include <math.h>
#include <nvector/nvector_serial.h> /* serial N_Vector types, fcts., macros */
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_types.h>   /* definition of type sunrealtype    */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix       */

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

/* User-supplied Functions Called by the Solver */
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private function to check function return values */
static int check_flag(void* flagvalue, const char* funcname, int opt);

/* Main Program */
int main(int argc, char* argv[])
{
  /* general problem parameters */
  sunrealtype T0     = SUN_RCONST(0.0);    /* initial time */
  sunrealtype Tf     = SUN_RCONST(10.0);   /* final time */
  sunrealtype dTout  = SUN_RCONST(1.0);    /* time between outputs */
  sunindextype NEQ   = 1;                  /* number of dependent vars. */
  sunrealtype reltol = SUN_RCONST(1.0e-6); /* tolerances */
  sunrealtype abstol = SUN_RCONST(1.0e-10);
  sunrealtype lambda = SUN_RCONST(-100.0); /* stiffness parameter */
  int stages         = 3;                  /* number of Radau IIA stages */

  /* general problem variables */
  int flag;                  /* reusable error-checking flag */
  N_Vector y         = NULL; /* empty vector for storing solution */
  SUNMatrix A        = NULL; /* empty matrix for linear solver */
  SUNLinearSolver LS = NULL; /* empty linear solver object */
  void* firkode_mem  = NULL; /* empty FIRKODE memory structure */
  FILE* UFID;
  sunrealtype t, tout, err, maxerr;
  long int nst, nst_a, nfe, nsetups, nje, nni, ncfn, netf, nsli;

  /* optional command-line inputs: number of stages and relative tolerance */
  if (argc > 1) { stages = atoi(argv[1]); }
  if (argc > 2) { reltol = (sunrealtype)atof(argv[2]); }

  /* Create the SUNDIALS context object for this simulation */
  SUNContext ctx;
  flag = SUNContext_Create(SUN_COMM_NULL, &ctx);
  if (check_flag(&flag, "SUNContext_Create", 1)) { return 1; }

  /* Initial diagnostics output */
  printf("\nAnalytical ODE test problem:\n");
  printf("    lambda = %" GSYM "\n", lambda);
  printf("    stages = %i (Radau IIA of order %i)\n", stages, 2 * stages - 1);
  printf("    reltol = %.1" ESYM "\n", reltol);
  printf("    abstol = %.1" ESYM "\n\n", abstol);

  /* Initialize data structures */
  y = N_VNew_Serial(NEQ, ctx); /* Create serial vector for solution */
  if (check_flag((void*)y, "N_VNew_Serial", 0)) { return 1; }
  N_VConst(SUN_RCONST(0.0), y); /* Specify initial condition */

  /* Call FIRKodeCreate and FIRKodeInit to initialize the integrator memory
     and specify the right-hand side function, the initial time T0, and
     the initial dependent variable vector y. */
  firkode_mem = FIRKodeCreate(ctx);
  if (check_flag((void*)firkode_mem, "FIRKodeCreate", 0)) { return 1; }
  flag = FIRKodeInit(firkode_mem, f, T0, y);
  if (check_flag(&flag, "FIRKodeInit", 1)) { return 1; }

  /* Set routines */
  flag = FIRKodeSetUserData(firkode_mem, (void*)&lambda);
  if (check_flag(&flag, "FIRKodeSetUserData", 1)) { return 1; }
  flag = FIRKodeSetNumStages(firkode_mem, stages);
  if (check_flag(&flag, "FIRKodeSetNumStages", 1)) { return 1; }
  flag = FIRKodeSStolerances(firkode_mem, reltol, abstol);
  if (check_flag(&flag, "FIRKodeSStolerances", 1)) { return 1; }
  flag = FIRKodeSetMaxNumSteps(firkode_mem, 100000);
  if (check_flag(&flag, "FIRKodeSetMaxNumSteps", 1)) { return 1; }

  /* Initialize dense matrix data structure and solver */
  A = SUNDenseMatrix(NEQ, NEQ, ctx);
  if (check_flag((void*)A, "SUNDenseMatrix", 0)) { return 1; }
  LS = SUNLinSol_Dense(y, A, ctx);
  if (check_flag((void*)LS, "SUNLinSol_Dense", 0)) { return 1; }

  /* Attach the linear solver and matrix, and set the Jacobian routine */
  flag = FIRKodeSetLinearSolver(firkode_mem, LS, A);
  if (check_flag(&flag, "FIRKodeSetLinearSolver", 1)) { return 1; }
  flag = FIRKodeSetJacFn(firkode_mem, Jac);
  if (check_flag(&flag, "FIRKodeSetJacFn", 1)) { return 1; }

  /* Open output stream for results, output comment line */
  UFID = fopen("solution.txt", "w");
  fprintf(UFID, "# t u\n");

  /* output initial condition to disk */
  fprintf(UFID, " %.16" ESYM " %.16" ESYM "\n", T0, NV_Ith_S(y, 0));

  /* Main time-stepping loop: calls FIRKodeEvolve to perform the integration,
     then prints results. Stops when the final time has been reached. */
  t      = T0;
  tout   = T0 + dTout;
  maxerr = SUN_RCONST(0.0);
  printf("        t           u          error\n");
  printf("   ---------------------------------------\n");
  while (Tf - t > SUN_RCONST(1.0e-15))
  {
    flag = FIRKodeEvolve(firkode_mem, tout, y, &t,
                         FIRK_NORMAL); /* call integrator */
    if (check_flag(&flag, "FIRKodeEvolve", 1)) { break; }
    err    = fabs(NV_Ith_S(y, 0) - atan(t));
    maxerr = (err > maxerr) ? err : maxerr;
    printf("  %10.6" FSYM "  %10.6" FSYM "  %10.3" ESYM "\n", t, NV_Ith_S(y, 0),
           err); /* access/print solution */
    fprintf(UFID, " %.16" ESYM " %.16" ESYM "\n", t, NV_Ith_S(y, 0));
    if (flag >= 0)
    { /* successful solve: update time */
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    }
    else
    { /* unsuccessful solve: break */
      fprintf(stderr, "Solver failure, stopping integration\n");
      break;
    }
  }
  printf("   ---------------------------------------\n");
  printf("  Max error = %.3" ESYM "\n", maxerr);
  fclose(UFID);

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
  flag = FIRKodeGetNumNonlinSolvIters(firkode_mem, &nni);
  check_flag(&flag, "FIRKodeGetNumNonlinSolvIters", 1);
  flag = FIRKodeGetNumNonlinSolvConvFails(firkode_mem, &ncfn);
  check_flag(&flag, "FIRKodeGetNumNonlinSolvConvFails", 1);
  flag = FIRKodeGetNumJacEvals(firkode_mem, &nje);
  check_flag(&flag, "FIRKodeGetNumJacEvals", 1);
  flag = FIRKodeGetNumStageLinIters(firkode_mem, &nsli);
  check_flag(&flag, "FIRKodeGetNumStageLinIters", 1);

  printf("\nFinal Solver Statistics:\n");
  printf("   Internal solver steps = %li (attempted = %li)\n", nst, nst_a);
  printf("   Total RHS evals = %li\n", nfe);
  printf("   Total linear solver setups = %li\n", nsetups);
  printf("   Total number of Jacobian evaluations = %li\n", nje);
  printf("   Total number of Newton iterations = %li\n", nni);
  printf("   Total number of stage linear iterations = %li\n", nsli);
  printf("   Total number of nonlinear solver convergence failures = %li\n",
         ncfn);
  printf("   Total number of error test failures = %li\n\n", netf);

  /* Clean up and return with successful completion */
  N_VDestroy(y);             /* Free y vector */
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
  sunrealtype* rdata = (sunrealtype*)user_data; /* cast user_data to sunrealtype */
  sunrealtype lambda = rdata[0]; /* set shortcut for stiffness parameter */
  sunrealtype u      = NV_Ith_S(y, 0); /* access current solution value */

  /* fill in the RHS function: "NV_Ith_S" accesses the 0th entry of ydot */
  NV_Ith_S(ydot, 0) = lambda * u + SUN_RCONST(1.0) / (SUN_RCONST(1.0) + t * t) -
                      lambda * atan(t);

  return 0; /* return with success */
}

/* Jacobian routine to compute J(t,y) = df/dy. */
static int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  sunrealtype* rdata = (sunrealtype*)user_data; /* cast user_data to sunrealtype */
  sunrealtype lambda = rdata[0]; /* set shortcut for stiffness parameter */
  sunrealtype* Jdata = SUNDenseMatrix_Data(J);

  /* Fill in Jacobian of f: set the first entry of the data array to set the (0,0) entry */
  Jdata[0] = lambda;

  return 0; /* return with success */
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
