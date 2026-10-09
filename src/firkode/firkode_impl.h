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
 * Implementation header file for the main FIRKODE integrator.
 * -----------------------------------------------------------------*/

#ifndef _FIRKODE_IMPL_H
#define _FIRKODE_IMPL_H

#include <stdarg.h>

#include <firkode/firkode.h>
#include <firkode/firkode_ls.h>
#include <sundials/priv/sundials_context_impl.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_math.h>

#include "sundials_logger_impl.h"
#include "sundials_macros.h"

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*===============================================================
  SHORTCUTS
  ===============================================================*/

#define FIRK_PROFILER firk_mem->sunctx->profiler
#define FIRK_LOGGER   firk_mem->sunctx->logger

/*===============================================================
  INTERNAL CONSTANTS
  ===============================================================*/

/* tolerance type */
#define FIRK_NN 0 /* no tolerances set          */
#define FIRK_SS 1 /* scalar rel, scalar abs     */
#define FIRK_SV 2 /* scalar rel, vector abs     */
#define FIRK_WF 3 /* user-supplied weight fn    */

/* the only number of stages that FIRKodeSetNumStages and FIRKodeSetOrder
   accept (design decision D7) */
#define FIRK_SUPPORTED_STAGES 3

/* integrator defaults */
#define HMIN_DEFAULT     SUN_RCONST(0.0)
#define HMAX_INV_DEFAULT SUN_RCONST(0.0)
#define MXHNIL_DEFAULT   10
#define MXSTEP_DEFAULT   500
#define MXNEF            7
#define MXNCF            10

/* linear solver setup heuristics */
#define MSBP_DEFAULT         20
#define DGMAX_LSETUP_DEFAULT SUN_RCONST(0.3)

/* Newton iteration defaults (RADAU5 conventions, design Table 5.9); the
   stopping tolerance eps_nls is RADAU5's rule
     max(NLSCOEF_UROUND_FAC u / reltol, min(NLSCOEF_MAX, sqrt(reltol)))
   unless the user fixes a value (firkSetNlsCoef) */
#define MAXCOR_DEFAULT     7
#define NLSCOEF_AUTO       SUN_RCONST(0.0)  /* nlscoef_user <= 0: the rule */
#define NLSCOEF_MAX        SUN_RCONST(0.03) /* cap; value without a reltol */
#define NLSCOEF_UROUND_FAC SUN_RCONST(10.0) /* roundoff floor 10 u / reltol */
#define THETA_MAX_DEFAULT  SUN_RCONST(0.99)
#define THETA_JBAD_DEFAULT SUN_RCONST(1.0) /* >= 1 disables the test */
#define THETA_JAC_DEFAULT  SUN_RCONST(1.0e-3)
#define RATE_EXP_DEFAULT   SUN_RCONST(0.8)
#define QNEWT_MIN          SUN_RCONST(1.0e-4)
#define QNEWT_MAX          SUN_RCONST(20.0)
#define PREDICT_FAIL_FAC   SUN_RCONST(0.8)
#define HFAC_NEWTON_FAIL   SUN_RCONST(0.5) /* divergence, iteration limit */

/* error estimate floor */
#define DSM_FLOOR SUN_RCONST(1.0e-10)

/* RADAU5's smallest step relative to |tn| (in units of the roundoff) */
#define HMIN_REL_ROUNDOFF SUN_RCONST(10.0)

/* step adaptivity defaults */
#define SAFETY    SUN_RCONST(0.9)
#define GROWTH    SUN_RCONST(8.0)
#define ETAMX1    SUN_RCONST(10000.0)
#define ETAMXF    SUN_RCONST(0.3)
#define ETAMIN    SUN_RCONST(0.2)
#define ETACF     SUN_RCONST(0.25)
#define SMALL_NEF 2
#define HFIXED_LB SUN_RCONST(1.0)
#define HFIXED_UB SUN_RCONST(1.2)

/* initial step size heuristics */
#define HLB_FACTOR SUN_RCONST(100.0)
#define HUB_FACTOR SUN_RCONST(0.1)
#define H_BIAS     SUN_RCONST(0.5)
#define MAX_ITERS  4

/* roundoff-related constants */
#define ONEPSM      SUN_RCONST(1.000001)
#define ONEMSM      SUN_RCONST(0.999999)
#define FUZZ_FACTOR SUN_RCONST(100.0)

/* control constants for firkStep (design Table 8.2, private) */
#define PREDICT_AGAIN  +3
#define TRY_AGAIN      +5
#define FIRST_CALL     +6
#define PREV_CONV_FAIL +7
#define PREV_ERR_FAIL  +8
#define RHSFUNC_RECVR  +10

/* kinds of the stages' stack (firkode_stack.c, design Section 8.3) */
#define FIRK_STACK_PLAIN 0 /* NVECTOR_MANYVECTOR    */
#define FIRK_STACK_MPI   1 /* NVECTOR_MPIMANYVECTOR */

/* rootfinding return values */
#define RTFOUND +1
#define CLOSERT +3

/* purposes of a solve with the block matrix M - gamma J (design Listing
   8.4); they differ in the tolerance given to an iterative block solver
   and in how its failures are treated */
typedef enum
{
  FIRKLS_BLOCK_NEWTON,   /* Newton system of the single-stage method */
  FIRKLS_BLOCK_ESTIMATE, /* filtered error estimate */
  FIRKLS_BLOCK_PRECOND   /* preconditioner application in the stacked solve */
} firkLsBlockSolve;

/* statuses of the Newton iteration (design Listing 8.3); the failures
   are positive, and negative values are FIRKODE return codes */
#define FIRK_NLS_SUCCESS        0
#define FIRK_NLS_MAXITER        11
#define FIRK_NLS_DIVERGED       12
#define FIRK_NLS_PREDICTED_FAIL 13
#define FIRK_NLS_RHS_RECVR      14
#define FIRK_NLS_LSOLVE_RECVR   15
#define FIRK_NLS_LSETUP_RECVR   16

/*===============================================================
  TIME STEP ADAPTIVITY MEMORY
  ===============================================================*/

typedef struct FIRKodeHAdaptMemRec
{
  sunrealtype etamax;         /* eta <= etamax                              */
  sunrealtype etamx1;         /* max step size change on first step         */
  sunrealtype etamxf;         /* h reduction factor on multiple error fails */
  sunrealtype etamin;         /* eta >= etamin on error test fail           */
  int small_nef;              /* bound to determine 'multiple' above        */
  sunrealtype etacf;          /* h reduction factor on nonlinear conv fail  */
  sunrealtype safety;         /* accuracy safety factor on h                */
  sunrealtype growth;         /* maximum step growth safety factor          */
  sunrealtype lbound;         /* eta lower bound to leave h unchanged       */
  sunrealtype ubound;         /* eta upper bound to leave h unchanged       */
  int p;                      /* error estimate order                       */
  int q;                      /* method order                               */
  sunbooleantype newt_safety; /* scale safety by Newton iteration count */

  SUNAdaptController hcontroller; /* temporal error controller        */
  sunbooleantype owncontroller;   /* flag indicating ownership        */
}* FIRKodeHAdaptMem;

/*===============================================================
  MAIN INTEGRATOR MEMORY BLOCK
  ===============================================================*/

typedef struct FIRKodeMemRec
{
  SUNContext sunctx;
  sunrealtype uround; /* machine unit roundoff */

  /*--------------------------
    Problem Specification Data
    --------------------------*/

  FIRKRhsFn f;     /* M y' = f(t,y(t))                          */
  void* user_data; /* user pointer passed to f                  */
  int itol;        /* FIRK_SS, FIRK_SV, FIRK_WF, or FIRK_NN     */

  sunrealtype reltol;       /* relative tolerance                */
  sunrealtype Sabstol;      /* scalar absolute tolerance         */
  N_Vector Vabstol;         /* vector absolute tolerance         */
  sunbooleantype atolmin0;  /* flag indicating min(abstol) = 0   */
  sunbooleantype user_efun; /* SUNTRUE if user sets efun         */
  FIRKEwtFn efun;           /* function to set ewt               */
  void* e_data;             /* user pointer passed to efun       */

  /*-------
    Method
    -------*/

  int mode;    /* FIRK_BEGINNER or FIRK_EXPERT                */
  int s;       /* requested number of stages                 */
  int s_alloc; /* number of stages the stage memory holds    */
  const struct FIRKodeTableMem* tbl; /* built-in Radau IIA table for s stages */
  int predictor;           /* stage predictor method                     */
  sunbooleantype refilter; /* recompute a large error estimate on the first
                              step or after a rejection (RADAU5)          */

  /*-------------------
    Vectors of length N
    -------------------*/

  N_Vector yn;     /* solution at the start of the current step   */
  N_Vector ycur;   /* candidate solution / scratch                */
  N_Vector fn;     /* f(tn, yn)                                   */
  N_Vector ydn;    /* M^{-1} f(tn, yn); aliases fn without mass   */
  N_Vector ewt;    /* error weight vector                         */
  N_Vector ele;    /* filtered local error estimate               */
  N_Vector ftemp;  /* temporary storage vector                    */
  N_Vector tempv1; /* temporary storage vector                    */
  N_Vector tempv2; /* temporary storage vector                    */
  N_Vector tempv3; /* temporary storage vector                    */

  sunbooleantype fn_is_current;  /* fn holds f(tn, yn)            */
  sunbooleantype ydn_is_current; /* ydn holds M^{-1} f(tn, yn)    */

  /*-------------------------------
    Stage storage (arrays of length s)
    -------------------------------*/

  N_Vector* Z;     /* stage increments z_i = Y_i - yn, predicted  */
  N_Vector* Zprev; /* stage increments of the last accepted step  */
  N_Vector* F;     /* stage right-hand sides f(tn + c_i h, Y_i)   */
  N_Vector* G;     /* Newton right-hand side -G_i(Z)              */
  N_Vector* delta; /* Newton update                               */

  int stack_kind;     /* FIRK_STACK_PLAIN or FIRK_STACK_MPI         */
  N_Vector g_stk;     /* stack of G                                */
  N_Vector delta_stk; /* stack of delta                            */
  N_Vector ewt_stk;   /* stack of ewt repeated s times (a view)    */

  sunbooleantype have_prev; /* Zprev, hold describe the last step  */

  /*-----------------------
    Fused vector operations
    -----------------------*/

  sunrealtype* cvals; /* array of scalars, length s + 2           */
  N_Vector* Xvecs;    /* array of vectors, length s + 2           */
  sunrealtype* wtmp;  /* scratch scalars, length s + 2            */

  /*-----------------
    Tstop information
    -----------------*/

  sunbooleantype tstopset;
  sunbooleantype tstopinterp;
  sunrealtype tstop;

  /*---------
    Step Data
    ---------*/

  sunrealtype hin;      /* initial step size                       */
  sunrealtype h;        /* current step size                       */
  sunrealtype hprime;   /* step size to be used on the next step   */
  sunrealtype next_h;   /* step size to be used on the next step   */
  sunrealtype eta;      /* eta = hprime / h                        */
  sunrealtype hold;     /* last successful step size               */
  sunrealtype h0u;      /* actual initial step size                */
  sunrealtype tn;       /* time at the start of the current step   */
  sunrealtype tcur;     /* current internal time                   */
  sunrealtype tretlast; /* last value of t returned by Evolve      */

  sunrealtype gamma;  /* gamma = h * gamma0                        */
  sunrealtype gammap; /* gamma at the last setup call              */
  sunrealtype gamrat; /* gamma / gammap                            */

  sunrealtype dsm; /* last local error test value                  */

  sunbooleantype fixedstep;  /* fixed step mode                    */
  sunbooleantype initsetup;  /* initial setup still pending        */
  sunbooleantype firststage; /* first step since (Re)Init/Reset    */

  /*-------------------
    Newton iteration
    -------------------*/

  int maxcor;               /* maximum Newton iterations m          */
  sunrealtype nlscoef_user; /* eps_nls as set; <= 0: RADAU5's rule  */
  sunrealtype nlscoef;      /* stopping tolerance eps_nls in use    */
  sunrealtype theta_max;    /* divergence threshold on the rate     */
  sunrealtype theta_jbad;   /* rate above which a setup is forced   */
  sunrealtype theta_jac;    /* rate at the start of each solve      */
  sunrealtype rate_exp;     /* exponent of the first-iteration factor */

  sunrealtype eta_nls;       /* theta/(1-theta) at the end of the last solve */
  sunrealtype theta;         /* rate at the end of the last solve        */
  sunrealtype nls_hfac;      /* step factor of the last failure          */
  long int nni_last;         /* iterations of the last Newton solve      */
  sunbooleantype jcur;       /* did the last setup evaluate J?           */
  sunbooleantype keep_setup; /* the dead band kept h (so gamma == gammap
                                at the next attempt): read by tests only */

  /*------------------------------
    Reuse policy (design D9)
    ------------------------------*/

  int reuse_policy;              /* FIRK_POLICY_RADAU5 or _PERIODIC        */
  sunbooleantype jac_current;    /* J was evaluated at the current (tn, yn) */
  sunbooleantype jac_stale_fail; /* J was not current at the last failure */
  sunrealtype dgamma_fail;       /* |gamma/gammap - 1| at the last failure */

  /*------
    Limits
    ------*/

  long int mxstep; /* maximum number of internal steps per call   */
  int mxhnil;      /* maximum number of t + h == t warnings       */
  int maxnef;      /* maximum number of error test failures       */
  int maxncf;      /* maximum number of nonlinear conv failures   */

  sunrealtype hmin;     /* |h| >= hmin                            */
  sunrealtype hmax_inv; /* |h| <= 1/hmax_inv                      */

  long int msbp;            /* max steps between lsetup calls     */
  sunrealtype dgmax_lsetup; /* gamma ratio threshold for lsetup   */
  long int nstlp;           /* step number of last setup call     */

  /*--------
    Counters
    --------*/

  long int nst;          /* number of internal steps taken        */
  long int nst_attempts; /* number of step attempts               */
  long int nfe;          /* number of f calls                     */
  long int ncfn;         /* number of nonlinear conv failures     */
  long int nni;          /* number of nonlinear iterations        */
  long int nnf;          /* number of nonlinear conv failures     */
  long int netf;         /* number of error test failures         */
  long int nsetups;      /* number of setup calls                 */
  int nhnil;             /* number of t + h == t warnings         */
  sunrealtype tolsf;     /* tolerance scale factor                */

  /*-----------
    Adaptivity
    -----------*/

  FIRKodeHAdaptMem hadapt_mem;

  /*------------------
    Linear Solver Data
    ------------------*/

  int (*linit)(struct FIRKodeMemRec* firk_mem);

  int (*lreinit)(struct FIRKodeMemRec* firk_mem);

  int (*lsetup)(struct FIRKodeMemRec* firk_mem, sunbooleantype jbad,
                sunbooleantype* jcur);

  int (*lsolve_stk)(struct FIRKodeMemRec* firk_mem, N_Vector b_stk,
                    N_Vector x_stk, int k);

  int (*lsolve_blk)(struct FIRKodeMemRec* firk_mem, N_Vector b, N_Vector x,
                    firkLsBlockSolve purpose, int k);

  int (*lfree)(struct FIRKodeMemRec* firk_mem);

  void* lmem; /* linear solver interface memory structure */

  /*----------------
    Mass Matrix Data
    ----------------*/

  int (*minit)(struct FIRKodeMemRec* firk_mem);

  int (*mreinit)(struct FIRKodeMemRec* firk_mem);

  int (*msetup)(struct FIRKodeMemRec* firk_mem, sunrealtype t, N_Vector vtemp1,
                N_Vector vtemp2, N_Vector vtemp3);

  int (*mmult)(struct FIRKodeMemRec* firk_mem, N_Vector v, N_Vector Mv);

  int (*msolve)(struct FIRKodeMemRec* firk_mem, N_Vector b, sunrealtype tol);

  int (*mfree)(struct FIRKodeMemRec* firk_mem);

  void* mass_mem;                /* mass matrix interface memory structure */
  sunbooleantype mass_set;       /* a mass matrix has been attached        */
  sunbooleantype mass_iterative; /* the mass linear solver is iterative    */

  /*----------------
    Rootfinding Data
    ----------------*/

  struct FIRKodeRootMemRec* root_mem;

  /*-----
    Flags
    -----*/

  sunbooleantype VabstolMallocDone;
  sunbooleantype MallocDone;

}* FIRKodeMem;

/*===============================================================
  INTERFACE TO LINEAR SOLVERS (design Listing 8.4)

  lsetup(jbad, &jcur) forms M - gamma J at (tn, yn, fn), with a new
  Jacobian if jbad, sets up the block solver, reports in jcur whether
  the Jacobian data is current, and records the setup (gammap,
  nsetups, nstlp) whatever its outcome. The caller decides when to
  set up and whether J is bad; lsetup is NULL when no setup is needed.

  lsolve_stk(b, x, k) solves the coupled stage system
    (I_s (x) M - h A (x) J) x = b
  for the k-th Newton iteration (k >= 1), on stacked (ManyVector)
  vectors.

  lsolve_blk(b, x, purpose, k) solves the block system
  (M - gamma J) x = b for the given purpose; k is the Newton iteration
  for FIRKLS_BLOCK_NEWTON.

  The solves accept x == b and may overwrite b. They return 0 on
  success, a positive value on a recoverable failure, and a negative
  value on an unrecoverable failure, as lsetup does.
  ===============================================================*/

/*===============================================================
  INTERNAL FUNCTION PROTOTYPES
  ===============================================================*/

/* firkode.c */

void firkProcessError(FIRKodeMem firk_mem, int error_code, int line,
                      const char* func, const char* file, const char* msgfmt,
                      ...);

int firkEwtSet(N_Vector ycur, N_Vector weight, void* data);

int firkRhs(FIRKodeMem firk_mem, sunrealtype t, N_Vector y, N_Vector ydot);

void firkSetNlsCoef(FIRKodeMem firk_mem);

/* firkode_stack.c */

int firkStackSetKind(FIRKodeMem firk_mem, N_Vector y0);

N_Vector firkStackNew(FIRKodeMem firk_mem, int s, N_Vector* vecs);

N_Vector firkStackWeightView(FIRKodeMem firk_mem, int s, N_Vector w);

N_Vector firkStage(FIRKodeMem firk_mem, N_Vector stk, int i);
N_Vector* firkCloneVectorArray(int count, N_Vector w);

/* firkode_ls.c: the Jacobian schedule that the PERIODIC policy reads */

void firkLsJacSchedule(FIRKodeMem firk_mem, sunbooleantype* due,
                       sunrealtype* dgmax_jbad);

int firkGetYdot(FIRKodeMem firk_mem);

int firkSetExpertDefaults(FIRKodeMem firk_mem);

/* firkode_io.c */

sunbooleantype firkExpertAllowed(FIRKodeMem firk_mem, const char* func);

int firkOrderToStages(int ord);

int firkSetNumStagesAny(void* firkode_mem, int s);

/* firkode_step.c */

int firkStep(FIRKodeMem firk_mem);

/* firkode_dense.c */

int firkPredict(FIRKodeMem firk_mem);

int firkDenseEval(FIRKodeMem firk_mem, sunrealtype t, int k, N_Vector dky);

/* firkode_adapt.c */

void firkReusePolicy(FIRKodeMem firk_mem, int nflag, sunbooleantype* jbad,
                     sunbooleantype* need_setup);

FIRKodeHAdaptMem firkAdaptInit(void);

sunrealtype firkHminEff(FIRKodeMem firk_mem);

void firkPrintAdaptMem(FIRKodeHAdaptMem hadapt_mem, FILE* outfile);

int firkAdapt(FIRKodeMem firk_mem, sunrealtype dsm, sunbooleantype accepted);

int firkReplaceAdaptController(FIRKodeMem firk_mem, SUNAdaptController C,
                               sunbooleantype take_ownership);

/* firkode_nls.c */

int firkNewton(FIRKodeMem firk_mem);

/*===============================================================
  ERROR MESSAGES
  ===============================================================*/

#define MSG_TIME   "t = " SUN_FORMAT_G
#define MSG_TIME_H "t = " SUN_FORMAT_G " and h = " SUN_FORMAT_G
#define MSG_TIME_INT                                                \
  "t = " SUN_FORMAT_G " is not between tcur - hold = " SUN_FORMAT_G \
  " and tcur = " SUN_FORMAT_G
#define MSG_TIME_TOUT  "tout = " SUN_FORMAT_G
#define MSG_TIME_TSTOP "tstop = " SUN_FORMAT_G

/* Initialization and I/O error messages */

#define MSGFIRK_NO_MEM        "firkode_mem = NULL illegal."
#define MSGFIRK_MEM_FAIL      "A memory request failed."
#define MSGFIRK_NULL_SUNCTX   "sunctx = NULL illegal."
#define MSGFIRK_NO_MALLOC     "Attempt to call before FIRKodeInit."
#define MSGFIRK_NEG_HMIN      "hmin < 0 illegal."
#define MSGFIRK_NEG_HMAX      "hmax < 0 illegal."
#define MSGFIRK_BAD_HMIN_HMAX "Inconsistent step size limits: hmin > hmax."
#define MSGFIRK_BAD_RELTOL    "reltol < 0 illegal."
#define MSGFIRK_BAD_ABSTOL    "abstol has negative component(s) (illegal)."
#define MSGFIRK_NULL_ABSTOL   "abstol = NULL illegal."
#define MSGFIRK_NULL_Y0       "y0 = NULL illegal."
#define MSGFIRK_NULL_F        "f = NULL illegal."
#define MSGFIRK_NULL_G        "g = NULL illegal."
#define MSGFIRK_BAD_NVECTOR   "A required vector operation is not implemented."
#define MSGFIRK_STACK_LOCAL_OPS                                     \
  "A vector with a communicator must implement nvdotprodlocal and " \
  "nvwsqrsumlocal."
#define MSGFIRK_BAD_K       "Illegal value for k."
#define MSGFIRK_NULL_DKY    "dky = NULL illegal."
#define MSGFIRK_BAD_T       "Illegal value for t. " MSG_TIME_INT
#define MSGFIRK_NO_ROOT     "Rootfinding was not initialized."
#define MSGFIRK_ONLY_STAGES "Only s = 3 stages are supported."
#define MSGFIRK_ONLY_ORDER                                                     \
  "Only the orders 4 and 5 (s = 3 stages) are supported, or ord <= 0 for the " \
  "default."
#define MSGFIRK_BAD_STAGES \
  "Illegal number of stages (must be between 1 and FIRK_MAX_STAGES)."
#define MSGFIRK_STAGES_LOCKED                                          \
  "The number of stages can only be changed before the first step or " \
  "after FIRKodeReInit."
#define MSGFIRK_BAD_MODE "Illegal mode (must be FIRK_BEGINNER or FIRK_EXPERT)."
#define MSGFIRK_MODE_LOCKED \
  "The mode can only be set before the first step (or after FIRKodeReInit)."
#define MSGFIRK_REQUIRES_EXPERT "%s requires FIRKodeSetMode(mem, FIRK_EXPERT)."
#define MSGFIRK_BAD_POLICY \
  "Illegal reuse policy (must be FIRK_POLICY_RADAU5 or FIRK_POLICY_PERIODIC)."
#define MSGFIRK_BAD_THETA_JAC \
  "The Jacobian reuse rate must be less than 1 (<= 0 for the default)."
#define MSGFIRK_BAD_CONTROLLER_TYPE \
  "The controller must be of type SUN_ADAPTCONTROLLER_H."
#define MSGFIRK_BAD_CONTROLLER_NAME                                          \
  "Unknown controller name (must be Soderlind, PID, PI, I, ExpGus, ImpGus, " \
  "H0211, H0321, H211 or H312)."
#define MSGFIRK_TABLE_FAIL \
  "No built-in Radau IIA table for this number of stages."
#define MSGFIRK_NO_LS \
  "FIRKODE requires a linear solver; call FIRKodeSetLinearSolver."

/* Evolve error messages */

#define MSGFIRK_NO_TOL    "No integration tolerances have been specified."
#define MSGFIRK_YOUT_NULL "yout = NULL illegal."
#define MSGFIRK_TRET_NULL "tret = NULL illegal."
#define MSGFIRK_BAD_EWT   "Initial ewt has component(s) equal to zero (illegal)."
#define MSGFIRK_EWT_NOW_BAD \
  "At " MSG_TIME ", a component of ewt has become <= 0."
#define MSGFIRK_BAD_ITASK  "Illegal value for itask."
#define MSGFIRK_BAD_H0     "h0 and tout - t0 inconsistent."
#define MSGFIRK_BAD_HFIXED "The fixed step size and tout - t0 are inconsistent."
#define MSGFIRK_BAD_TOUT                    \
  "Trouble interpolating at " MSG_TIME_TOUT \
  ". tout too far back in direction of integration"
#define MSGFIRK_EWT_FAIL "The user-provided EwtSet function failed."
#define MSGFIRK_EWT_NOW_FAIL \
  "At " MSG_TIME ", the user-provided EwtSet function failed."
#define MSGFIRK_LINIT_FAIL     "The linear solver's init routine failed."
#define MSGFIRK_MASSINIT_FAIL  "The mass matrix solver's init routine failed."
#define MSGFIRK_MASSSETUP_FAIL "The mass matrix solver's setup routine failed."
#define MSGFIRK_HNIL_DONE                                                  \
  "The above warning has been issued mxhnil times and will not be issued " \
  "again for this problem."
#define MSGFIRK_TOO_CLOSE "tout too close to t0 to start integration."
#define MSGFIRK_MAX_STEPS \
  "At " MSG_TIME ", mxstep steps taken before reaching tout."
#define MSGFIRK_TOO_MUCH_ACC "At " MSG_TIME ", too much accuracy requested."
#define MSGFIRK_HNIL                                                       \
  "Internal " MSG_TIME_H " are such that t + h = t on the next step. The " \
  "solver will continue anyway."
#define MSGFIRK_ERR_FAILS \
  "At " MSG_TIME_H ", the error test failed repeatedly or with |h| = hmin."
#define MSGFIRK_CONV_FAILS \
  "At " MSG_TIME_H         \
  ", the corrector convergence test failed repeatedly or with |h| = hmin."
#define MSGFIRK_SETUP_FAILED \
  "At " MSG_TIME ", the setup routine failed in an unrecoverable manner."
#define MSGFIRK_SOLVE_FAILED \
  "At " MSG_TIME ", the solve routine failed in an unrecoverable manner."
#define MSGFIRK_RHSFUNC_FAILED \
  "At " MSG_TIME               \
  ", the right-hand side routine failed in an unrecoverable manner."
#define MSGFIRK_RHSFUNC_UNREC                                                 \
  "At " MSG_TIME ", the right-hand side failed in a recoverable manner, but " \
  "no recovery is possible."
#define MSGFIRK_RHSFUNC_REPTD \
  "At " MSG_TIME " repeated recoverable right-hand side function errors."
#define MSGFIRK_RHSFUNC_FIRST \
  "The right-hand side routine failed at the first call."
#define MSGFIRK_RTFUNC_FAILED                                            \
  "At " MSG_TIME ", the rootfinding routine failed in an unrecoverable " \
  "manner."
#define MSGFIRK_CLOSE_ROOTS "Root found at and very near " MSG_TIME "."
#define MSGFIRK_BAD_TSTOP                                    \
  "The value " MSG_TIME_TSTOP " is behind current " MSG_TIME \
  " in the direction of integration."
#define MSGFIRK_INACTIVE_ROOTS                                         \
  "At the end of the first step, there are still some root functions " \
  "identically 0. This warning will not be issued again."
#define MSGFIRK_MASSSOLVE_FAIL \
  "At " MSG_TIME ", the mass matrix solve failed in an unrecoverable manner."
#define MSGFIRK_MASSMULT_FAIL \
  "At " MSG_TIME ", the mass matrix-vector product failed."
#define MSGFIRK_CONTROLLER_ERR \
  "At " MSG_TIME ", the time step controller failed."

#ifdef __cplusplus
}
#endif

#endif
