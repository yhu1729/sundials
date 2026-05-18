/* -----------------------------------------------------------------
 * Programmer(s): SUNDIALS team
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
 * This example shows how to enable the optional SUNContext stack trace
 * feature and inspect the logical SUNDIALS call path after an error.
 *
 * By default the example prints the full stack trace, including source file
 * and line locations. Run with --no-print-trace to print a stable summary.
 * -----------------------------------------------------------------*/

#include <cvode/cvode.h>
#include <nvector/nvector_manyvector.h>
#include <nvector/nvector_serial.h>
#include <stdio.h>
#include <string.h>
#include <sundials/sundials_context.h>

static void IgnoreError(int line, const char* func, const char* file,
                        const char* msg, SUNErrCode err_code,
                        void* err_user_data, SUNContext sunctx);
static int RunDeepStackTraceExample(SUNContext sunctx,
                                    sunbooleantype print_trace);
static int check_sunerr(SUNErrCode err, const char* funcname);

int main(int argc, char* argv[])
{
  SUNContext sunctx                = NULL;
  void* cvode_mem                  = NULL;
  const SUNStackTraceFrame* frames = NULL;
  int nframes                      = 0;
  int retval                       = 0;
  int status                       = 1;
  sunbooleantype print_trace       = SUNTRUE;

  if ((argc > 1) && (strcmp(argv[1], "--no-print-trace") == 0))
  {
    print_trace = SUNFALSE;
  }

  retval = SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (check_sunerr(retval, "SUNContext_Create")) { goto cleanup; }

  retval = SUNContext_SetStackTraceEnabled(sunctx, SUNTRUE);
  if (check_sunerr(retval, "SUNContext_SetStackTraceEnabled")) { goto cleanup; }

  /*
   * Keep this example's regression output stable. Applications that want the
   * default SUNDIALS error log should leave the default handler installed.
   */
  retval = SUNContext_PopErrHandler(sunctx);
  if (check_sunerr(retval, "SUNContext_PopErrHandler")) { goto cleanup; }

  retval = SUNContext_PushErrHandler(sunctx, IgnoreError, NULL);
  if (check_sunerr(retval, "SUNContext_PushErrHandler")) { goto cleanup; }

  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (cvode_mem == NULL)
  {
    fprintf(stderr, "CVodeCreate failed\n");
    SUNContext_PrintStackTrace(sunctx, stderr);
    goto cleanup;
  }

  /* This is intentionally invalid before CVodeInit. */
  retval = CVodeSStolerances(cvode_mem, 1.0e-4, 1.0e-8);
  if (retval != CV_NO_MALLOC)
  {
    fprintf(stderr, "Expected CV_NO_MALLOC, but got %d\n", retval);
    goto cleanup;
  }

  retval = SUNContext_GetStackTrace(sunctx, &frames, &nframes);
  if (check_sunerr(retval, "SUNContext_GetStackTrace")) { goto cleanup; }

  printf("CVodeSStolerances returned CV_NO_MALLOC as expected.\n");
  printf("Recorded %d stack trace frame(s).\n", nframes);

  if (nframes > 0)
  {
    printf("Leaf frame: %s\n", frames[0].func);
    printf("Error code: %s\n",
           (frames[0].code == CV_NO_MALLOC) ? "CV_NO_MALLOC" : "unexpected");
    printf("Message: %s\n", frames[0].msg);
  }

  if (print_trace) { SUNContext_PrintStackTrace(sunctx, stdout); }

  if (RunDeepStackTraceExample(sunctx, print_trace)) { goto cleanup; }

  retval = SUNContext_ClearStackTrace(sunctx);
  if (check_sunerr(retval, "SUNContext_ClearStackTrace")) { goto cleanup; }

  status = 0;

cleanup:
  CVodeFree(&cvode_mem);
  SUNContext_Free(&sunctx);
  return status;
}

static void IgnoreError(int line, const char* func, const char* file,
                        const char* msg, SUNErrCode err_code,
                        void* err_user_data, SUNContext sunctx)
{
  (void)line;
  (void)func;
  (void)file;
  (void)msg;
  (void)err_code;
  (void)err_user_data;
  (void)sunctx;
}

static int RunDeepStackTraceExample(SUNContext sunctx,
                                    sunbooleantype print_trace)
{
  N_Vector subvec                  = NULL;
  N_Vector manyvec                 = NULL;
  N_Vector subvectors[1]           = {NULL};
  const SUNStackTraceFrame* frames = NULL;
  sunindextype buffer_size         = 0;
  int nframes                      = 0;
  SUNErrCode retval                = SUN_SUCCESS;
  int status                       = 1;

  retval = SUNContext_ClearStackTrace(sunctx);
  if (check_sunerr(retval, "SUNContext_ClearStackTrace")) { goto cleanup; }

  subvec = N_VNew_Serial(1, sunctx);
  if (subvec == NULL)
  {
    fprintf(stderr, "N_VNew_Serial failed\n");
    goto cleanup;
  }

  subvectors[0] = subvec;
  manyvec       = N_VNew_ManyVector(1, subvectors, sunctx);
  if (manyvec == NULL)
  {
    fprintf(stderr, "N_VNew_ManyVector failed\n");
    goto cleanup;
  }

  /*
   * Force a missing optional operation so this example records a deep logical
   * trace through N_VBufSize_ManyVector into the serial subvector.
   */
  subvec->ops->nvbufsize = NULL;

  retval = N_VBufSize(manyvec, &buffer_size);
  if (retval != SUN_ERR_NOT_IMPLEMENTED)
  {
    fprintf(stderr, "Expected SUN_ERR_NOT_IMPLEMENTED, but got %d\n", retval);
    goto cleanup;
  }

  retval = SUNContext_GetStackTrace(sunctx, &frames, &nframes);
  if (check_sunerr(retval, "SUNContext_GetStackTrace")) { goto cleanup; }

  printf("N_VBufSize returned SUN_ERR_NOT_IMPLEMENTED as expected.\n");
  printf("Deep stack trace recorded %d frame(s).\n", nframes);

  if (nframes != 2)
  {
    fprintf(stderr, "Expected 2 deep stack trace frames.\n");
    goto cleanup;
  }

  if ((frames[0].code != SUN_ERR_NOT_IMPLEMENTED) ||
      (frames[1].code != SUN_ERR_NOT_IMPLEMENTED) ||
      (strcmp(frames[0].func, "N_VBufSize") != 0) ||
      (strcmp(frames[1].func, "N_VBufSize_ManyVector") != 0))
  {
    fprintf(stderr, "Unexpected deep stack trace frame order.\n");
    goto cleanup;
  }

  printf("Deep leaf frame: %s\n", frames[0].func);
  printf("Deep propagation frame: %s\n", frames[1].func);

  if (print_trace) { SUNContext_PrintStackTrace(sunctx, stdout); }

  status = 0;

cleanup:
  N_VDestroy(manyvec);
  N_VDestroy(subvec);
  return status;
}

static int check_sunerr(SUNErrCode err, const char* funcname)
{
  if (err == SUN_SUCCESS) { return 0; }

  fprintf(stderr, "%s failed with error code %d\n", funcname, err);
  return 1;
}
