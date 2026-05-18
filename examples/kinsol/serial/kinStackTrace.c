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
 * feature and print a logical SUNDIALS stack trace after a KINSOL error.
 *
 * By default the example prints the full stack trace, including source file
 * and line locations. Run with --no-print-trace to print a stable summary.
 * -----------------------------------------------------------------*/

#include <kinsol/kinsol.h>
#include <stdio.h>
#include <string.h>
#include <sundials/sundials_context.h>

static void IgnoreError(int line, const char* func, const char* file,
                        const char* msg, SUNErrCode err_code,
                        void* err_user_data, SUNContext sunctx);
static int check_sunerr(SUNErrCode err, const char* funcname);

int main(int argc, char* argv[])
{
  SUNContext sunctx                = NULL;
  void* kinmem                     = NULL;
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

  retval = SUNContext_PopErrHandler(sunctx);
  if (check_sunerr(retval, "SUNContext_PopErrHandler")) { goto cleanup; }

  retval = SUNContext_PushErrHandler(sunctx, IgnoreError, NULL);
  if (check_sunerr(retval, "SUNContext_PushErrHandler")) { goto cleanup; }

  kinmem = KINCreate(sunctx);
  if (kinmem == NULL)
  {
    fprintf(stderr, "KINCreate failed\n");
    SUNContext_PrintStackTrace(sunctx, stderr);
    goto cleanup;
  }

  /* This is intentionally invalid before KINInit. */
  retval = KINSol(kinmem, NULL, KIN_NONE, NULL, NULL);
  if (retval != KIN_NO_MALLOC)
  {
    fprintf(stderr, "Expected KIN_NO_MALLOC, but got %d\n", retval);
    goto cleanup;
  }

  retval = SUNContext_GetStackTrace(sunctx, &frames, &nframes);
  if (check_sunerr(retval, "SUNContext_GetStackTrace")) { goto cleanup; }

  printf("KINSol returned KIN_NO_MALLOC as expected.\n");
  printf("Recorded %d stack trace frame(s).\n", nframes);

  if (nframes > 0)
  {
    printf("Leaf frame: %s\n", frames[0].func);
    printf("Error code: %s\n",
           (frames[0].code == KIN_NO_MALLOC) ? "KIN_NO_MALLOC" : "unexpected");
    printf("Message: %s\n", frames[0].msg);
  }

  if (print_trace) { SUNContext_PrintStackTrace(sunctx, stdout); }

  status = 0;

cleanup:
  KINFree(&kinmem);
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

static int check_sunerr(SUNErrCode err, const char* funcname)
{
  if (err == SUN_SUCCESS) { return 0; }

  fprintf(stderr, "%s failed with error code %d\n", funcname, err);
  return 1;
}
