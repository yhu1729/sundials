/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
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
 * SUNDIALS context class. A context object holds data that all
 * SUNDIALS objects in a simulation share.
 * ----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sundials/priv/sundials_context_impl.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_errors.h>
#include <sundials/sundials_logger.h>
#include <sundials/sundials_profiler.h>
#include <sundials/sundials_types.h>

#include "sundials_adiak_metadata.h"
#include "sundials_macros.h"

#define SUNCTX_DEFAULT_STACK_TRACE_DEPTH 32

/* Forward declaration of function used to destroy any data allocated for Python */
#if defined(SUNDIALS_ENABLE_PYTHON)
void SUNContextFunctionTable_Destroy(void* ptr);
#endif

static char* sunContext_CopyStackTraceMessage(const char* msg, SUNErrCode code)
{
  const char* src = (msg != NULL) ? msg : SUNGetErrMsg(code);
  size_t len;
  char* copy;

  if (src == NULL) { return NULL; }

  len  = strlen(src);
  copy = (char*)malloc(len + 1);
  if (copy == NULL) { return NULL; }

  memcpy(copy, src, len + 1);
  return copy;
}

static void sunContext_FreeStackTraceMessages(SUNContext sunctx)
{
  int i;

  if (sunctx == NULL || sunctx->stack_trace == NULL) { return; }

  for (i = 0; i < sunctx->stack_trace_count; i++)
  {
    free((void*)sunctx->stack_trace[i].msg);
    sunctx->stack_trace[i].msg = NULL;
  }
}

static void sunContext_ResetStackTrace(SUNContext sunctx)
{
  if (sunctx == NULL) { return; }

  sunContext_FreeStackTraceMessages(sunctx);
  sunctx->stack_trace_count     = 0;
  sunctx->stack_trace_err       = SUN_SUCCESS;
  sunctx->stack_trace_truncated = SUNFALSE;
}

static void sunContext_AppendStackTraceFrame(SUNContext sunctx, int line,
                                             const char* func, const char* file,
                                             const char* msg, SUNErrCode code)
{
  SUNStackTraceFrame* frame;

  if (sunctx == NULL || sunctx->stack_trace == NULL ||
      sunctx->stack_trace_count >= sunctx->stack_trace_capacity)
  {
    if (sunctx != NULL) { sunctx->stack_trace_truncated = SUNTRUE; }
    return;
  }

  frame       = &sunctx->stack_trace[sunctx->stack_trace_count];
  frame->func = func;
  frame->file = file;
  frame->msg  = sunContext_CopyStackTraceMessage(msg, code);
  frame->line = line;
  frame->code = code;

  sunctx->stack_trace_count++;
}

void sunContext_TraceRaise(SUNContext sunctx, int line, const char* func,
                           const char* file, const char* msg, SUNErrCode code)
{
  if (sunctx == NULL || !sunctx->stack_trace_enabled || code >= 0) { return; }

  sunContext_ResetStackTrace(sunctx);
  sunctx->stack_trace_err = code;
  sunContext_AppendStackTraceFrame(sunctx, line, func, file, msg, code);
}

void sunContext_TracePropagate(SUNContext sunctx, int line, const char* func,
                               const char* file, const char* msg,
                               SUNErrCode code)
{
  if (sunctx == NULL || !sunctx->stack_trace_enabled || code >= 0) { return; }

  if (sunctx->stack_trace_count == 0 || sunctx->stack_trace_err != code)
  {
    sunContext_ResetStackTrace(sunctx);
    sunctx->stack_trace_err = code;
  }

  sunContext_AppendStackTraceFrame(sunctx, line, func, file, msg, code);
}

SUNErrCode SUNContext_Create(SUNComm comm, SUNContext* sunctx_out)
{
  SUNErrCode err       = SUN_SUCCESS;
  SUNProfiler profiler = NULL;
  SUNLogger logger     = NULL;
  SUNContext sunctx    = NULL;
  SUNErrHandler eh     = NULL;

  *sunctx_out = NULL;
  sunctx      = (SUNContext)malloc(sizeof(struct SUNContext_));

  /* SUNContext_Create cannot assert or log since the SUNContext is not yet
   * created */
  if (!sunctx) { return SUN_ERR_MALLOC_FAIL; }

  sunctx->python       = NULL;
  sunctx->logger       = NULL;
  sunctx->own_logger   = SUNFALSE;
  sunctx->profiler     = NULL;
  sunctx->own_profiler = SUNFALSE;
  sunctx->last_err     = SUN_SUCCESS;
  sunctx->err_handler  = NULL;
  sunctx->comm         = comm;
  sunctx->stack_trace_enabled   = SUNFALSE;
  sunctx->stack_trace           = NULL;
  sunctx->stack_trace_count     = 0;
  sunctx->stack_trace_capacity  = SUNCTX_DEFAULT_STACK_TRACE_DEPTH;
  sunctx->stack_trace_err       = SUN_SUCCESS;
  sunctx->stack_trace_truncated = SUNFALSE;

  SUNFunctionBegin(sunctx);

#ifdef SUNDIALS_ADIAK_ENABLED
  adiak_init(&comm);
  sunAdiakCollectMetadata();
#endif

  do {
#if SUNDIALS_LOGGING_LEVEL > 0
#if SUNDIALS_MPI_ENABLED
    err = SUNLogger_CreateFromEnv(comm, &logger);
    SUNCheckCallNoRet(err);
    if (err) { break; }
#else
    err = SUNLogger_CreateFromEnv(SUN_COMM_NULL, &logger);
    SUNCheckCallNoRet(err);
    if (err) { break; }
#endif
#else
    err = SUNLogger_Create(SUN_COMM_NULL, 0, &logger);
    SUNCheckCallNoRet(err);
    if (err) { break; }
    err = SUNLogger_SetErrorFilename(logger, "");
    SUNCheckCallNoRet(err);
    if (err) { break; }
    err = SUNLogger_SetWarningFilename(logger, "");
    SUNCheckCallNoRet(err);
    if (err) { break; }
    err = SUNLogger_SetInfoFilename(logger, "");
    SUNCheckCallNoRet(err);
    if (err) { break; }
    err = SUNLogger_SetDebugFilename(logger, "");
    SUNCheckCallNoRet(err);
    if (err) { break; }
#endif

#if defined(SUNDIALS_ENABLE_PROFILING) && !defined(SUNDIALS_CALIPER_ENABLED)
    err = SUNProfiler_Create(comm, "SUNContext Default", &profiler);
    SUNCheckCallNoRet(err);
    if (err) { break; }
#endif

    err = SUNErrHandler_Create(SUNLogErrHandlerFn, NULL, &eh);
    SUNCheckCallNoRet(err);
    if (err) { break; }

    sunctx->python       = NULL;
    sunctx->logger       = logger;
    sunctx->own_logger   = logger != NULL;
    sunctx->profiler     = profiler;
    sunctx->own_profiler = profiler != NULL;
    sunctx->last_err     = SUN_SUCCESS;
    sunctx->err_handler  = eh;
    sunctx->comm         = comm;
    sunctx->stack_trace_enabled   = SUNFALSE;
    sunctx->stack_trace           = NULL;
    sunctx->stack_trace_count     = 0;
    sunctx->stack_trace_capacity  = SUNCTX_DEFAULT_STACK_TRACE_DEPTH;
    sunctx->stack_trace_err       = SUN_SUCCESS;
    sunctx->stack_trace_truncated = SUNFALSE;
  }
  while (0);

  if (err)
  {
#if defined(SUNDIALS_ENABLE_PROFILING) && !defined(SUNDIALS_CALIPER_ENABLED)
    SUNCheckCallNoRet(SUNProfiler_Free(&profiler));
#endif
    SUNCheckCallNoRet(SUNLogger_Destroy(&logger));
    free(sunctx);
  }
  else { *sunctx_out = sunctx; }

  return err;
}

SUNErrCode SUNContext_GetLastError(SUNContext sunctx)
{
  if (!sunctx) { return SUN_ERR_SUNCTX_CORRUPT; }

  SUNFunctionBegin(sunctx);
  SUNErrCode err   = sunctx->last_err;
  sunctx->last_err = SUN_SUCCESS;
  return err;
}

SUNErrCode SUNContext_PeekLastError(SUNContext sunctx)
{
  if (!sunctx) { return SUN_ERR_SUNCTX_CORRUPT; }

  SUNFunctionBegin(sunctx);
  return sunctx->last_err;
}

SUNErrCode SUNContext_PushErrHandler(SUNContext sunctx, SUNErrHandlerFn err_fn,
                                     void* err_user_data)
{
  if (!sunctx || !err_fn) { return SUN_ERR_SUNCTX_CORRUPT; }

  SUNFunctionBegin(sunctx);
  SUNErrHandler new_err_handler = NULL;
  if (SUNErrHandler_Create(err_fn, err_user_data, &new_err_handler))
  {
    return SUN_ERR_CORRUPT;
  }
  new_err_handler->previous = sunctx->err_handler;
  sunctx->err_handler       = new_err_handler;
  return SUN_SUCCESS;
}

SUNErrCode SUNContext_PopErrHandler(SUNContext sunctx)
{
  if (!sunctx) { return SUN_ERR_SUNCTX_CORRUPT; }

  SUNFunctionBegin(sunctx);
  if (sunctx->err_handler)
  {
    SUNErrHandler eh = sunctx->err_handler;
    if (sunctx->err_handler->previous)
    {
      sunctx->err_handler = sunctx->err_handler->previous;
    }
    else { sunctx->err_handler = NULL; }
    SUNErrHandler_Destroy(&eh);
  }
  return SUN_SUCCESS;
}

SUNErrCode SUNContext_ClearErrHandlers(SUNContext sunctx)
{
  if (!sunctx) { return SUN_ERR_SUNCTX_CORRUPT; }

  SUNFunctionBegin(sunctx);
  while (sunctx->err_handler != NULL)
  {
    SUNCheckCall(SUNContext_PopErrHandler(sunctx));
  }
  return SUN_SUCCESS;
}

SUNErrCode SUNContext_SetStackTraceEnabled(SUNContext sunctx,
                                           sunbooleantype enabled)
{
  if (!sunctx) { return SUN_ERR_SUNCTX_CORRUPT; }

  SUNFunctionBegin(sunctx);

  if (enabled)
  {
    if (sunctx->stack_trace_enabled) { return SUN_SUCCESS; }

    sunctx->stack_trace = (SUNStackTraceFrame*)calloc(
      (size_t)sunctx->stack_trace_capacity, sizeof(SUNStackTraceFrame));
    if (sunctx->stack_trace == NULL) { return SUN_ERR_MALLOC_FAIL; }

    sunctx->stack_trace_count     = 0;
    sunctx->stack_trace_err       = SUN_SUCCESS;
    sunctx->stack_trace_truncated = SUNFALSE;
    sunctx->stack_trace_enabled   = SUNTRUE;
    return SUN_SUCCESS;
  }

  sunContext_ResetStackTrace(sunctx);
  free(sunctx->stack_trace);
  sunctx->stack_trace         = NULL;
  sunctx->stack_trace_enabled = SUNFALSE;
  return SUN_SUCCESS;
}

SUNErrCode SUNContext_SetStackTraceMaxDepth(SUNContext sunctx, int max_depth)
{
  SUNStackTraceFrame* new_stack_trace = NULL;
  sunbooleantype old_truncated;
  int new_count;
  int old_count;
  int i;

  if (!sunctx) { return SUN_ERR_SUNCTX_CORRUPT; }
  if (max_depth <= 0) { return SUN_ERR_ARG_OUTOFRANGE; }

  SUNFunctionBegin(sunctx);

  if (!sunctx->stack_trace_enabled)
  {
    sunctx->stack_trace_capacity = max_depth;
    return SUN_SUCCESS;
  }

  new_stack_trace = (SUNStackTraceFrame*)calloc((size_t)max_depth,
                                                sizeof(SUNStackTraceFrame));
  if (new_stack_trace == NULL) { return SUN_ERR_MALLOC_FAIL; }

  old_count = sunctx->stack_trace_count;
  old_truncated = sunctx->stack_trace_truncated;
  new_count = (old_count < max_depth) ? old_count : max_depth;
  for (i = 0; i < new_count; i++)
  {
    new_stack_trace[i] = sunctx->stack_trace[i];
  }

  for (i = new_count; i < old_count; i++)
  {
    free((void*)sunctx->stack_trace[i].msg);
  }

  free(sunctx->stack_trace);
  sunctx->stack_trace           = new_stack_trace;
  sunctx->stack_trace_count     = new_count;
  sunctx->stack_trace_capacity  = max_depth;
  sunctx->stack_trace_truncated = old_truncated || (new_count < old_count);

  return SUN_SUCCESS;
}

SUNErrCode SUNContext_GetStackTrace(SUNContext sunctx,
                                    const SUNStackTraceFrame** frames,
                                    int* count)
{
  if (!sunctx) { return SUN_ERR_SUNCTX_CORRUPT; }
  if (frames == NULL || count == NULL) { return SUN_ERR_ARG_CORRUPT; }

  SUNFunctionBegin(sunctx);

  *frames = sunctx->stack_trace;
  *count  = sunctx->stack_trace_count;
  return SUN_SUCCESS;
}

SUNErrCode SUNContext_ClearStackTrace(SUNContext sunctx)
{
  if (!sunctx) { return SUN_ERR_SUNCTX_CORRUPT; }

  SUNFunctionBegin(sunctx);
  sunContext_ResetStackTrace(sunctx);
  return SUN_SUCCESS;
}

SUNErrCode SUNContext_PrintStackTrace(SUNContext sunctx, FILE* fp)
{
  int i;

  if (!sunctx) { return SUN_ERR_SUNCTX_CORRUPT; }
  if (fp == NULL) { return SUN_ERR_ARG_CORRUPT; }

  SUNFunctionBegin(sunctx);

  fprintf(fp, "SUNDIALS stack trace (%d frame%s):\n",
          sunctx->stack_trace_count,
          (sunctx->stack_trace_count == 1) ? "" : "s");
  for (i = 0; i < sunctx->stack_trace_count; i++)
  {
    SUNStackTraceFrame* frame = &sunctx->stack_trace[i];
    fprintf(fp, "  [%d] %s at %s:%d: %s\n", i,
            (frame->func != NULL) ? frame->func : "<unknown>",
            (frame->file != NULL) ? frame->file : "<unknown>", frame->line,
            (frame->msg != NULL) ? frame->msg : SUNGetErrMsg(frame->code));
  }
  if (sunctx->stack_trace_truncated)
  {
    fprintf(fp, "  ... stack trace truncated at %d frames\n",
            sunctx->stack_trace_capacity);
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNContext_GetProfiler(SUNContext sunctx, SUNProfiler* profiler)
{
  if (!sunctx) { return SUN_ERR_SUNCTX_CORRUPT; }

  SUNFunctionBegin(sunctx);

#ifdef SUNDIALS_ENABLE_PROFILING
  /* get profiler */
  *profiler = sunctx->profiler;
#else
  *profiler = NULL;
#endif

  return SUN_SUCCESS;
}

SUNErrCode SUNContext_SetProfiler(SUNContext sunctx, SUNProfiler profiler)
{
  if (!sunctx) { return SUN_ERR_SUNCTX_CORRUPT; }

  SUNFunctionBegin(sunctx);

#ifdef SUNDIALS_ENABLE_PROFILING
  /* free any existing profiler */
  if (sunctx->profiler && sunctx->own_profiler)
  {
    SUNCheckCall(SUNProfiler_Free(&(sunctx->profiler)));
    sunctx->profiler = NULL;
  }

  /* set profiler */
  sunctx->profiler     = profiler;
  sunctx->own_profiler = SUNFALSE;
#else
  /* silence warnings when profiling is disabled */
  ((void)profiler);
#endif

  return SUN_SUCCESS;
}

SUNErrCode SUNContext_GetLogger(SUNContext sunctx, SUNLogger* logger)
{
  if (!sunctx) { return SUN_ERR_SUNCTX_CORRUPT; }

  SUNFunctionBegin(sunctx);

  /* get logger */
  *logger = sunctx->logger;
  return SUN_SUCCESS;
}

SUNErrCode SUNContext_SetLogger(SUNContext sunctx, SUNLogger logger)
{
  if (!sunctx) { return SUN_ERR_SUNCTX_CORRUPT; }

  SUNFunctionBegin(sunctx);

  /* free any existing logger */
  if (sunctx->logger && sunctx->own_logger)
  {
    if (SUNLogger_Destroy(&(sunctx->logger))) { return SUN_ERR_DESTROY_FAIL; }
    sunctx->logger = NULL;
  }

  /* set logger */
  sunctx->logger     = logger;
  sunctx->own_logger = SUNFALSE;

  return SUN_SUCCESS;
}

SUNErrCode SUNContext_Free(SUNContext* sunctx)
{
#ifdef SUNDIALS_ADIAK_ENABLED
  adiak_fini();
#endif

  if (!sunctx || !(*sunctx)) { return SUN_SUCCESS; }

#if defined(SUNDIALS_ENABLE_PROFILING) && !defined(SUNDIALS_CALIPER_ENABLED)
  /* Find out where we are printing to */
  FILE* fp                    = NULL;
  char* sunprofiler_print_env = getenv("SUNPROFILER_PRINT");
  fp                          = NULL;
  if (sunprofiler_print_env)
  {
    if (!strcmp(sunprofiler_print_env, "0")) { fp = NULL; }
    else if (!strcmp(sunprofiler_print_env, "1") ||
             !strcmp(sunprofiler_print_env, "TRUE") ||
             !strcmp(sunprofiler_print_env, "stdout"))
    {
      fp = stdout;
    }
    else { fp = fopen(sunprofiler_print_env, "a"); }
  }

  /* Enforce that the profiler is freed before finalizing,
     if it is not owned by the sunctx. */
  if ((*sunctx)->profiler)
  {
    if (fp) { SUNProfiler_Print((*sunctx)->profiler, fp); }
    if (fp) { fclose(fp); }
    if ((*sunctx)->own_profiler) { SUNProfiler_Free(&(*sunctx)->profiler); }
  }
#endif

  if ((*sunctx)->logger && (*sunctx)->own_logger)
  {
    SUNLogger_Destroy(&(*sunctx)->logger);
  }

  SUNContext_ClearErrHandlers(*sunctx);

  sunContext_ResetStackTrace(*sunctx);
  free((*sunctx)->stack_trace);
  (*sunctx)->stack_trace = NULL;

#if defined(SUNDIALS_ENABLE_PYTHON)
  SUNContextFunctionTable_Destroy((*sunctx)->python);
#endif
  (*sunctx)->python = NULL;

  free(*sunctx);
  *sunctx = NULL;

  return SUN_SUCCESS;
}
