/* -----------------------------------------------------------------
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
 * -----------------------------------------------------------------*/

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <kinsol/kinsol.h>
#include <sundials/sundials_core.hpp>

#include "../../utilities/dumpstderr.hpp"
#include "kinsol/kinsol_impl.h"
#include "sundials/sundials_context.hpp"

static const std::string errfile{"test_error_handling.err"};

class KINErrConditionTest : public testing::Test
{
protected:
  KINErrConditionTest()
  {
    SUNContext_ClearErrHandlers(sunctx);
    SUNContext_PushErrHandler(sunctx, SUNLogErrHandlerFn, NULL);
    SUNContext_GetLogger(sunctx, &logger);
    kinmem = KINCreate(sunctx);
  }

  ~KINErrConditionTest() { KINFree(&kinmem); }

  void* kinmem;
  SUNLogger logger;
  sundials::Context sunctx;
};

TEST_F(KINErrConditionTest, WarningIsPrinted)
{
  const SUNStackTraceFrame* frames = nullptr;
  int count                       = -1;

  ASSERT_EQ(SUNContext_SetStackTraceEnabled(sunctx, SUNTRUE), SUN_SUCCESS);
  SUNErrCode err = SUNLogger_SetWarningFilename(logger, errfile.c_str());
  ASSERT_EQ(err, SUN_SUCCESS);
  KINMemRec* kin_mem = (KINMemRec*)kinmem;
  KINProcessError(kin_mem, KIN_WARNING, __LINE__, __func__, __FILE__, "test");
  err = SUNLogger_Flush(logger, SUN_LOGLEVEL_WARNING);
  ASSERT_EQ(err, SUN_SUCCESS);
  std::string output = dumpstderr(sunctx, errfile);
#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_WARNING
  EXPECT_THAT(output, testing::AllOf(testing::StartsWith("[WARNING]"),
                                     testing::HasSubstr("[rank 0]"),
                                     testing::HasSubstr("test")));
#else
  EXPECT_EQ(output, "");
#endif
  ASSERT_EQ(SUNContext_GetStackTrace(sunctx, &frames, &count), SUN_SUCCESS);
  EXPECT_NE(frames, nullptr);
  EXPECT_EQ(count, 0);
}

TEST_F(KINErrConditionTest, ErrorIsPrinted)
{
  const SUNStackTraceFrame* frames = nullptr;
  int count                       = 0;

  ASSERT_EQ(SUNContext_SetStackTraceEnabled(sunctx, SUNTRUE), SUN_SUCCESS);
  SUNErrCode err = SUNLogger_SetErrorFilename(logger, errfile.c_str());
  ASSERT_EQ(err, SUN_SUCCESS);
  // -1 is an illegal value
  int ierr = KINSetNumMaxIters(kinmem, -1);
  ASSERT_NE(ierr, KIN_SUCCESS);
  std::string output = dumpstderr(sunctx, errfile);
#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_ERROR
  EXPECT_THAT(output, testing::AllOf(testing::StartsWith("[ERROR]"),
                                     testing::HasSubstr("[rank 0]"),
                                     testing::HasSubstr(MSG_BAD_MXITER)));
#else
  EXPECT_EQ(output, "");
#endif
  ASSERT_EQ(SUNContext_GetStackTrace(sunctx, &frames, &count), SUN_SUCCESS);
  ASSERT_NE(frames, nullptr);
  ASSERT_EQ(count, 1);
  EXPECT_STREQ(frames[0].func, "KINSetNumMaxIters");
  EXPECT_EQ(frames[0].code, KIN_ILL_INPUT);
  ASSERT_NE(frames[0].msg, nullptr);
  EXPECT_THAT(frames[0].msg, testing::HasSubstr(MSG_BAD_MXITER));
}
