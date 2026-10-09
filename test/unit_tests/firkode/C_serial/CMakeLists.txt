# ---------------------------------------------------------------
# Programmer(s): Yifan Hu @ UMBC
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025-2026, Lawrence Livermore National Security,
# University of Maryland Baltimore County, and the SUNDIALS contributors.
# Copyright (c) 2013-2025, Lawrence Livermore National Security
# and Southern Methodist University.
# Copyright (c) 2002-2013, Lawrence Livermore National Security.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# FIRKODE C serial unit tests
# ---------------------------------------------------------------

# List of test tuples of the form "name\;args"
set(FIRKODE_unit_tests
    "firk_test_adapt\;"
    "firk_test_dahlquist\;1"
    "firk_test_dahlquist\;2"
    "firk_test_dahlquist\;3"
    "firk_test_dense_output\;1"
    "firk_test_dense_output\;2"
    "firk_test_dense_output\;3"
    "firk_test_errest\;"
    "firk_test_ewt\;"
    "firk_test_firkls\;"
    "firk_test_getuserdata\;"
    "firk_test_hnil\;"
    "firk_test_mode\;"
    "firk_test_newton\;"
    "firk_test_nlscoef\;"
    "firk_test_options\;"
    "firk_test_recoverable\;"
    "firk_test_reset\;"
    "firk_test_return_codes\;"
    "firk_test_reuse\;"
    "firk_test_rootfind\;"
    "firk_test_setorder\;"
    "firk_test_stack_ops\;"
    "firk_test_tables\;"
    "firk_test_tstop\;"
    "firk_test_userdata\;"
    "firk_test_workspace\;")

# Tests that use the shared mock objects of ../firk_test_mocks.c
set(FIRKODE_mock_tests firk_test_firkls)

# Add the build and install targets for each test
foreach(test_tuple ${FIRKODE_unit_tests})

  # parse the test tuple
  list(GET test_tuple 0 test)
  list(GET test_tuple 1 test_args)

  # check if this test has already been added, only need to add test source
  # files once for testing with different inputs
  if(NOT TARGET ${test})

    if(test IN_LIST FIRKODE_mock_tests)
      firkode_add_unit_test_target(${test} ${test}.c ../firk_test_mocks.c)
    else()
      firkode_add_unit_test_target(${test} ${test}.c)
    endif()

  endif()

  # check if test args are provided and set the test name
  if("${test_args}" STREQUAL "")
    set(test_name ${test})
  else()
    string(REPLACE " " "_" test_name "${test}_${test_args}")
    string(REPLACE " " ";" test_args "${test_args}")
  endif()

  # add test to regression tests
  add_test(NAME ${test_name} COMMAND ${test} ${test_args})

endforeach()

# The log records (H31), at logging level 4 or more: the test writes a log and
# checks its balance, then suntools.logs parses the same file
if(SUNDIALS_LOGGING_LEVEL GREATER_EQUAL 4)
  firkode_add_unit_test_target(firk_test_logging firk_test_logging.c)
  set(_firk_log ${CMAKE_CURRENT_BINARY_DIR}/firk_test_logging.log)
  add_test(NAME firk_test_logging COMMAND firk_test_logging ${_firk_log})
  set_tests_properties(
    firk_test_logging PROPERTIES FIXTURES_SETUP firk_test_logging_log LABELS
                                 "logging")
  # the parser needs Python, which the test runner may not have looked for
  if(NOT Python3_FOUND)
    find_package(Python3 COMPONENTS Interpreter)
  endif()
  if(Python3_FOUND)
    add_test(
      NAME firk_test_logging_suntools
      COMMAND
        ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/firk_test_logging.py
        ${_firk_log} ${SUNDIALS_SOURCE_DIR}/suntools/src)
    set_tests_properties(
      firk_test_logging_suntools
      PROPERTIES FIXTURES_REQUIRED firk_test_logging_log LABELS "logging")
  endif()
endif()

message(STATUS "Added FIRKODE C serial unit tests")
