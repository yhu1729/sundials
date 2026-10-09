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
# FIRKODE C++ serial unit tests
# ---------------------------------------------------------------

# List of test tuples of the form "name\;args\;type"
set(unit_tests
    "firk_test_order.cpp\;pr 1 -1\;exclude-single"
    "firk_test_order.cpp\;pr 2 -1\;exclude-single"
    "firk_test_order.cpp\;pr 3 -1\;exclude-single"
    "firk_test_order.cpp\;pr 1 -1e4\;exclude-single"
    "firk_test_order.cpp\;pr 2 -1e4\;exclude-single"
    "firk_test_order.cpp\;pr 3 -1e4\;exclude-single"
    "firk_test_order.cpp\;kpr 1\;exclude-single"
    "firk_test_order.cpp\;kpr 2\;exclude-single"
    "firk_test_order.cpp\;kpr 3\;exclude-single"
    "firk_test_krylov.cpp\;\;exclude-single"
    "firk_test_mass.cpp\;\;exclude-single"
    "firk_test_mass_fem.cpp\;\;exclude-single"
    "firk_test_mass_fem_krylov.cpp\;\;exclude-single"
    "firk_test_mass_fem_krylov.cpp\;3\;exclude-single"
    "firk_test_mass_fem_krylov.cpp\;3 step\;exclude-single"
    "firk_test_robertson.cpp\;\;exclude-single"
    "firk_test_root.cpp\;\;exclude-single"
    "firk_test_vdp.cpp\;\;exclude-single")

# Tests of Krylov paths, whose iteration counts differ across platforms, are
# judged by their exit code only (design Table 11.1)
set(nodiff_tests firk_test_krylov firk_test_mass firk_test_mass_fem_krylov)

# Add the build and install targets for each test
foreach(test_tuple ${unit_tests})

  # parse the test tuple
  list(GET test_tuple 0 test)
  list(GET test_tuple 1 test_args)
  list(GET test_tuple 2 test_type)

  # extract the file name without extension
  get_filename_component(test_target ${test} NAME_WE)

  # check if this test has already been added, only need to add test source
  # files once for testing with different inputs
  if(NOT TARGET ${test_target})

    firkode_add_unit_test_target(${test_target} ${test})

  endif()

  # check if test args are provided and set the test name
  if("${test_args}" STREQUAL "")
    set(test_name ${test_target})
  else()
    string(REGEX REPLACE " " "_" test_name ${test_target}_${test_args})
  endif()

  # add test to regression tests
  if(test_target IN_LIST nodiff_tests)
    sundials_add_test(
      ${test_name} ${test_target}
      TEST_ARGS ${test_args}
      EXAMPLE_TYPE ${test_type}
      NODIFF)
  else()
    sundials_add_test(
      ${test_name} ${test_target}
      TEST_ARGS ${test_args}
      ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
      ANSWER_FILE ${test_name}.out
      EXAMPLE_TYPE ${test_type})
  endif()

endforeach()

# The comparison with radau5.f (design Section 10.3) reads the reference files
# of ../radau5; its accept/reject decisions flip at roundoff across platforms,
# so it is judged by its exit code
firkode_add_unit_test_target(firk_test_radau5_steps firk_test_radau5_steps.cpp)
sundials_add_test(
  firk_test_radau5_steps firk_test_radau5_steps
  TEST_ARGS ${CMAKE_CURRENT_SOURCE_DIR}/../radau5
  EXAMPLE_TYPE exclude-single
  NODIFF)

message(STATUS "Added FIRKODE CXX serial unit tests")
