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
# FIRKODE C parallel unit tests
# ---------------------------------------------------------------

# List of test tuples of the form "name\;tasks"
set(FIRKODE_parallel_tests "firk_test_stack_mpi\;1" "firk_test_stack_mpi\;4")

foreach(test_tuple ${FIRKODE_parallel_tests})

  list(GET test_tuple 0 test)
  list(GET test_tuple 1 number_of_tasks)

  if(NOT TARGET ${test})
    firkode_add_unit_test_target(${test} ${test}.c)
    target_link_libraries(${test} PRIVATE MPI::MPI_C sundials_nvecparallel)
  endif()

  # the tests check themselves; the output differs with the number of tasks
  sundials_add_test(
    ${test}_${number_of_tasks} ${test}
    MPI_NPROCS ${number_of_tasks}
    NODIFF)

endforeach()

message(STATUS "Added FIRKODE C parallel unit tests")
