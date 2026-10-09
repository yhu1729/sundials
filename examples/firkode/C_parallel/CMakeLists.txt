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
# CMakeLists.txt file for FIRKODE C parallel examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;tasks\;type". The examples check themselves
# against a serial run on rank 0 and are judged by their exit code (their Krylov
# iteration counts differ across platforms, design Table 11.1).
set(FIRKODE_examples "firk_robertson_batched_p\;1\;exclude-single"
                     "firk_robertson_batched_p\;4\;exclude-single")

foreach(example_tuple ${FIRKODE_examples})

  list(GET example_tuple 0 example)
  list(GET example_tuple 1 number_of_tasks)
  list(GET example_tuple 2 example_type)

  if(NOT TARGET ${example})
    sundials_add_executable(${example} ${example}.c)
    set_target_properties(${example} PROPERTIES FOLDER "Examples")
    target_link_libraries(${example} sundials_firkode sundials_nvecparallel
                          sundials_nvecserial MPI::MPI_C ${EXE_EXTRA_LINK_LIBS})
  endif()

  sundials_add_test(
    ${example}_${number_of_tasks} ${example}
    MPI_NPROCS ${number_of_tasks}
    EXAMPLE_TYPE ${example_type}
    NODIFF)

  # install example source
  if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
    install(FILES ${example}.c
            DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/firkode/C_parallel)
  endif()

endforeach()

# create Makefile and CMakeLists.txt for examples
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)

  # Install the README file
  install(FILES README
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/firkode/C_parallel)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER "FIRKODE")
  set(SOLVER_LIB "sundials_firkode")

  examples2string(FIRKODE_examples EXAMPLES)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then be used
  # as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_parallel_C_ex.in
    ${PROJECT_BINARY_DIR}/examples/firkode/C_parallel/CMakeLists.txt @ONLY)

  # install CMakelists.txt
  install(FILES ${PROJECT_BINARY_DIR}/examples/firkode/C_parallel/CMakeLists.txt
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/firkode/C_parallel)

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template for the
  # user's own programs.

  if(UNIX)
    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_parallel_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/firkode/C_parallel/Makefile_ex @ONLY)
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/firkode/C_parallel/Makefile_ex
      DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/firkode/C_parallel
      RENAME Makefile)
  endif()

  # add test_install target
  sundials_add_test_install(firkode C_parallel)

endif()
