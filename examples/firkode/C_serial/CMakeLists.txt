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
# CMakeLists.txt file for FIRKODE C serial examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is 'develop' for
# examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(FIRKODE_examples
    # tests that are always run
    "firk_analytic.c\;\;"
    # develop tests
    "firk_analytic.c\;1\;develop"
    "firk_analytic.c\;2\;develop"
    "firk_analytic.c\;5\;develop"
    "firk_diurnal_kry.c\;\;develop"
    "firk_heat1D_fem_mass.c\;\;develop"
    "firk_robertson.c\;\;exclude-single"
    "firk_robertson_root.c\;\;exclude-single")

# Auxiliary files to install
set(FIRKODE_extras)

# Add the build and install targets for each FIRKODE example
foreach(example_tuple ${FIRKODE_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 example_type)

  # extract the file name without extension
  get_filename_component(example_target ${example} NAME_WE)

  if(NOT TARGET ${example_target})
    # example source files
    sundials_add_executable(${example_target} ${example})

    # folder for IDEs
    set_target_properties(${example_target} PROPERTIES FOLDER "Examples")

    # libraries to link against
    target_link_libraries(${example_target} sundials_firkode
                          sundials_nvecserial ${EXE_EXTRA_LINK_LIBS})
  endif()

  # check if example args are provided and set the test name
  if("${example_args}" STREQUAL "")
    set(test_name ${example_target})
  else()
    string(REGEX REPLACE " " "_" test_name ${example_target}_${example_args})
  endif()

  # add example to regression tests
  sundials_add_test(
    ${test_name} ${example_target}
    TEST_ARGS ${example_args}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${test_name}.out
    EXAMPLE_TYPE ${example_type})

  # install example source and out files
  if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
    install(FILES ${example} ${test_name}.out
            DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/firkode/C_serial)
  endif()

endforeach(example_tuple ${FIRKODE_examples})

if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  sundials_install_examples(
    firkode FIRKODE_examples
    CMAKE_TEMPLATE cmakelists_C_ex.in
    MAKE_TEMPLATE makefile_serial_C_ex.in
    SOLVER_LIBRARY sundials_firkode
    SUNDIALS_TARGETS nvecserial firkode
    DESTINATION firkode/C_serial
    EXTRA_FILES README ${FIRKODE_extras}
    TEST_INSTALL C_serial)
endif()
