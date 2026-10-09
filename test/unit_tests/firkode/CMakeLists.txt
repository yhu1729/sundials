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
# unit_tests/firkode level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

# Creates the executable of a FIRKODE unit test. The test links the FIRKODE
# object files and the object libraries bundled into the FIRKODE library rather
# than the library itself, so that it can call internal functions on every
# platform (a Windows DLL exports only the public API).
function(firkode_add_unit_test_target target)
  sundials_add_executable(${target} ${ARGN})
  set_target_properties(${target} PROPERTIES FOLDER "unit_tests")
  target_include_directories(
    ${target}
    PRIVATE $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
            ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/src
            ${CMAKE_SOURCE_DIR}/test/unit_tests)
  get_property(_objlibs GLOBAL PROPERTY SUNDIALS_FIRKODE_OBJECT_LIBRARIES)
  target_link_libraries(${target} PRIVATE $<TARGET_OBJECTS:sundials_firkode_obj>
                                          ${_objlibs} ${EXE_EXTRA_LINK_LIBS})
  # $<TARGET_OBJECTS:...> does not add a dependency on the object library
  add_dependencies(${target} sundials_firkode_obj)
endfunction()

# C unit tests
add_subdirectory(C_serial)
if(SUNDIALS_ENABLE_MPI AND SUNDIALS_ENABLE_NVECTOR_PARALLEL)
  add_subdirectory(C_parallel)
endif()

# C++ unit tests
if(CXX_FOUND)
  add_subdirectory(CXX_serial)
endif()
