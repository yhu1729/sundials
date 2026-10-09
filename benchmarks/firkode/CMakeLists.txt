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
# CMakeLists.txt file for the FIRKODE benchmarks: firk_klin (the
# stacked Krylov iteration counts) and firkode_benchmark (the serial
# harness comparing FIRKODE with CVODE, ARKODE and IDA).
# ---------------------------------------------------------------

if(SUNDIALS_ENABLE_FIRKODE AND CXX_FOUND)

  # ---- firk_klin ----

  sundials_add_executable(firk_klin firk_klin.cpp)
  set_target_properties(firk_klin PROPERTIES FOLDER "Benchmarks")
  target_link_libraries(
    firk_klin
    PRIVATE sundials_firkode
            sundials_nvecserial
            sundials_sunmatrixband
            sundials_sunlinsolband
            sundials_sunlinsoldense
            sundials_sunmatrixdense
            sundials_sunlinsolpcg
            sundials_sunlinsolspgmr
            ${EXE_EXTRA_LINK_LIBS})

  install(TARGETS firk_klin
          DESTINATION "${SUNDIALS_BENCHMARKS_INSTALL_PATH}/firkode")

  # ---- firkode_benchmark ----

  # The harness compares against tight references in double precision; the
  # extended-precision build only computes such references.
  if(SUNDIALS_PRECISION MATCHES "DOUBLE|EXTENDED")

    set(sources
        firkode_benchmark.cpp
        problems_small.cpp
        problem_diurnal.cpp
        problem_heat2d.cpp
        problem_heat1d_fem.cpp
        problem_sector.cpp
        problem_wave1d_fem.cpp
        problems_testset.cpp
        problem_bruss3d.cpp
        integrator_firkode.cpp)

    set(libs
        sundials_firkode
        sundials_nvecserial
        sundials_sunlinsoldense
        sundials_sunlinsolband
        sundials_sunlinsolspgmr
        sundials_sunlinsolpcg
        sundials_sunmatrixdense
        sundials_sunmatrixband)

    set(defs)

    if(SUNDIALS_ENABLE_CVODE)
      list(APPEND sources integrator_cvode.cpp)
      list(APPEND libs sundials_cvode)
      list(APPEND defs BENCHMARK_HAVE_CVODE)
    endif()

    if(SUNDIALS_ENABLE_ARKODE)
      list(APPEND sources integrator_arkode.cpp)
      list(APPEND libs sundials_arkode)
      list(APPEND defs BENCHMARK_HAVE_ARKODE)
    endif()

    if(SUNDIALS_ENABLE_IDA)
      list(APPEND sources integrator_ida.cpp)
      list(APPEND libs sundials_ida)
      list(APPEND defs BENCHMARK_HAVE_IDA)
    endif()

    if(SUNDIALS_ENABLE_SUNLINSOL_KLU)
      list(APPEND libs sundials_sunlinsolklu)
      list(APPEND defs BENCHMARK_HAVE_KLU)
    endif()

    sundials_add_executable(firkode_benchmark ${sources})

    set_target_properties(firkode_benchmark PROPERTIES FOLDER "Benchmarks")

    target_compile_definitions(firkode_benchmark PRIVATE ${defs})

    # problems/kpr.hpp and problems/pr.hpp are shared with the unit tests
    target_include_directories(
      firkode_benchmark PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
                            ${SUNDIALS_SOURCE_DIR}/test/unit_tests)

    target_link_libraries(firkode_benchmark PRIVATE ${libs} ${EXE_EXTRA_LINK_LIBS})

    add_dependencies(benchmark firkode_benchmark)

    install(TARGETS firkode_benchmark
            DESTINATION "${SUNDIALS_BENCHMARKS_INSTALL_PATH}/firkode")

    message(STATUS "Added FIRKODE benchmark harness firkode_benchmark")

  else()
    message(STATUS "Skipping firkode_benchmark (requires double or extended precision)")
  endif()

  install(FILES README.md
          DESTINATION "${SUNDIALS_BENCHMARKS_INSTALL_PATH}/firkode")

endif()
