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
# CMakeLists.txt file for the FIRKODE library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall FIRKODE\n\")")

# Add variable firkode_SOURCES with the sources for the FIRKODE library
set(firkode_SOURCES firkode.c firkode_io.c firkode_ls.c firkode_nls.c
                    firkode_root.c firkode_tables.c)

# Add variable firkode_HEADERS with the exported FIRKODE header files
set(firkode_HEADERS firkode.h firkode.hpp firkode_ls.h firkode_tables.h)

# Add prefix with complete path to the FIRKODE header files
add_prefix(${SUNDIALS_SOURCE_DIR}/include/firkode/ firkode_HEADERS)

# Create the library
sundials_add_library(
  sundials_firkode
  SOURCES ${firkode_SOURCES}
  HEADERS ${firkode_HEADERS}
  INCLUDE_SUBDIR firkode
  LINK_LIBRARIES PUBLIC sundials_core
  OBJECT_LIBRARIES
    sundials_sunmemsys_obj
    sundials_nvecserial_obj
    sundials_nvecmanyvector_obj
    sundials_sunadaptcontrollersoderlind_obj
    sundials_sunadaptcontrollerimexgus_obj
    sundials_sunmatrixband_obj
    sundials_sunmatrixdense_obj
    sundials_sunmatrixsparse_obj
    sundials_sunlinsolband_obj
    sundials_sunlinsoldense_obj
    sundials_sunlinsolspbcgs_obj
    sundials_sunlinsolspfgmr_obj
    sundials_sunlinsolspgmr_obj
    sundials_sunlinsolsptfqmr_obj
    sundials_sunlinsolpcg_obj
    sundials_sunnonlinsolauto_obj
    sundials_sunnonlinsolnewton_obj
    sundials_sunnonlinsolfixedpoint_obj
  OUTPUT_NAME sundials_firkode
  VERSION ${firkodelib_VERSION}
  SOVERSION ${firkodelib_SOVERSION})

# Finished FIRKODE
message(STATUS "Added FIRKODE module")
