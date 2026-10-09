.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2025-2026, Lawrence Livermore National Security,
   University of Maryland Baltimore County, and the SUNDIALS contributors.
   Copyright (c) 2013-2025, Lawrence Livermore National Security
   and Southern Methodist University.
   Copyright (c) 2002-2013, Lawrence Livermore National Security.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _FIRKODE.Usage:

*************
Using FIRKODE
*************

This chapter is concerned with the use of FIRKODE for the solution of initial
value problems (IVPs) in C and C++ applications. The chapter builds upon
:numref:`SUNDIALS`. The following sections treat the header files and the
layout of the user's main program, and provide descriptions of the FIRKODE
user-callable functions and user-supplied functions.

The example programs in the ``examples/firkode`` directory of the SUNDIALS
source tree may also be helpful. Those codes may be used as templates (with
the removal of some lines used in testing) and are included in the FIRKODE
package.

The user should be aware that not all ``SUNLinearSolver`` and ``SUNMatrix``
modules are compatible with all ``N_Vector`` implementations. Details on
compatibility are given in the documentation for each ``SUNMatrix`` module
(:numref:`SUNMatrix`) and each ``SUNLinearSolver`` module
(:numref:`SUNLinSol`). For example, ``NVECTOR_PARALLEL`` is not compatible
with the dense, banded, or sparse ``SUNMatrix`` types, or with the
corresponding dense, banded, or sparse ``SUNLinearSolver`` modules. Please
check :numref:`SUNMatrix` and :numref:`SUNLinSol` to verify compatibility
between these modules. In addition, FIRKODE stores the stage increments of the
implicit method in an :ref:`NVECTOR_MANYVECTOR <NVectors.ManyVector>` object
built from clones of the user's vector, so the user's vector implementation
must provide the operations listed in :numref:`NVectors.FIRKODE`.

FIRKODE uses various constants for both input and output. These are defined
as needed in this chapter, but for convenience are also listed separately in
:numref:`FIRKODE.Constants`.

.. toctree::
   :maxdepth: 1

   General
   Skeleton
   User_callable
   User_supplied
   Tables
