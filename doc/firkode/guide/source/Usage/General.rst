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

.. _FIRKODE.Usage.General:

Access to library and header files
==================================

At this point, it is assumed that the installation of FIRKODE, following the
procedure described in :numref:`Installation`, has been completed
successfully. In the proceeding text, the directories ``libdir`` and
``incdir`` are the installation library and include directories,
respectively. For a default installation, these are ``instdir/lib`` and
``instdir/include``, respectively, where ``instdir`` is the directory
where SUNDIALS was installed.

Regardless of where the user's application program resides, its associated
compilation and load commands must make reference to the appropriate
locations for the library and header files required by FIRKODE. FIRKODE
symbols are found in ``libdir/libsundials_firkode.lib``. Thus, in addition to
linking to ``libdir/libsundials_core.lib``, FIRKODE users need to link to the
FIRKODE library. Symbols for additional SUNDIALS modules, vectors and
algebraic solvers, are found in

.. code-block:: none

   <libdir>/libsundials_nvec*.lib
   <libdir>/libsundials_sunmat*.lib
   <libdir>/libsundials_sunlinsol*.lib
   <libdir>/libsundials_sunadaptcontroller*.lib

The file extension ``.lib`` is typically ``.so`` for shared libraries and
``.a`` for static libraries. Applications built with CMake may instead link
to the imported target ``SUNDIALS::firkode`` provided by the SUNDIALS CMake
package configuration.

The relevant header files for FIRKODE are located in the subdirectory
``incdir/include/firkode``. To use FIRKODE the application needs to include
the header file for the main integrator, which also includes the header file
of the linear solver interface and of the coefficient tables:

.. code-block:: c

   #include <firkode/firkode.h>

The calling program must also include an :c:type:`N_Vector` implementation
header file, of the form ``nvector/nvector_*.h``. See :numref:`NVectors`
for the appropriate name. A linear solver module header file of the form
``sunlinsol/sunlinsol_*.h``, where ``*`` is the name of the linear
solver module, is always required since FIRKODE requires a linear solver (see
:numref:`SUNLinSol` for the appropriate name). If the linear solver is
matrix-based, the linear solver header will also include a header file of the
form ``sunmatrix/sunmatrix_*.h`` where ``*`` is the name of the matrix
implementation compatible with the linear solver. If a non-default time step
adaptivity controller is used, the corresponding header file of the form
``sunadaptcontroller/sunadaptcontroller_*.h`` is also required.

.. _FIRKODE.Usage.DataTypes:

Data types
==========

FIRKODE uses the SUNDIALS data types ``sunrealtype``, ``sunindextype``, and
``sunbooleantype`` described in :numref:`SUNDIALS.DataTypes`. All FIRKODE
functions take a pointer to the FIRKODE memory block, of type ``void*``, as
their first argument; this pointer is returned by :c:func:`FIRKodeCreate`.
