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

.. _FIRKODE.Organization:

*****************
Code Organization
*****************

The FIRKODE package is written in ANSI C. The following summarizes the basic
structure of the package, although knowledge of this structure is not
necessary for its use.

The central integration module, implemented in the files ``firkode.h``,
``firkode_impl.h``, and ``firkode.c``, deals with the evaluation of the
stage equations, the prediction of the stage values, the estimation of the
local error, the selection of the step size, and the interpolation to user
output points, among other issues. The Radau IIA coefficient tables are
constructed in ``firkode_tables.c`` and exposed through the public header
``firkode_tables.h`` (see :numref:`FIRKODE.Usage.Tables`).

The nonlinear stage equations are solved by a simplified Newton iteration
implemented in ``firkode_nls.c``. The iteration operates on all stage
increments simultaneously, which are stored in a vector of the
:ref:`NVECTOR_MANYVECTOR <NVectors.ManyVector>` type built from
:math:`s` clones of the user's vector, and uses the
``SUNNONLINSOL_NEWTON`` module internally with the
convergence test of RADAU5. The nonlinear solver is not user-replaceable in
this version of FIRKODE.

The linear solver interface, FIRKLS, is implemented in ``firkode_ls.h``,
``firkode_ls_impl.h``, and ``firkode_ls.c``. It supports both direct and
iterative linear solvers built using the generic ``SUNLinearSolver`` API
(see :numref:`SUNLinSol`). These solvers may utilize a ``SUNMatrix``
object (see :numref:`SUNMatrix`) for storing Jacobian information, or they
may be matrix-free. The user's linear solver operates on the *block* system
:math:`M - \gamma J` of the size of the ODE system, exactly as in CVODE or
ARKODE. FIRKLS uses this block solver in two ways: as a block-diagonal
preconditioner for an internal :ref:`SUNLINSOL_SPFGMR <SUNLinSol.SPFGMR>`
iteration on the stacked Newton system of all stages, and directly for the
filtered error estimate. For users employing :ref:`SUNMATRIX_DENSE
<SUNMatrix.Dense>` or :ref:`SUNMATRIX_BAND <SUNMatrix.Band>` Jacobian
matrices, FIRKODE includes algorithms for their approximation through
difference quotients, although the user also has the option of supplying a
routine to compute the Jacobian (or an approximation to it) directly. This
user-supplied routine is required when using sparse or user-supplied Jacobian
matrices. For users employing matrix-free iterative linear solvers, FIRKODE
includes an algorithm for the approximation by difference quotients of the
product :math:`Jv`; again, the user has the option of providing routines for
this operation. For preconditioned iterative methods, the preconditioning must
be supplied by the user, in two phases: setup and solve.

FIRKLS also provides an interface for a constant, non-identity mass matrix
:math:`M`, which may be supplied either as a ``SUNMatrix`` object with a
matrix-based linear solver or through a user-supplied matrix-vector product
routine with a matrix-free linear solver.

The rootfinding module, implemented in ``firkode_root.c``, is adapted from
ARKODE and locates roots of user-supplied functions using the dense output
of the collocation polynomial.

All state information used by FIRKODE to solve a given problem is saved in a
structure, and a pointer to that structure is returned to the user. There is
no global data in the FIRKODE package, and so, in this respect, it is
reentrant. State information specific to the linear solver is saved in a
separate structure, a pointer to which resides in the FIRKODE memory
structure.
