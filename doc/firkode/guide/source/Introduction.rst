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

.. _FIRKODE.Introduction:

************
Introduction
************

FIRKODE is a solver for stiff initial value problems (IVPs) for systems of
ordinary differential equations (ODEs) based on the fully implicit Radau IIA
Runge--Kutta methods. It is part of a software family called SUNDIALS: SUite
of Nonlinear and DIfferential/ALgebraic equation Solvers
:cite:p:`HBGLSSW:05`. This suite consists of CVODE, ARKODE, KINSOL, IDA,
FIRKODE, and variants of these with sensitivity analysis capabilities.

.. _FIRKODE.Introduction.background:

Background
==========

The Radau IIA methods :cite:p:`HaWa:91` are the collocation methods based on
the Radau quadrature nodes. The :math:`s`-stage method has classical order
:math:`2s-1` and stage order :math:`s`, and is A-stable, L-stable, and
stiffly accurate. Because the stage order does not degrade for very stiff
problems and the methods damp stiff components in a single step, they are
among the most robust choices for very stiff ODEs and for problems that are
close to differential-algebraic systems. The most widely used implementation
is the Fortran code RADAU5 :cite:p:`HaWa:91,HaWa:99`, which implements the
three-stage, fifth-order method with a simplified Newton iteration, a
filtered error estimate, and a predictive step size controller.

FIRKODE brings these methods into the SUNDIALS framework. The stage equations
are solved with a simplified Newton iteration on the *stacked* system of all
stages, formulated on an :c:type:`N_Vector` of the SUNDIALS many-vector
type, and the resulting linear systems are solved by a flexible Krylov
iteration in which the user's :c:type:`SUNLinearSolver` is applied
block-wise as a preconditioner. In this way all SUNDIALS linear solvers
(dense, banded, sparse, and Krylov) and all vector implementations may be
used with the fully implicit methods, and the additional work compared to a
diagonally implicit method is confined to a small number of matrix-vector
products with the Jacobian. FIRKODE provides the Radau IIA methods with
:math:`s = 1, 2, 3` stages in closed form and constructs the tables for
:math:`4 \le s \le 9` at run time.

The organization of FIRKODE follows that of CVODE and ARKODE: the linear
solver interface, the rootfinding module, and the general structure of the
integrator are adapted from those packages. As with all SUNDIALS packages,
FIRKODE is written in terms of the generic :c:type:`N_Vector`,
:c:type:`SUNMatrix`, and :c:type:`SUNLinearSolver` operations, so that it
may be used with any of the vector, matrix, and linear solver implementations
supplied with SUNDIALS or with user-supplied implementations.

Changes to SUNDIALS in release X.Y.Z
====================================

.. include:: ../../../shared/RecentChanges.rst

For changes in prior versions of SUNDIALS see :numref:`Changelog`.

.. _FIRKODE.Introduction.reading:

Reading this User Guide
=======================

This user guide is a combination of general usage instructions. Specific
example programs are provided in the ``examples/firkode`` directory of the
SUNDIALS source tree. We expect that some readers will want to concentrate on
the general instructions, while others will refer mostly to the examples, and
the organization is intended to accommodate both styles.

The structure of this document is as follows:

-  In :numref:`FIRKODE.Mathematics`, we give short descriptions of the
   numerical methods implemented by FIRKODE for the solution of initial
   value problems for systems of ODEs, including the Radau IIA methods, the
   nonlinear and linear solution strategies, error estimation and step size
   control, dense output, rootfinding, and the treatment of non-identity mass
   matrices.

-  The following chapter describes the software organization of the FIRKODE
   solver (:numref:`FIRKODE.Organization`).

-  :numref:`FIRKODE.Usage` is the main usage document for FIRKODE for C
   applications. It includes a complete description of the user interface
   for the integration of ODE initial value problems.

-  :numref:`NVectors` gives a brief overview of the generic ``N_Vector``
   module shared among the various components of SUNDIALS, and details on the
   ``N_Vector`` implementations provided with SUNDIALS.

-  :numref:`SUNMatrix` gives a brief overview of the generic ``SUNMatrix``
   module shared among the various components of SUNDIALS, and details on the
   ``SUNMatrix`` implementations provided with SUNDIALS.

-  :numref:`SUNLinSol` gives a brief overview of the generic
   ``SUNLinearSolver`` module shared among the various components of
   SUNDIALS, and details on the ``SUNLinearSolver`` implementations
   provided with SUNDIALS.

-  :numref:`SUNAdaptController` describes the time step adaptivity
   controllers shared among the SUNDIALS time integrators.

-  Finally, in the appendices, we provide detailed instructions for the
   installation of FIRKODE, within the structure of SUNDIALS
   (:numref:`Installation`), as well as a list of all the constants used for
   input to and output from FIRKODE functions (:numref:`FIRKODE.Constants`).

Finally, the reader should be aware of the following notational conventions
in this user guide: program listings and identifiers (such as
:c:func:`FIRKodeInit`) within textual explanations are hyperlinked to their
definitions directly; fields in C structures (such as *content*) appear in
italics; and packages or modules, such as FIRKLS, are written in all
capitals.


SUNDIALS License and Notices
============================

.. ifconfig:: package_name != 'super'

   .. include:: ../../../shared/LicenseReleaseNumbers.rst

.. ifconfig:: package_name == 'super'

   All SUNDIALS packages are released open source, under the BSD 3-Clause
   license for more details see the LICENSE and NOTICE files provided with all
   SUNDIALS packages.
