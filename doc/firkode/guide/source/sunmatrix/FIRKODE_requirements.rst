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

.. _SUNMatrix.FIRKODE:

SUNMatrix functions used by FIRKODE
===================================

In :numref:`SUNMatrix.FIRKODE.Table`, we list the matrix functions in the
``SUNMatrix`` module used within the FIRKODE package. The main FIRKODE
integrator does not call any ``SUNMatrix`` functions directly, so the table
is specific to the FIRKLS interface. We further note that the FIRKLS
interface only utilizes these routines when supplied with a *matrix-based*
linear solver, i.e., the ``SUNMatrix`` object passed to
:c:func:`FIRKodeSetLinearSolver` or :c:func:`FIRKodeSetMassLinearSolver` was
not ``NULL``.

At this point, we should emphasize that the FIRKODE user does not need to
know anything about the usage of matrix functions by the FIRKODE code modules
in order to use FIRKODE. The information is presented as an implementation
detail for the interested reader.

.. _SUNMatrix.FIRKODE.Table:
.. table:: List of matrix functions usage by FIRKODE code modules
   :align: center

   +-------------------------------+-----------------+
   |                               |      FIRKLS     |
   +===============================+=================+
   | :c:func:`SUNMatClone`         | x               |
   +-------------------------------+-----------------+
   | :c:func:`SUNMatDestroy`       | x               |
   +-------------------------------+-----------------+
   | :c:func:`SUNMatZero`          | x               |
   +-------------------------------+-----------------+
   | :c:func:`SUNMatGetID`         | x               |
   +-------------------------------+-----------------+
   | :c:func:`SUNMatCopy`          | x               |
   +-------------------------------+-----------------+
   | :c:func:`SUNMatScaleAddI`     | 1               |
   +-------------------------------+-----------------+
   | :c:func:`SUNMatScaleAdd`      | 2               |
   +-------------------------------+-----------------+
   | :c:func:`SUNMatMatvec`        | 3               |
   +-------------------------------+-----------------+
   | :c:func:`SUNMatMatvecSetup`   | :math:`\dagger` |
   +-------------------------------+-----------------+

The matrix functions listed with a :math:`\dagger` symbol are optionally
used, in that these are only called if they are implemented in the
``SUNMatrix`` module that is being used (i.e. their function pointers are
non-``NULL``). Special cases (numbers match markings in table):

1. This routine is used to form the block system matrix
   :math:`I - \gamma J` when no mass matrix is supplied.

2. This routine is used to form the block system matrix
   :math:`M - \gamma J` when a mass matrix is supplied; the Jacobian and
   mass matrices must then have the same ``SUNMatrix`` type.

3. This routine is used for the Jacobian-vector products
   :math:`J v` required by the stacked stage system when the method has
   more than one stage, and for mass matrix-vector products :math:`M v`
   when a matrix-based mass linear solver is used, unless the user supplies
   the corresponding routine (:c:func:`FIRKodeSetJacTimes` or
   :c:func:`FIRKodeSetMassTimes`).
