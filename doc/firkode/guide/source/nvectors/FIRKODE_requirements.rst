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

.. _NVectors.FIRKODE:

NVECTOR functions used by FIRKODE
=================================

In :numref:`NVectors.FIRKODE.Table` below, we list the vector functions in
the NVECTOR module used within the FIRKODE package. The table also shows, for
each function, which of the code modules uses the function. The FIRKODE
column shows function usage within the main integrator module (including the
internal Newton iteration and the rootfinding module), while the FIRKLS
column shows function usage within the FIRKODE linear solver interface.

At this point, we should emphasize that the FIRKODE user does not need to
know anything about the usage of vector functions by the FIRKODE code modules
in order to use FIRKODE. The information is presented as an implementation
detail for the interested reader.

.. _NVectors.FIRKODE.Table:
.. table:: List of vector functions usage by FIRKODE code modules

   +--------------------------------+-------------+------------+
   |                                | FIRKODE     | FIRKLS     |
   +================================+=============+============+
   | :c:func:`N_VGetLength`         | x           | x          |
   +--------------------------------+-------------+------------+
   | :c:func:`N_VClone`             | x           | x          |
   +--------------------------------+-------------+------------+
   | :c:func:`N_VCloneEmpty`        |             | 1          |
   +--------------------------------+-------------+------------+
   | :c:func:`N_VDestroy`           | x           | x          |
   +--------------------------------+-------------+------------+
   | :c:func:`N_VGetArrayPointer`   |             | 1          |
   +--------------------------------+-------------+------------+
   | :c:func:`N_VSetArrayPointer`   |             | 1          |
   +--------------------------------+-------------+------------+
   | :c:func:`N_VLinearSum`         | x           | x          |
   +--------------------------------+-------------+------------+
   | :c:func:`N_VConst`             | x           | x          |
   +--------------------------------+-------------+------------+
   | :c:func:`N_VProd`              | x           |            |
   +--------------------------------+-------------+------------+
   | :c:func:`N_VDiv`               | x           |            |
   +--------------------------------+-------------+------------+
   | :c:func:`N_VScale`             | x           | x          |
   +--------------------------------+-------------+------------+
   | :c:func:`N_VAbs`               | x           |            |
   +--------------------------------+-------------+------------+
   | :c:func:`N_VInv`               | x           |            |
   +--------------------------------+-------------+------------+
   | :c:func:`N_VAddConst`          | x           |            |
   +--------------------------------+-------------+------------+
   | :c:func:`N_VDotProd`           | 2           | 2          |
   +--------------------------------+-------------+------------+
   | :c:func:`N_VMaxNorm`           | x           |            |
   +--------------------------------+-------------+------------+
   | :c:func:`N_VWrmsNorm`          | x           | x          |
   +--------------------------------+-------------+------------+
   | :c:func:`N_VMin`               | 3           |            |
   +--------------------------------+-------------+------------+
   | :c:func:`N_VLinearCombination` | x           | x          |
   +--------------------------------+-------------+------------+

Special cases (numbers match markings in table):

1. These routines are only required if an internal difference-quotient
   routine for constructing :ref:`SUNMATRIX_DENSE <SUNMatrix.Dense>` or
   :ref:`SUNMATRIX_BAND <SUNMatrix.Band>` Jacobian matrices is used.

2. This routine is required by the internal :ref:`SUNLINSOL_SPFGMR
   <SUNLinSol.SPFGMR>` iteration on the stacked stage system when the
   method has more than one stage.

3. This routine is only used when a vector of absolute tolerances with a
   zero component is supplied.

The stage increments are stored in an :ref:`NVECTOR_MANYVECTOR
<NVectors.ManyVector>` object whose subvectors are clones of the user's
vector. Consequently the user's vector implementation must provide the
operations required by that module (in particular :c:func:`N_VGetLength`),
and the fused operation :c:func:`N_VLinearCombination` is used whenever it is
provided (the generic fallback is used otherwise). Each SUNLINSOL object may
require additional NVECTOR routines not listed in the table above. Please see
the relevant descriptions of these modules in :numref:`SUNLinSol` for
additional detail on their NVECTOR requirements.
