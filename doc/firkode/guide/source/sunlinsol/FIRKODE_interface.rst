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

.. _SUNLinSol.FIRKODE:

FIRKODE SUNLinearSolver interface
=================================

:numref:`SUNLinSol.FIRKODE.Table` below lists the ``SUNLinearSolver`` module
linear solver functions used within the FIRKLS interface. As with the
``SUNMatrix`` module, we emphasize that the FIRKODE user does not need to
know detailed usage of linear solver functions by the FIRKODE code modules in
order to use FIRKODE. The information is presented as an implementation
detail for the interested reader.

The linear solver functions listed below are marked with 'x' to indicate that
they are required, or with :math:`\dagger` to indicate that they are only
called if they are non-``NULL`` in the ``SUNLinearSolver`` implementation
that is being used. Note:

#. ``SUNLinSolNumIters`` is only used to accumulate overall iterative linear
   solver statistics. If it is not implemented by the ``SUNLinearSolver``
   module, then FIRKLS will consider all solves as requiring zero iterations.

#. Although FIRKLS does not call ``SUNLinSolLastFlag`` directly, this
   routine is available for users to query linear solver issues directly.

#. Although FIRKLS does not call ``SUNLinSolFree`` directly, this routine
   should be available for users to call when cleaning up from a simulation.

.. _SUNLinSol.FIRKODE.Table:
.. table:: List of linear solver function usage in the FIRKLS interface

   +----------------------------------------+-----------------+-----------------+------------------+
   |                                        |     DIRECT      |    ITERATIVE    | MATRIX_ITERATIVE |
   +========================================+=================+=================+==================+
   | :c:func:`SUNLinSolGetType`             | x               | x               | x                |
   +----------------------------------------+-----------------+-----------------+------------------+
   | :c:func:`SUNLinSolSetATimes`           | :math:`\dagger` | x               | :math:`\dagger`  |
   +----------------------------------------+-----------------+-----------------+------------------+
   | :c:func:`SUNLinSolSetPreconditioner`   | :math:`\dagger` | :math:`\dagger` | :math:`\dagger`  |
   +----------------------------------------+-----------------+-----------------+------------------+
   | :c:func:`SUNLinSolSetScalingVectors`   | :math:`\dagger` | :math:`\dagger` | :math:`\dagger`  |
   +----------------------------------------+-----------------+-----------------+------------------+
   | :c:func:`SUNLinSolSetZeroGuess`        | :math:`\dagger` | :math:`\dagger` | :math:`\dagger`  |
   +----------------------------------------+-----------------+-----------------+------------------+
   | :c:func:`SUNLinSolInitialize`          | x               | x               | x                |
   +----------------------------------------+-----------------+-----------------+------------------+
   | :c:func:`SUNLinSolSetup`               | x               | x               | x                |
   +----------------------------------------+-----------------+-----------------+------------------+
   | :c:func:`SUNLinSolSolve`               | x               | x               | x                |
   +----------------------------------------+-----------------+-----------------+------------------+
   | :math:`^1` :c:func:`SUNLinSolNumIters` |                 | :math:`\dagger` | :math:`\dagger`  |
   +----------------------------------------+-----------------+-----------------+------------------+
   | :math:`^2` :c:func:`SUNLinSolLastFlag` |                 |                 |                  |
   +----------------------------------------+-----------------+-----------------+------------------+
   | :math:`^3` :c:func:`SUNLinSolFree`     |                 |                 |                  |
   +----------------------------------------+-----------------+-----------------+------------------+

Since there are a wide range of potential ``SUNLinearSolver`` use cases, the
following subsections describe some details of the FIRKLS interface, in the
case that interested users wish to develop custom ``SUNLinearSolver`` modules.

.. _SUNLinSol.FIRKODE.Block:

Block systems and the stacked stage system
------------------------------------------

The user's ``SUNLinearSolver`` object is used by FIRKLS to solve *block*
systems of the form :math:`(M - \gamma J)\, x = b`, where :math:`M` is the
mass matrix (the identity by default), :math:`J` is the Jacobian of the ODE
right-hand side, and :math:`\gamma = h \gamma_0` with :math:`\gamma_0`
a method constant (see :numref:`FIRKODE.Mathematics.Newton`). These are
systems of the size of the ODE, exactly as in CVODE and ARKODE. The block
solves are used in two ways:

* as a block-diagonal preconditioner for the internal
  :ref:`SUNLINSOL_SPFGMR <SUNLinSol.SPFGMR>` iteration that solves the
  stacked Newton system for all stages simultaneously, and

* directly, for the filtered local error estimate.

For a one-stage method (backward Euler) the stacked system coincides with
the block system and the user's solver is used directly. The flexible GMRES
iteration is used for the stacked system precisely because the block solves
may be inexact and may vary between iterations when the user's solver is
itself iterative. When the user's linear solver is matrix-free, its
preconditioner solve routine is applied block-wise directly rather than
nesting a second Krylov iteration.

.. _SUNLinSol.FIRKODE.Lagged:

Lagged matrix information
-------------------------

If the ``SUNLinearSolver`` object supplied to FIRKODE is matrix-based, then
FIRKLS maintains a copy of the Jacobian :math:`J` separately from the linear
system matrix :math:`M - \gamma J`, so that the system matrix may be
rebuilt when only :math:`\gamma` has changed. The Jacobian is evaluated at
the beginning of a step, at the point :math:`(t_n, y_n)`, and is reused
across subsequent steps until FIRKODE determines that it should be updated
(see :c:func:`FIRKodeSetJacEvalFrequency` and
:c:func:`FIRKodeSetDeltaGammaMaxBadJac`). The saved Jacobian is also used for
the Jacobian-vector products required by the stacked stage system when the
matrix implements :c:func:`SUNMatMatvec` and the user has not supplied a
Jacobian-times-vector routine.

.. _SUNLinSol.FIRKODE.Iterative.Tolerance:

Iterative linear solver tolerance
---------------------------------

If the ``SUNLinearSolver`` object supplied to FIRKODE is iterative, then
the block solves are performed to a tolerance :math:`\epsilon_L = \epsilon\,
\epsilon_{nls}\, n_f` in the weighted 2-norm, where :math:`\epsilon_{nls}`
is the nonlinear solver convergence coefficient (see
:c:func:`FIRKodeSetNonlinConvCoef`), :math:`\epsilon` is the factor set by
:c:func:`FIRKodeSetEpsLin`, and :math:`n_f` is the norm conversion factor
set by :c:func:`FIRKodeSetLSNormFactor`. For the error estimate, whose
right-hand side :math:`b` may be small, the tolerance is
:math:`\epsilon\, \min(\epsilon_{nls}, \|b\|)\, n_f`, so that a small
estimate is computed rather than taken as zero; block solves that
precondition the stacked stage system use the tolerance
:math:`\epsilon\, \|b\|\, n_f` relative to their right-hand side. The
corresponding tolerance for the stacked stage system is controlled by
:c:func:`FIRKodeSetStageEpsLin`.
