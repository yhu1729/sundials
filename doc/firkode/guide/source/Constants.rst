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

.. _FIRKODE.Constants:

*****************
FIRKODE Constants
*****************

Below we list all input and output constants used by the main solver and
linear solver modules, together with their numerical values and a short
description of their meaning. :numref:`FIRKODE.Constants.in_constants`
contains the FIRKODE input constants, and
:numref:`FIRKODE.Constants.out_constants` contains the FIRKODE output
constants.


.. _FIRKODE.Constants.in_constants:
.. table:: FIRKODE input constants

   +------------------------------+-----+-------------------------------------------------------+
   | ``FIRK_NORMAL``              | 1   | Solver returns at specified output time.              |
   +------------------------------+-----+-------------------------------------------------------+
   | ``FIRK_ONE_STEP``            | 2   | Solver returns after each successful step.            |
   +------------------------------+-----+-------------------------------------------------------+
   | ``FIRK_BEGINNER``            | 0   | BEGINNER mode: the essential options only.            |
   +------------------------------+-----+-------------------------------------------------------+
   | ``FIRK_EXPERT``              | 1   | EXPERT mode: every option.                            |
   +------------------------------+-----+-------------------------------------------------------+
   | ``FIRK_PREDICT_EXTRAPOLATE`` | 0   | Predict the stage increments by extrapolation of the  |
   |                              |     | previous collocation polynomial.                      |
   +------------------------------+-----+-------------------------------------------------------+
   | ``FIRK_PREDICT_TRIVIAL``     | 1   | Predict zero stage increments.                        |
   +------------------------------+-----+-------------------------------------------------------+
   | ``FIRK_POLICY_RADAU5``       | 0   | RADAU5's Jacobian and factorization reuse policy.     |
   +------------------------------+-----+-------------------------------------------------------+
   | ``FIRK_POLICY_PERIODIC``     | 1   | CVODE's periodic reuse policy.                        |
   +------------------------------+-----+-------------------------------------------------------+
   | ``FIRK_MAX_STAGES``          | 9   | Maximum number of stages.                             |
   +------------------------------+-----+-------------------------------------------------------+
   | ``FIRK_DEFAULT_STAGES``      | 3   | Default number of stages.                             |
   +------------------------------+-----+-------------------------------------------------------+
   | ``SUN_PREC_NONE``            | 0   | No preconditioning                                    |
   +------------------------------+-----+-------------------------------------------------------+
   | ``SUN_PREC_LEFT``            | 1   | Preconditioning on the left only.                     |
   +------------------------------+-----+-------------------------------------------------------+
   | ``SUN_PREC_RIGHT``           | 2   | Preconditioning on the right only.                    |
   +------------------------------+-----+-------------------------------------------------------+
   | ``SUN_PREC_BOTH``            | 3   | Preconditioning on both the left and the right.       |
   +------------------------------+-----+-------------------------------------------------------+
   | ``SUN_MODIFIED_GS``          | 1   | Use modified Gram-Schmidt procedure.                  |
   +------------------------------+-----+-------------------------------------------------------+
   | ``SUN_CLASSICAL_GS``         | 2   | Use classical Gram-Schmidt procedure.                 |
   +------------------------------+-----+-------------------------------------------------------+


.. _FIRKODE.Constants.out_constants:
.. table:: FIRKODE output constants
   :widths: 30 5 60

   +-------------------------------+-----+---------------------------------------------------------------------+
   | **General outputs**                                                                                       |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRK_SUCCESS``              | 0   | Successful function return.                                         |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRK_TSTOP_RETURN``         | 1   | :c:func:`FIRKodeEvolve` succeeded by reaching the specified         |
   |                               |     | stopping point.                                                     |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRK_ROOT_RETURN``          | 2   | :c:func:`FIRKodeEvolve` succeeded and found one or more roots.      |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRK_WARNING``              | 99  | The code of a warning sent to the logger (a step with               |
   |                               |     | :math:`t + h = t`, root functions identically zero); never returned.|
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRK_TOO_MUCH_WORK``        | -1  | The solver took ``mxstep`` internal steps but could not reach       |
   |                               |     | ``tout``.                                                           |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRK_TOO_MUCH_ACC``         | -2  | The solver could not satisfy the accuracy demanded by the user for  |
   |                               |     | some internal step.                                                 |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRK_ERR_FAILURE``          | -3  | Error test failures occurred too many times during one internal     |
   |                               |     | time step, or with :math:`|h| = h_{min}`.                           |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRK_CONV_FAILURE``         | -4  | Convergence test failures occurred too many times during one        |
   |                               |     | internal time step, or with :math:`|h| = h_{min}`.                  |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRK_LINIT_FAIL``           | -5  | The linear solver's initialization function failed.                 |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRK_LSETUP_FAIL``          | -6  | The linear solver's setup function failed in an unrecoverable       |
   |                               |     | manner.                                                             |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRK_LSOLVE_FAIL``          | -7  | The linear solver's solve function failed in an unrecoverable       |
   |                               |     | manner.                                                             |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRK_RHSFUNC_FAIL``         | -8  | The right-hand side function failed in an unrecoverable manner.     |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRK_FIRST_RHSFUNC_ERR``    | -9  | The right-hand side function failed at the first call.              |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRK_REPTD_RHSFUNC_ERR``    | -10 | The right-hand side function had repeated recoverable errors.       |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRK_UNREC_RHSFUNC_ERR``    | -11 | The right-hand side function had a recoverable error, but no        |
   |                               |     | recovery is possible.                                               |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRK_RTFUNC_FAIL``          | -12 | The rootfinding function failed in an unrecoverable manner.         |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRK_MEM_FAIL``             | -20 | A memory allocation failed.                                         |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRK_MEM_NULL``             | -21 | The ``firkode_mem`` argument was ``NULL``.                          |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRK_ILL_INPUT``            | -22 | One of the function inputs is illegal.                              |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRK_NO_MALLOC``            | -23 | The FIRKODE memory block was not allocated by a call to             |
   |                               |     | :c:func:`FIRKodeInit`.                                              |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRK_BAD_K``                | -24 | The derivative order :math:`k` is larger than allowed.              |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRK_BAD_T``                | -25 | The time :math:`t` is outside the last step taken.                  |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRK_BAD_DKY``              | -26 | The output derivative vector is ``NULL``.                           |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRK_TOO_CLOSE``            | -27 | The output and initial times are too close to each other.           |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRK_VECTOROP_ERR``         | -28 | A vector operation failed.                                          |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRK_CONTEXT_ERR``          | -32 | Reserved (:c:func:`FIRKodeCreate` returns ``NULL`` for a ``NULL``   |
   |                               |     | ``SUNContext``).                                                    |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRK_MASSINIT_FAIL``        | -40 | The mass matrix solver's initialization function failed.            |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRK_MASSSETUP_FAIL``       | -41 | The mass matrix solver's setup function failed.                     |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRK_MASSSOLVE_FAIL``       | -42 | The mass matrix solve failed in an unrecoverable manner.            |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRK_MASSFREE_FAIL``        | -43 | Reserved for the mass matrix solver's free function.                |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRK_MASSMULT_FAIL``        | -44 | The mass matrix-vector product failed.                              |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRK_INVALID_TABLE``        | -45 | A coefficient table is invalid (reserved for user tables).          |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRK_TABLE_FAIL``           | -46 | No built-in coefficient table exists for the number of stages.      |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRK_CONTROLLER_ERR``       | -47 | A call to the time step adaptivity controller failed.               |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRK_STEPPER_ERR``          | -48 | An error occurred in the SUNStepper interface (reserved).           |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRK_ORDER_CHANGE_FAIL``    | -49 | A change of the number of stages failed (reserved for order         |
   |                               |     | selection).                                                         |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRK_UNRECOGNIZED_ERR``     | -99 | An unrecognized error occurred.                                     |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | **FIRKLS linear solver interface outputs**                                                                |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRKLS_SUCCESS``            | 0   | Successful function return.                                         |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRKLS_MEM_NULL``           | -1  | The ``firkode_mem`` argument was ``NULL``.                          |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRKLS_LMEM_NULL``          | -2  | The FIRKLS linear solver interface has not been initialized.        |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRKLS_ILL_INPUT``          | -3  | The FIRKLS solver is not compatible with the current ``N_Vector``   |
   |                               |     | module, or an input value was illegal.                              |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRKLS_MEM_FAIL``           | -4  | A memory allocation request failed.                                 |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRKLS_MASSMEM_NULL``       | -6  | The mass matrix solver interface has not been initialized.          |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRKLS_JACFUNC_UNRECVR``    | -7  | The Jacobian function failed in an unrecoverable manner.            |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRKLS_JACFUNC_RECVR``      | -8  | The Jacobian function had a recoverable error.                      |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRKLS_MASSFUNC_UNRECVR``   | -9  | The mass matrix function failed in an unrecoverable manner.         |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRKLS_MASSFUNC_RECVR``     | -10 | The mass matrix function had a recoverable error.                   |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRKLS_SUNMAT_FAIL``        | -11 | An error occurred with the current ``SUNMatrix`` module.            |
   +-------------------------------+-----+---------------------------------------------------------------------+
   | ``FIRKLS_SUNLS_FAIL``         | -12 | An error occurred with the current ``SUNLinearSolver`` module.      |
   +-------------------------------+-----+---------------------------------------------------------------------+
