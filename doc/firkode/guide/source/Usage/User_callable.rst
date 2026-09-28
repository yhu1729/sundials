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

.. _FIRKODE.Usage.user_callable:

User-callable functions
=======================

This section describes the FIRKODE functions that are called by the user to
setup and then solve an IVP. Some of these are required. However, starting
with :numref:`FIRKODE.Usage.optional_input`, the functions listed involve
optional inputs/outputs or restarting, and those paragraphs may be skipped for
a casual use of FIRKODE. In any case, refer to :numref:`FIRKODE.Usage.Skeleton`
for the correct order of these calls.

On an error, each user-callable function returns a negative value and sends an
error message to the error handler routine, which prints the message on
``stderr`` by default. However, the user can set a file as error output or can
provide their own error handler function (see :numref:`SUNDIALS.Errors`).

.. _FIRKODE.Usage.init:

FIRKODE initialization and deallocation functions
-------------------------------------------------

The following three functions must be called in the order listed. The last
one is to be called only after the IVP solution is complete, as it frees the
FIRKODE memory block created and allocated by the first two calls.

.. c:function:: void* FIRKodeCreate(SUNContext sunctx)

   The function :c:func:`FIRKodeCreate` instantiates a FIRKODE solver object.

   **Arguments:**
      - ``sunctx`` -- the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   **Return value:**
      - If successful, :c:func:`FIRKodeCreate` returns a pointer to the newly
        created FIRKODE memory block (of type ``void *``). Otherwise, it
        returns ``NULL``.

.. c:function:: int FIRKodeInit(void* firkode_mem, FIRKRhsFn f, sunrealtype t0, N_Vector y0)

   The function :c:func:`FIRKodeInit` provides required problem and solution
   specifications, allocates internal memory, and initializes FIRKODE.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block returned by :c:func:`FIRKodeCreate`.
      - ``f`` -- is the C function which computes the right-hand side function :math:`f` in the ODE :math:`M \dot{y} = f(t,y)`. This function has the form ``f(t, y, ydot, user_data)`` (for full details see :c:type:`FIRKRhsFn`).
      - ``t0`` -- is the initial value of :math:`t`.
      - ``y0`` -- is the initial value of :math:`y`.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The call was successful.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized through a previous call to :c:func:`FIRKodeCreate`.
      - ``FIRK_MEM_FAIL`` -- A memory allocation request has failed.
      - ``FIRK_ILL_INPUT`` -- An input argument to :c:func:`FIRKodeInit` has an illegal value, or the vector ``y0`` does not implement a required operation (see :numref:`NVectors.FIRKODE`).

   **Notes:**
      If an error occurred, :c:func:`FIRKodeInit` also sends an error message
      to the error handler function.

.. c:function:: void FIRKodeFree(void** firkode_mem)

   The function :c:func:`FIRKodeFree` frees the memory allocated by a previous
   call to :c:func:`FIRKodeCreate`, including the memory of the FIRKLS linear
   solver interface, the mass matrix interface, and the rootfinding module
   (but not the ``SUNLinearSolver`` and ``SUNMatrix`` objects supplied by the
   user).

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.

   **Return value:**
      - The function :c:func:`FIRKodeFree` has no return value.


.. _FIRKODE.Usage.tolerances:

FIRKODE tolerance specification functions
-----------------------------------------

One of the following three functions must be called to specify the
integration tolerances (or directly specify the weights used in evaluating
WRMS vector norms). Note that this call must be made after the call to
:c:func:`FIRKodeInit`.

.. c:function:: int FIRKodeSStolerances(void* firkode_mem, sunrealtype reltol, sunrealtype abstol)

   The function :c:func:`FIRKodeSStolerances` specifies scalar relative and
   absolute tolerances.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``reltol`` -- is the scalar relative error tolerance.
      - ``abstol`` -- is the scalar absolute error tolerance.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The call was successful.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.
      - ``FIRK_NO_MALLOC`` -- The allocation function :c:func:`FIRKodeInit` has not been called.
      - ``FIRK_ILL_INPUT`` -- One of the input tolerances was negative.

.. c:function:: int FIRKodeSVtolerances(void* firkode_mem, sunrealtype reltol, N_Vector abstol)

   The function :c:func:`FIRKodeSVtolerances` specifies a scalar relative
   tolerance and vector absolute tolerances.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``reltol`` -- is the scalar relative error tolerance.
      - ``abstol`` -- is the vector of absolute error tolerances.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The call was successful.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.
      - ``FIRK_NO_MALLOC`` -- The allocation function :c:func:`FIRKodeInit` has not been called.
      - ``FIRK_ILL_INPUT`` -- The relative error tolerance was negative or the absolute tolerance had a negative component.

   **Notes:**
      This choice of tolerances is important when the absolute error
      tolerance needs to be different for each component of the state vector
      :math:`y`.

.. c:function:: int FIRKodeWFtolerances(void* firkode_mem, FIRKEwtFn efun)

   The function :c:func:`FIRKodeWFtolerances` specifies a user-supplied
   function ``efun`` that sets the multiplicative error weights :math:`W_i`
   for use in the weighted RMS norm, which are normally defined by
   :eq:`FIRKODE_errwt`.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``efun`` -- is the C function which defines the ``ewt`` vector (see :c:type:`FIRKEwtFn`).

   **Return value:**
      - ``FIRK_SUCCESS`` -- The call was successful.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.
      - ``FIRK_NO_MALLOC`` -- The allocation function :c:func:`FIRKodeInit` has not been called.
      - ``FIRK_ILL_INPUT`` -- The function ``efun`` is ``NULL``.

**General advice on the choice of tolerances.** For many users, the
appropriate choices for tolerance values in ``reltol`` and ``abstol`` are a
concern. The relative tolerance ``reltol`` controls the relative error in
each solution component and should not be smaller than about
:math:`10^{-15}` in double precision; the absolute tolerances control the
error in components that are small relative to their typical magnitude, and
should be chosen to reflect the level at which each component becomes
physically insignificant. The tolerances are used both in the local error
test and in the convergence test of the Newton iteration, where the
nonlinear residual is required to be smaller than the local error tolerance
by the factor set with :c:func:`FIRKodeSetNonlinConvCoef`. Since the
tolerances control local errors, the global error at the end of a long
integration is typically larger than the tolerances suggest; the actual
error can be assessed by comparing runs with different tolerances.


.. _FIRKODE.Usage.method:

Method selection functions
--------------------------

By default FIRKODE uses the three-stage Radau IIA method of order five. A
method with a different number of stages may be selected with one of the
following functions, which must be called before the first call to
:c:func:`FIRKodeEvolve` or directly after :c:func:`FIRKodeReInit`.

.. c:function:: int FIRKodeSetNumStages(void* firkode_mem, int s)

   The function :c:func:`FIRKodeSetNumStages` selects the :math:`s`-stage
   Radau IIA method, of order :math:`2s-1`.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``s`` -- the number of stages, :math:`1 \le s \le` ``FIRK_MAX_STAGES`` (= 9).

   **Return value:**
      - ``FIRK_SUCCESS`` -- The call was successful.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.
      - ``FIRK_ILL_INPUT`` -- The number of stages is out of range, or the integration has already started and the number of stages would change.

   **Notes:**
      The tables for :math:`s = 1, 2, 3` are hard-coded; the tables for
      :math:`s \ge 4` are constructed at run time (see
      :numref:`FIRKODE.Usage.Tables`). The default is :math:`s = 3`. Larger
      numbers of stages increase the storage of the internal Krylov solver
      for the coupled stage system (see :c:func:`FIRKodeSetStageSolverMaxl`).

.. c:function:: int FIRKodeSetOrder(void* firkode_mem, int ord)

   The function :c:func:`FIRKodeSetOrder` selects the Radau IIA method with
   the smallest number of stages whose order is at least ``ord``, i.e.,
   :math:`s = \lceil (\mathrm{ord}+1)/2 \rceil`.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``ord`` -- the requested order of accuracy; a non-positive value selects the default method.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The call was successful.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.
      - ``FIRK_ILL_INPUT`` -- The requested order exceeds that of the method with ``FIRK_MAX_STAGES`` stages, or the integration has already started and the number of stages would change.

.. c:function:: int FIRKodeGetNumStages(void* firkode_mem, int* s)

   The function :c:func:`FIRKodeGetNumStages` returns the number of stages
   of the current method.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``s`` -- the number of stages.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The call was successful.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

.. c:function:: int FIRKodeGetCurrentTable(void* firkode_mem, FIRKodeTable* T)

   The function :c:func:`FIRKodeGetCurrentTable` returns a pointer to the
   coefficient table of the current method (see :numref:`FIRKODE.Usage.Tables`).

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``T`` -- the coefficient table.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The call was successful.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.
      - ``FIRK_ILL_INPUT`` -- The table has not been constructed yet, i.e., :c:func:`FIRKodeEvolve` has not been called since the method was selected.

   **Notes:**
      The returned table is owned by FIRKODE and must not be modified or
      freed by the user.


.. _FIRKODE.Usage.linear_solver:

Linear solver interface function
--------------------------------

FIRKODE solves the coupled nonlinear stage equations with a simplified Newton
iteration (see :numref:`FIRKODE.Mathematics.Newton`), whose linear systems are
solved by an internal flexible Krylov iteration that is preconditioned
block-wise with the user's linear solver, applied to matrices of the form

.. math::
   \mathcal{A} = M - \gamma J \, , \qquad J = \frac{\partial f}{\partial y} \, ,
   :label: FIRKODE_Newton_block

where :math:`M` is the mass matrix and :math:`\gamma` is a scalar
proportional to the step size (see :numref:`FIRKODE.Mathematics.linear`).
The same block solver is used for the filtered local error estimate.
Consequently a ``SUNLinearSolver`` object *must* be attached to FIRKODE
before the first call to :c:func:`FIRKodeEvolve`. FIRKODE supports all of
the linear solver modules described in :numref:`SUNLinSol`: direct solvers
(which require a ``SUNMatrix`` object for :math:`\mathcal{A}`), matrix-based
iterative solvers, and matrix-free iterative solvers (for which the user may
supply a Jacobian-times-vector routine and a preconditioner). The user must
create the linear solver object (and, if applicable, the matrix object) and
then attach it with the following function.

.. c:function:: int FIRKodeSetLinearSolver(void* firkode_mem, SUNLinearSolver LS, SUNMatrix A)

   The function :c:func:`FIRKodeSetLinearSolver` attaches a generic
   ``SUNLinearSolver`` object ``LS`` and corresponding template Jacobian
   ``SUNMatrix`` object ``A`` (if applicable) to FIRKODE, initializing the
   FIRKLS linear solver interface.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``LS`` -- ``SUNLinearSolver`` object to use for solving linear systems of the form :eq:`FIRKODE_Newton_block`.
      - ``A`` -- ``SUNMatrix`` object used as a template for the Jacobian (or ``NULL`` if not applicable).

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The FIRKLS initialization was successful.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_ILL_INPUT`` -- The FIRKLS interface is not compatible with the ``LS`` or ``A`` input objects or is incompatible with the current ``N_Vector`` module, or :c:func:`FIRKodeInit` has not been called.
      - ``FIRKLS_SUNLS_FAIL`` -- A call to the ``LS`` object failed.
      - ``FIRKLS_MEM_FAIL`` -- A memory allocation request failed.

   **Notes:**
      If ``LS`` is a matrix-based linear solver, then the template Jacobian
      matrix ``A`` will be used in the solve process, so if additional storage
      is required within the ``SUNMatrix`` object (e.g. for factorization of a
      banded matrix), ensure that the input object is allocated with
      sufficient size (see :numref:`SUNMatrix` for further information).

      When using sparse linear solvers, it is typically much more efficient
      to supply ``A`` so that it includes the full sparsity pattern of the
      matrices :math:`M - \gamma J`, even if :math:`J` itself has zeros in
      nonzero locations of :math:`M`. The reasoning for this is that
      :math:`\mathcal{A}` is constructed in-place, on top of the
      user-specified values of ``A``.

      For methods with more than one stage the internal Krylov iteration on
      the coupled stage system requires matrix-vector products with the
      Jacobian. When ``LS`` is matrix-based these are computed with
      :c:func:`SUNMatMatvec` on the saved Jacobian matrix if the matrix
      implementation provides this operation, and by a difference quotient
      otherwise; in both cases a user-supplied Jacobian-times-vector routine
      may be attached with :c:func:`FIRKodeSetJacTimes` instead.

      The default difference-quotient Jacobian approximation is only
      available for the :ref:`SUNMATRIX_DENSE <SUNMatrix.Dense>` and
      :ref:`SUNMATRIX_BAND <SUNMatrix.Band>` matrix types; for other matrix
      types a Jacobian routine must be supplied with :c:func:`FIRKodeSetJacFn`.


.. _FIRKODE.Usage.mass_solver:

Mass matrix solver interface function
-------------------------------------

For problems :math:`M \dot{y} = f(t,y)` with a non-identity mass matrix, a
second ``SUNLinearSolver`` object for the mass matrix must be attached with
the following function (see :numref:`FIRKODE.Mathematics.mass`). The user
must also supply either a routine that fills the mass matrix
(:c:func:`FIRKodeSetMassFn`, for matrix-based linear solvers) or a routine
that computes mass matrix-vector products (:c:func:`FIRKodeSetMassTimes`,
for matrix-free linear solvers). Mass matrices are never approximated by
difference quotients.

.. c:function:: int FIRKodeSetMassLinearSolver(void* firkode_mem, SUNLinearSolver LS, SUNMatrix M, sunbooleantype time_dep)

   The function :c:func:`FIRKodeSetMassLinearSolver` attaches a generic
   ``SUNLinearSolver`` object ``LS`` and corresponding template mass
   ``SUNMatrix`` object ``M`` (if applicable) to FIRKODE for the solution of
   linear systems with the mass matrix.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``LS`` -- ``SUNLinearSolver`` object to use for solving linear systems with the mass matrix.
      - ``M`` -- ``SUNMatrix`` object used as a template for the mass matrix (or ``NULL`` for a matrix-free linear solver).
      - ``time_dep`` -- flag denoting whether the mass matrix depends on the independent variable. Only ``SUNFALSE`` is supported.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The call was successful.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_ILL_INPUT`` -- The FIRKLS interface is not compatible with the ``LS`` or ``M`` input objects or with the current ``N_Vector`` module, ``time_dep`` is ``SUNTRUE``, or :c:func:`FIRKodeInit` has not been called.
      - ``FIRKLS_SUNLS_FAIL`` -- A call to the ``LS`` object failed.
      - ``FIRKLS_MEM_FAIL`` -- A memory allocation request failed.

   **Notes:**
      The mass matrix linear solver and the Jacobian linear solver attached
      with :c:func:`FIRKodeSetLinearSolver` must either both be matrix-based
      (in which case ``M`` and the Jacobian template ``A`` must have the same
      ``SUNMatrix`` type, since the block system matrix :math:`M - \gamma J`
      is assembled from both) or both be matrix-free; mixed configurations
      are rejected at the first call to :c:func:`FIRKodeEvolve`.

      If ``LS`` is a matrix-based linear solver, then the template mass
      matrix ``M`` will be used in the solve process, so if additional storage
      is required within the ``SUNMatrix`` object (e.g. for factorization of a
      banded matrix), ensure that the input object is allocated with
      sufficient size. For direct linear solvers FIRKLS keeps a separate copy
      of ``M`` for the factorization so that ``M`` itself remains available
      for matrix-vector products.

      The mass matrix is constructed and factored once, at the first call to
      :c:func:`FIRKodeEvolve` after :c:func:`FIRKodeInit` or
      :c:func:`FIRKodeReInit`, and reused for the rest of the integration.


.. _FIRKODE.Usage.rootinit:

Rootfinding initialization function
-----------------------------------

While solving the IVP, FIRKODE has the capability to find the roots of a set
of user-defined functions. To activate the root finding algorithm, call the
following function. This is normally called only once, prior to the first
call to :c:func:`FIRKodeEvolve`, but if the rootfinding problem is to be
changed during the solution, :c:func:`FIRKodeRootInit` can also be called
prior to a continuation call to :c:func:`FIRKodeEvolve`.

.. c:function:: int FIRKodeRootInit(void* firkode_mem, int nrtfn, FIRKRootFn g)

   The function :c:func:`FIRKodeRootInit` specifies that the roots of a set
   of functions :math:`g_i(t,y)` are to be found while the IVP is being
   solved.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nrtfn`` -- is the number of root functions :math:`g_i`.
      - ``g`` -- is the C function which defines the ``nrtfn`` functions :math:`g_i(t,y)` whose roots are sought. See :c:type:`FIRKRootFn` for details.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The call was successful.
      - ``FIRK_MEM_NULL`` -- The ``firkode_mem`` argument was ``NULL``.
      - ``FIRK_MEM_FAIL`` -- A memory allocation failed.
      - ``FIRK_ILL_INPUT`` -- The function ``g`` is ``NULL``, but ``nrtfn`` :math:`> 0`.

   **Notes:**
      If a new IVP is to be solved with a call to :c:func:`FIRKodeReInit`,
      where the new IVP has no rootfinding problem but the prior one did,
      then call :c:func:`FIRKodeRootInit` with ``nrtfn = 0``.


.. _FIRKODE.Usage.evolve:

FIRKODE solver function
-----------------------

This is the central step in the solution process -- the call to perform the
integration of the IVP. One of the input arguments (``itask``) specifies one
of two modes as to where FIRKODE is to return a solution. But these modes
are modified if the user has set a stop time (with
:c:func:`FIRKodeSetStopTime`) or requested rootfinding.

.. c:function:: int FIRKodeEvolve(void* firkode_mem, sunrealtype tout, N_Vector yout, sunrealtype* tret, int itask)

   The function :c:func:`FIRKodeEvolve` integrates the ODE over an interval in :math:`t`.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``tout`` -- the next time at which a computed solution is desired.
      - ``yout`` -- the computed solution vector.
      - ``tret`` -- the time reached by the solver (output).
      - ``itask`` -- a flag indicating the job of the solver for the next user step. The ``FIRK_NORMAL`` option causes the solver to take internal steps until it has reached or just passed the user-specified ``tout`` parameter. The solver then interpolates in order to return an approximate value of :math:`y(t_{out})`. The ``FIRK_ONE_STEP`` option tells the solver to take just one internal step and then return the solution at the point reached by that step.

   **Return value:**
      - ``FIRK_SUCCESS`` -- :c:func:`FIRKodeEvolve` succeeded and no roots were found.
      - ``FIRK_TSTOP_RETURN`` -- :c:func:`FIRKodeEvolve` succeeded by reaching the stopping point specified through the optional input function :c:func:`FIRKodeSetStopTime`.
      - ``FIRK_ROOT_RETURN`` -- :c:func:`FIRKodeEvolve` succeeded and found one or more roots. In this case, ``tret`` is the location of the root. If ``nrtfn`` :math:`> 1`, call :c:func:`FIRKodeGetRootInfo` to see which :math:`g_i` were found to have a root.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized through a previous call to :c:func:`FIRKodeCreate`.
      - ``FIRK_NO_MALLOC`` -- The FIRKODE memory was not allocated by a call to :c:func:`FIRKodeInit`.
      - ``FIRK_ILL_INPUT`` -- One of the inputs to :c:func:`FIRKodeEvolve` was illegal, or some other input to the solver was illegal or missing. The latter category includes the following situations: the tolerances have not been set; no linear solver has been attached; a component of the error weight vector became zero during internal time-stepping; the direction of a fixed step size disagrees with ``tout``; a root of one of the root functions was found both at a point :math:`t` and also very near :math:`t`.
      - ``FIRK_TOO_CLOSE`` -- The initial time :math:`t_0` and the output time :math:`t_{out}` are too close to each other and the user did not specify an initial step size.
      - ``FIRK_TOO_MUCH_WORK`` -- The solver took ``mxstep`` internal steps but still could not reach ``tout``. The default value for ``mxstep`` is 500.
      - ``FIRK_TOO_MUCH_ACC`` -- The solver could not satisfy the accuracy demanded by the user for some internal step.
      - ``FIRK_ERR_FAILURE`` -- Either error test failures occurred too many times (default 7) during one internal time step, or with :math:`|h| = h_{min}`.
      - ``FIRK_CONV_FAILURE`` -- Either convergence test failures occurred too many times (default 10) during one internal time step, or with :math:`|h| = h_{min}`, or the Newton iteration failed in fixed-step mode.
      - ``FIRK_LINIT_FAIL`` -- The linear solver interface's initialization function failed.
      - ``FIRK_LSETUP_FAIL`` -- The linear solver interface's setup function failed in an unrecoverable manner.
      - ``FIRK_LSOLVE_FAIL`` -- The linear solver interface's solve function failed in an unrecoverable manner.
      - ``FIRK_RHSFUNC_FAIL`` -- The right-hand side function failed in an unrecoverable manner.
      - ``FIRK_FIRST_RHSFUNC_ERR`` -- The right-hand side function had a recoverable error at the first call.
      - ``FIRK_REPTD_RHSFUNC_ERR`` -- Convergence test failures occurred too many times due to repeated recoverable errors in the right-hand side function. This flag will also be returned if the right-hand side function had repeated recoverable errors during the estimation of an initial step size.
      - ``FIRK_UNREC_RHSFUNC_ERR`` -- The right-hand function had a recoverable error, but no recovery was possible. This occurs when the right-hand side function fails recoverably when evaluated at the beginning of a step.
      - ``FIRK_RTFUNC_FAIL`` -- The rootfinding function failed.
      - ``FIRK_NLS_INIT_FAIL`` -- The nonlinear solver initialization failed.
      - ``FIRK_NLS_SETUP_FAIL`` -- The nonlinear solver setup failed.
      - ``FIRK_NLS_FAIL`` -- The nonlinear solver failed in an unrecoverable manner.
      - ``FIRK_MASSINIT_FAIL`` -- The mass matrix solver interface's initialization function failed.
      - ``FIRK_MASSSETUP_FAIL`` -- The mass matrix solver interface's setup function failed.
      - ``FIRK_MASSSOLVE_FAIL`` -- The mass matrix solve failed in an unrecoverable manner.
      - ``FIRK_MASSMULT_FAIL`` -- The mass matrix-vector product failed.
      - ``FIRK_TABLE_FAIL`` -- The coefficient table of the selected method could not be constructed.
      - ``FIRK_CONTROLLER_ERR`` -- A call to the time step adaptivity controller failed.
      - ``FIRK_VECTOROP_ERR`` -- A vector operation failed.

   **Notes:**
      The vector ``yout`` can occupy the same space as the vector ``y0`` of
      initial conditions that was passed to :c:func:`FIRKodeInit`.

      In the ``FIRK_ONE_STEP`` mode, ``tout`` is used only on the first call,
      and only to get the direction and a rough scale of the independent
      variable.

      If a stop time is enabled (through a call to
      :c:func:`FIRKodeSetStopTime`), then :c:func:`FIRKodeEvolve` returns the
      solution at ``tstop``. Once the integrator returns at a stop time, any
      future testing for ``tstop`` is disabled (and can be re-enabled only
      through a new call to :c:func:`FIRKodeSetStopTime`).

      All failure return values are negative and so the test ``flag < 0``
      will trap all :c:func:`FIRKodeEvolve` failures.

      On any error return in which one or more internal steps were taken by
      :c:func:`FIRKodeEvolve`, the returned values of ``tret`` and ``yout``
      correspond to the farthest point reached in the integration. On all
      other error returns, ``tret`` and ``yout`` are left unchanged from the
      previous :c:func:`FIRKodeEvolve` return.


.. _FIRKODE.Usage.optional_input:

Optional input functions
------------------------

There are numerous optional input parameters that control the behavior of
the FIRKODE solver. FIRKODE provides functions that can be used to change
these optional input parameters from their default values. The main inputs
are divided into the following categories:

* :numref:`FIRKODE.Usage.optional_input.main_table` lists the main FIRKODE
  optional input functions,

* :numref:`FIRKODE.Usage.optional_input.adapt_table` lists the time step
  adaptivity optional input functions,

* :numref:`FIRKODE.Usage.optional_input.nls_table` lists the nonlinear solver
  optional input functions,

* :numref:`FIRKODE.Usage.optional_input.ls_table` lists the FIRKLS linear
  solver interface optional input functions,

* :numref:`FIRKODE.Usage.optional_input.mass_table` lists the mass matrix
  solver interface optional input functions, and

* :numref:`FIRKODE.Usage.optional_input.root_table` lists the rootfinding
  optional input functions.

These optional inputs are described in detail in the remainder of this
section. For the most casual use of FIRKODE, the reader can skip to
:numref:`FIRKODE.Usage.user_supplied`.

We note that, on an error return, all of the optional input functions send an
error message to the error handler function. All error return values are
negative, so the test ``flag < 0`` will catch all errors.

The optional input calls can, unless otherwise noted, be executed in any
order. A call to a ``FIRKodeSet***`` function can, unless otherwise noted, be
made at any time from the user's calling program and, if successful, takes
effect immediately.

.. _FIRKODE.Usage.optional_input.main:

Main solver optional input functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _FIRKODE.Usage.optional_input.main_table:

.. table:: Optional inputs for FIRKODE

   +-------------------------------+----------------------------------------------+----------------+
   |      **Optional input**       |              **Function name**               |  **Default**   |
   +===============================+==============================================+================+
   | User data                     | :c:func:`FIRKodeSetUserData`                 | ``NULL``       |
   +-------------------------------+----------------------------------------------+----------------+
   | Maximum no. of internal steps | :c:func:`FIRKodeSetMaxNumSteps`              | 500            |
   | before :math:`t_{out}`        |                                              |                |
   +-------------------------------+----------------------------------------------+----------------+
   | Maximum no. of warnings for   | :c:func:`FIRKodeSetMaxHnilWarns`             | 10             |
   | :math:`t_n+h=t_n`             |                                              |                |
   +-------------------------------+----------------------------------------------+----------------+
   | Initial step size             | :c:func:`FIRKodeSetInitStep`                 | estimated      |
   +-------------------------------+----------------------------------------------+----------------+
   | Minimum absolute step size    | :c:func:`FIRKodeSetMinStep`                  | 0.0            |
   +-------------------------------+----------------------------------------------+----------------+
   | Maximum absolute step size    | :c:func:`FIRKodeSetMaxStep`                  | :math:`\infty` |
   +-------------------------------+----------------------------------------------+----------------+
   | Fixed step size               | :c:func:`FIRKodeSetFixedStep`                | disabled       |
   +-------------------------------+----------------------------------------------+----------------+
   | Value of :math:`t_{stop}`     | :c:func:`FIRKodeSetStopTime`                 | undefined      |
   +-------------------------------+----------------------------------------------+----------------+
   | Interpolate at                | :c:func:`FIRKodeSetInterpolateStopTime`      | ``SUNFALSE``   |
   | :math:`t_{stop}`              |                                              |                |
   +-------------------------------+----------------------------------------------+----------------+
   | Disable the stop time         | :c:func:`FIRKodeClearStopTime`               | N/A            |
   +-------------------------------+----------------------------------------------+----------------+
   | Maximum no. of error test     | :c:func:`FIRKodeSetMaxErrTestFails`          | 7              |
   | failures                      |                                              |                |
   +-------------------------------+----------------------------------------------+----------------+
   | Maximum no. of convergence    | :c:func:`FIRKodeSetMaxConvFails`             | 10             |
   | failures                      |                                              |                |
   +-------------------------------+----------------------------------------------+----------------+
   | Linear solver setup frequency | :c:func:`FIRKodeSetLSetupFrequency`          | 20             |
   +-------------------------------+----------------------------------------------+----------------+
   | Max change in :math:`\gamma`  | :c:func:`FIRKodeSetDeltaGammaMaxLSetup`      | 0.2            |
   | without a linear solver setup |                                              |                |
   +-------------------------------+----------------------------------------------+----------------+
   | Stage predictor               | :c:func:`FIRKodeSetPredictorMethod`          | extrapolation  |
   +-------------------------------+----------------------------------------------+----------------+
   | Second filtering of the error | :c:func:`FIRKodeSetErrorRefilter`            | ``SUNTRUE``    |
   | estimate                      |                                              |                |
   +-------------------------------+----------------------------------------------+----------------+

.. c:function:: int FIRKodeSetUserData(void* firkode_mem, void* user_data)

   The function :c:func:`FIRKodeSetUserData` specifies the user data block
   ``user_data`` and attaches it to the main FIRKODE memory block.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``user_data`` -- pointer to the user data.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

   **Notes:**
      If specified, the pointer to ``user_data`` is passed to all
      user-supplied functions that have it as an argument. Otherwise, a
      ``NULL`` pointer is passed. If ``user_data`` is needed in user linear
      solver or preconditioner functions, the call to
      :c:func:`FIRKodeSetUserData` must be made *before* the call to specify
      the linear solver.

.. c:function:: int FIRKodeSetMaxNumSteps(void* firkode_mem, long int mxsteps)

   The function :c:func:`FIRKodeSetMaxNumSteps` specifies the maximum number
   of steps to be taken by the solver in its attempt to reach the next
   output time.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``mxsteps`` -- maximum allowed number of steps.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

   **Notes:**
      Passing ``mxsteps = 0`` results in FIRKODE using the default value
      (500). Passing ``mxsteps < 0`` disables the test (not recommended).

.. c:function:: int FIRKodeSetMaxHnilWarns(void* firkode_mem, int mxhnil)

   The function :c:func:`FIRKodeSetMaxHnilWarns` specifies the maximum number
   of messages issued by the solver warning that :math:`t + h = t` on the
   next internal step.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``mxhnil`` -- maximum number of warning messages (:math:`> 0`).

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

   **Notes:**
      The default value is 10. A negative value for ``mxhnil`` indicates that
      no warning messages should be issued.

.. c:function:: int FIRKodeSetInitStep(void* firkode_mem, sunrealtype hin)

   The function :c:func:`FIRKodeSetInitStep` specifies the initial step size.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``hin`` -- value of the initial step size to be attempted. Pass 0.0 to have FIRKODE use the default value.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.
      - ``FIRK_ILL_INPUT`` -- A fixed step size is in use.

   **Notes:**
      By default, FIRKODE estimates the initial step size to be the solution
      :math:`h` of the equation :math:`\| h^2 \ddot{y} / 2 \|_{\text{WRMS}} =
      1`, where :math:`\ddot{y}` is an estimated value of the second
      derivative of the solution at :math:`t_0` (see
      :numref:`FIRKODE.Mathematics.initial_step`).

.. c:function:: int FIRKodeSetMinStep(void* firkode_mem, sunrealtype hmin)

   The function :c:func:`FIRKodeSetMinStep` specifies a lower bound on the
   magnitude of the step size.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``hmin`` -- minimum absolute value of the step size (:math:`\ge 0.0`).

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.
      - ``FIRK_ILL_INPUT`` -- Either ``hmin`` is negative or it exceeds the maximum allowable step size.

   **Notes:**
      The default value is 0.0. Pass ``hmin = 0.0`` to obtain the default.

.. c:function:: int FIRKodeSetMaxStep(void* firkode_mem, sunrealtype hmax)

   The function :c:func:`FIRKodeSetMaxStep` specifies an upper bound on the
   magnitude of the step size.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``hmax`` -- maximum absolute value of the step size (:math:`\ge 0.0`).

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.
      - ``FIRK_ILL_INPUT`` -- Either ``hmax`` is negative or it is smaller than the minimum allowable step size.

   **Notes:**
      Pass ``hmax = 0.0`` to obtain the default value :math:`\infty`.

.. c:function:: int FIRKodeSetFixedStep(void* firkode_mem, sunrealtype hfixed)

   The function :c:func:`FIRKodeSetFixedStep` disables time step adaptivity
   and specifies the fixed step size to be used for all internal steps.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``hfixed`` -- the fixed step size (with sign indicating the direction of integration). Pass 0.0 to re-enable time step adaptivity.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

   **Notes:**
      In fixed-step mode the local error estimate is not computed, the error
      test is bypassed, and a failure of the Newton iteration immediately
      terminates the integration with ``FIRK_CONV_FAILURE``. The minimum and
      maximum step sizes are set to ``hfixed``; when adaptivity is re-enabled
      these bounds revert to their defaults. Fixed-step mode is intended for
      testing and for problems where the user has independent knowledge of a
      suitable step size.

.. c:function:: int FIRKodeSetStopTime(void* firkode_mem, sunrealtype tstop)

   The function :c:func:`FIRKodeSetStopTime` specifies the value of the
   independent variable :math:`t` past which the solution is not to proceed.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``tstop`` -- value of the independent variable past which the solution is not to proceed.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.
      - ``FIRK_ILL_INPUT`` -- The value of ``tstop`` is behind the current :math:`t` value, in the direction of integration.

   **Notes:**
      The default, if this routine is not called, is that no stop time is
      imposed. Once the integrator returns at a stop time, any future testing
      for ``tstop`` is disabled (and can be re-enabled only through a new call
      to :c:func:`FIRKodeSetStopTime`). A stop time not reached before a call
      to :c:func:`FIRKodeReInit` or :c:func:`FIRKodeReset` will remain active
      but can be disabled by calling :c:func:`FIRKodeClearStopTime`.

.. c:function:: int FIRKodeSetInterpolateStopTime(void* firkode_mem, sunbooleantype interp)

   The function :c:func:`FIRKodeSetInterpolateStopTime` specifies that the
   output solution should be interpolated when the current :math:`t` equals
   the specified ``tstop`` (instead of merely copying the internal solution
   :math:`y_n`).

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``interp`` -- flag indicating to use interpolation (``SUNTRUE``) or copy (``SUNFALSE``).

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

.. c:function:: int FIRKodeClearStopTime(void* firkode_mem)

   Disables the stop time set with :c:func:`FIRKodeSetStopTime`.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The stop time has been disabled.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

   **Notes:**
      The stop time can be re-enabled though a new call to
      :c:func:`FIRKodeSetStopTime`.

.. c:function:: int FIRKodeSetMaxErrTestFails(void* firkode_mem, int maxnef)

   The function :c:func:`FIRKodeSetMaxErrTestFails` specifies the maximum
   number of error test failures permitted in attempting one step.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``maxnef`` -- maximum number of error test failures allowed on one step (:math:`> 0`).

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

   **Notes:**
      The default value is 7. A non-positive value restores the default.

.. c:function:: int FIRKodeSetMaxConvFails(void* firkode_mem, int maxncf)

   The function :c:func:`FIRKodeSetMaxConvFails` specifies the maximum number
   of nonlinear solver convergence failures permitted during one step.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``maxncf`` -- maximum number of allowable nonlinear solver convergence failures per step (:math:`> 0`).

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

   **Notes:**
      The default value is 10. A non-positive value restores the default.

.. c:function:: int FIRKodeSetLSetupFrequency(void* firkode_mem, long int msbp)

   The function :c:func:`FIRKodeSetLSetupFrequency` specifies the frequency
   of calls to the linear solver setup function (which updates the matrix
   :math:`M - \gamma J` and its factorization or preconditioner).

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``msbp`` -- the number of steps between linear solver setup calls (:math:`\ge 0`).

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.
      - ``FIRK_ILL_INPUT`` -- The value of ``msbp`` is negative.

   **Notes:**
      The default value is 20; passing 0 restores the default. The linear
      solver setup is additionally called on the first step, after a
      nonlinear solver or error test failure, and whenever :math:`\gamma`
      has changed by more than the fraction set with
      :c:func:`FIRKodeSetDeltaGammaMaxLSetup` (see
      :numref:`FIRKODE.Mathematics.linear`).

.. c:function:: int FIRKodeSetDeltaGammaMaxLSetup(void* firkode_mem, sunrealtype dgmax_lsetup)

   The function :c:func:`FIRKodeSetDeltaGammaMaxLSetup` specifies the maximum
   relative change in :math:`\gamma` since the last linear solver setup that
   does not trigger a new setup.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``dgmax_lsetup`` -- the relative change threshold (:math:`\ge 0`).

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

   **Notes:**
      The default value is 0.2. A negative value restores the default.

.. c:function:: int FIRKodeSetPredictorMethod(void* firkode_mem, int method)

   The function :c:func:`FIRKodeSetPredictorMethod` selects the method used
   to compute the initial guess for the stage increments of each step.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``method`` -- ``FIRK_PREDICT_EXTRAPOLATE`` (extrapolate the collocation polynomial of the previous step, the default) or ``FIRK_PREDICT_TRIVIAL`` (start from the current solution, i.e., zero increments).

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.
      - ``FIRK_ILL_INPUT`` -- The value of ``method`` is illegal.

   **Notes:**
      The trivial predictor is always used on the first step and after a
      reinitialization. The trivial predictor can be more robust when the
      solution is strongly non-smooth, at the price of more Newton iterations.

.. c:function:: int FIRKodeSetErrorRefilter(void* firkode_mem, sunbooleantype onoff)

   The function :c:func:`FIRKodeSetErrorRefilter` enables or disables the
   second filtering of the local error estimate on the first step and after
   a rejected step (see :numref:`FIRKODE.Mathematics.error`).

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``onoff`` -- ``SUNTRUE`` to enable (the default) or ``SUNFALSE`` to disable the second filtering.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.


.. _FIRKODE.Usage.optional_input.adapt:

Time step adaptivity optional input functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following functions control the selection of the step size after each
step attempt (see :numref:`FIRKODE.Mathematics.adaptivity`). None of them
has any effect in fixed-step mode.

.. _FIRKODE.Usage.optional_input.adapt_table:

.. table:: Optional time step adaptivity inputs for FIRKODE

   +-------------------------------+----------------------------------------------+----------------+
   |      **Optional input**       |              **Function name**               |  **Default**   |
   +===============================+==============================================+================+
   | Time step controller object   | :c:func:`FIRKodeSetAdaptController`          | ImpGus         |
   +-------------------------------+----------------------------------------------+----------------+
   | Safety factor                 | :c:func:`FIRKodeSetSafetyFactor`             | 0.9            |
   +-------------------------------+----------------------------------------------+----------------+
   | Error bias                    | :c:func:`FIRKodeSetErrorBias`                | controller     |
   +-------------------------------+----------------------------------------------+----------------+
   | Maximum step growth           | :c:func:`FIRKodeSetMaxGrowth`                | 8.0            |
   +-------------------------------+----------------------------------------------+----------------+
   | Minimum step reduction        | :c:func:`FIRKodeSetMinReduction`             | 0.2            |
   +-------------------------------+----------------------------------------------+----------------+
   | Step size dead band           | :c:func:`FIRKodeSetFixedStepBounds`          | [1.0, 1.2]     |
   +-------------------------------+----------------------------------------------+----------------+
   | Maximum first step growth     | :c:func:`FIRKodeSetMaxFirstGrowth`           | 10000.0        |
   +-------------------------------+----------------------------------------------+----------------+
   | Maximum growth after error    | :c:func:`FIRKodeSetMaxEFailGrowth`           | 0.3            |
   | test failures                 |                                              |                |
   +-------------------------------+----------------------------------------------+----------------+
   | No. of error test failures    | :c:func:`FIRKodeSetSmallNumEFails`           | 2              |
   | before limiting growth        |                                              |                |
   +-------------------------------+----------------------------------------------+----------------+
   | Step reduction after a        | :c:func:`FIRKodeSetMaxCFailGrowth`           | 0.25           |
   | convergence failure           |                                              |                |
   +-------------------------------+----------------------------------------------+----------------+
   | Newton iteration count        | :c:func:`FIRKodeSetNewtonCountSafety`        | ``SUNTRUE``    |
   | safety adjustment             |                                              |                |
   +-------------------------------+----------------------------------------------+----------------+

.. c:function:: int FIRKodeSetAdaptController(void* firkode_mem, SUNAdaptController C)

   Sets a user-supplied time step controller object.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``C`` -- user-supplied time adaptivity controller. If ``C`` is ``NULL`` then the default controller is created.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.
      - ``FIRK_MEM_FAIL`` -- ``C`` was ``NULL`` and the default controller could not be allocated.
      - ``FIRK_ILL_INPUT`` -- The controller is not of type ``SUN_ADAPTCONTROLLER_H``.

   **Notes:**
      The default controller is the implicit Gustafsson controller
      :c:func:`SUNAdaptController_ImpGus` (see
      :numref:`SUNAdaptController.ImExGus`), which is the predictive
      controller of RADAU5; the I controller of
      :numref:`SUNAdaptController.Soderlind` reproduces the classical
      (non-predictive) RADAU5 step size selection. The controller is
      supplied with the order :math:`p = s` of the error estimate. The
      user-supplied controller is not freed by FIRKODE.

.. c:function:: int FIRKodeSetSafetyFactor(void* firkode_mem, sunrealtype safety)

   Specifies the safety factor to be applied to the step size proposed by the
   time step controller.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``safety`` -- safety factor (:math:`0 < \mathrm{safety} \le 1`).

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.
      - ``FIRK_ILL_INPUT`` -- The safety factor exceeds 1.

   **Notes:**
      The default value is 0.9; a non-positive value restores the default.

.. c:function:: int FIRKodeSetErrorBias(void* firkode_mem, sunrealtype bias)

   Specifies the bias to be applied to the error estimates within the time
   step controller.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``bias`` -- bias applied to the error estimate (:math:`\ge 1`); any value below 1 restores the controller's default.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized, or no controller is attached.
      - ``FIRK_CONTROLLER_ERR`` -- The call to the controller failed.

   **Notes:**
      This function must be called *after* :c:func:`FIRKodeSetAdaptController`
      since the bias is stored in the controller object.

.. c:function:: int FIRKodeSetMaxGrowth(void* firkode_mem, sunrealtype mx_growth)

   Specifies the maximum growth factor of the step size between two steps.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``mx_growth`` -- maximum allowed growth factor (:math:`> 1`).

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

   **Notes:**
      The default value is 8.0; any value :math:`\le 1` restores the default.

.. c:function:: int FIRKodeSetMinReduction(void* firkode_mem, sunrealtype eta_min)

   Specifies the minimum reduction factor of the step size between two steps.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``eta_min`` -- minimum allowed reduction factor (:math:`0 < \eta_{\min} < 1`).

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

   **Notes:**
      The default value is 0.2; any value outside :math:`(0, 1)` restores
      the default.

.. c:function:: int FIRKodeSetFixedStepBounds(void* firkode_mem, sunrealtype lb, sunrealtype ub)

   Specifies the interval of step size ratios :math:`h_{\text{new}}/h` within
   which the step size is left unchanged after a successful step, so that
   the matrix :math:`M - \gamma J` need not be rebuilt.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``lb`` -- lower bound of the interval (:math:`\le 1`).
      - ``ub`` -- upper bound of the interval (:math:`\ge 1`).

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

   **Notes:**
      The default interval is :math:`[1.0, 1.2]`; an interval that does not
      contain 1 restores the default.

.. c:function:: int FIRKodeSetMaxFirstGrowth(void* firkode_mem, sunrealtype etamx1)

   Specifies the maximum growth factor of the step size after the first
   step.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``etamx1`` -- maximum allowed growth factor after the first step (:math:`> 1`).

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

   **Notes:**
      The default value is 10000.0; any value :math:`\le 1` restores the
      default.

.. c:function:: int FIRKodeSetMaxEFailGrowth(void* firkode_mem, sunrealtype etamxf)

   Specifies the maximum step size ratio to be used after a number of
   consecutive error test failures.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``etamxf`` -- maximum step size ratio after repeated error test failures (:math:`0 < \eta_{\max,f} \le 1`).

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

   **Notes:**
      The default value is 0.3; any value outside :math:`(0, 1]` restores
      the default. The limit takes effect after the number of failures set
      with :c:func:`FIRKodeSetSmallNumEFails`.

.. c:function:: int FIRKodeSetSmallNumEFails(void* firkode_mem, int small_nef)

   Specifies the number of consecutive error test failures after which the
   step size ratio is limited by the value set with
   :c:func:`FIRKodeSetMaxEFailGrowth`.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``small_nef`` -- number of error test failures (:math:`> 0`).

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

   **Notes:**
      The default value is 2; a non-positive value restores the default.

.. c:function:: int FIRKodeSetMaxCFailGrowth(void* firkode_mem, sunrealtype etacf)

   Specifies the step size reduction factor applied after a nonlinear solver
   convergence failure (unless the failure was predicted by the convergence
   test, in which case the reduction is computed from the observed
   convergence rate, see :numref:`FIRKODE.Mathematics.Newton`).

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``etacf`` -- step size reduction factor (:math:`0 < \eta_{cf} \le 1`).

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

   **Notes:**
      The default value is 0.25; any value outside :math:`(0, 1]` restores
      the default.

.. c:function:: int FIRKodeSetNewtonCountSafety(void* firkode_mem, sunbooleantype onoff)

   Enables or disables the reduction of the safety factor according to the
   number of Newton iterations of the last nonlinear solve (see
   :numref:`FIRKODE.Mathematics.adaptivity`).

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``onoff`` -- ``SUNTRUE`` to enable (the default) or ``SUNFALSE`` to disable the adjustment.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.


.. _FIRKODE.Usage.optional_input.nls:

Nonlinear solver optional input functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following functions control the simplified Newton iteration used to
solve the stage equations (see :numref:`FIRKODE.Mathematics.Newton`).

.. _FIRKODE.Usage.optional_input.nls_table:

.. table:: Optional nonlinear solver inputs for FIRKODE

   +-------------------------------+----------------------------------------------+----------------+
   |      **Optional input**       |              **Function name**               |  **Default**   |
   +===============================+==============================================+================+
   | Maximum no. of Newton         | :c:func:`FIRKodeSetMaxNonlinIters`           | 7              |
   | iterations                    |                                              |                |
   +-------------------------------+----------------------------------------------+----------------+
   | Nonlinear convergence         | :c:func:`FIRKodeSetNonlinConvCoef`           | 0.1            |
   | coefficient                   |                                              |                |
   +-------------------------------+----------------------------------------------+----------------+
   | Newton divergence rate        | :c:func:`FIRKodeSetNonlinDivergenceRate`     | 0.99           |
   +-------------------------------+----------------------------------------------+----------------+
   | Convergence rate triggering a | :c:func:`FIRKodeSetJacBadConvRate`           | 1.0 (disabled) |
   | Jacobian update               |                                              |                |
   +-------------------------------+----------------------------------------------+----------------+

.. c:function:: int FIRKodeSetMaxNonlinIters(void* firkode_mem, int maxcor)

   The function :c:func:`FIRKodeSetMaxNonlinIters` specifies the maximum
   number of Newton iterations permitted per step.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``maxcor`` -- maximum number of Newton iterations per step (:math:`> 0`).

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.
      - ``FIRK_ILL_INPUT`` -- The nonlinear solver rejected the value.

   **Notes:**
      The default value is 7; a non-positive value restores the default.

.. c:function:: int FIRKodeSetNonlinConvCoef(void* firkode_mem, sunrealtype nlscoef)

   The function :c:func:`FIRKodeSetNonlinConvCoef` specifies the safety
   factor :math:`\epsilon_{nls}` used in the nonlinear convergence test
   :eq:`FIRKODE_nlstest`.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nlscoef`` -- coefficient in the nonlinear convergence test (:math:`> 0`).

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

   **Notes:**
      The default value is 0.1; a non-positive value restores the default.
      The value also scales the tolerances of iterative linear solvers (see
      :c:func:`FIRKodeSetEpsLin` and :c:func:`FIRKodeSetStageEpsLin`).

.. c:function:: int FIRKodeSetNonlinDivergenceRate(void* firkode_mem, sunrealtype theta_max)

   The function :c:func:`FIRKodeSetNonlinDivergenceRate` specifies the
   contraction rate above which the Newton iteration is declared divergent.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``theta_max`` -- the divergence threshold (:math:`> 0`).

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

   **Notes:**
      The default value is 0.99; a non-positive value restores the default.

.. c:function:: int FIRKodeSetJacBadConvRate(void* firkode_mem, sunrealtype theta_jbad)

   The function :c:func:`FIRKodeSetJacBadConvRate` specifies a contraction
   rate of the Newton iteration above which the Jacobian is re-evaluated at
   the next step even though the iteration converged.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``theta_jbad`` -- the contraction rate threshold; values :math:`\ge 1` disable the test.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

   **Notes:**
      The test is disabled by default (the value 1.0); a non-positive value
      restores the default. A threshold of, e.g., 0.1 mimics the strategy of
      RADAU5 for problems with rapidly varying Jacobians.


.. _FIRKODE.Usage.optional_input.ls:

Linear solver interface optional input functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The mathematical explanation of the linear solver methods available to
FIRKODE is provided in :numref:`FIRKODE.Mathematics.linear`. We group the
user-callable routines into four categories: general routines concerning the
update frequency of the Jacobian, routines for matrix-based linear solvers,
routines for matrix-free linear solvers, and routines for the internal
Krylov iteration on the coupled stage system. All of these functions must be
called *after* the FIRKLS linear solver interface has been initialized with
:c:func:`FIRKodeSetLinearSolver`.

.. _FIRKODE.Usage.optional_input.ls_table:

.. table:: Optional inputs for the FIRKLS linear solver interface

   +-------------------------------------+----------------------------------------------+----------------+
   |      **Optional input**             |              **Function name**               |  **Default**   |
   +=====================================+==============================================+================+
   | Jacobian function                   | :c:func:`FIRKodeSetJacFn`                    | DQ             |
   +-------------------------------------+----------------------------------------------+----------------+
   | Linear system function              | :c:func:`FIRKodeSetLinSysFn`                 | internal       |
   +-------------------------------------+----------------------------------------------+----------------+
   | Jacobian evaluation frequency       | :c:func:`FIRKodeSetJacEvalFrequency`         | 51             |
   +-------------------------------------+----------------------------------------------+----------------+
   | Max change in :math:`\gamma`        | :c:func:`FIRKodeSetDeltaGammaMaxBadJac`      | 0.2            |
   | attributed to a bad Jacobian        |                                              |                |
   +-------------------------------------+----------------------------------------------+----------------+
   | Block linear solve tolerance factor | :c:func:`FIRKodeSetEpsLin`                   | 0.05           |
   +-------------------------------------+----------------------------------------------+----------------+
   | Newton linear solve tolerance       | :c:func:`FIRKodeSetLSNormFactor`             | vector length  |
   | conversion factor                   |                                              |                |
   +-------------------------------------+----------------------------------------------+----------------+
   | Preconditioner functions            | :c:func:`FIRKodeSetPreconditioner`           | ``NULL``       |
   +-------------------------------------+----------------------------------------------+----------------+
   | Jacobian-times-vector functions     | :c:func:`FIRKodeSetJacTimes`                 | DQ             |
   +-------------------------------------+----------------------------------------------+----------------+
   | Jacobian-times-vector DQ RHS        | :c:func:`FIRKodeSetJacTimesRhsFn`            | ``NULL``       |
   | function                            |                                              |                |
   +-------------------------------------+----------------------------------------------+----------------+
   | Stage Krylov subspace dimension     | :c:func:`FIRKodeSetStageSolverMaxl`          | max(5, 3s)     |
   +-------------------------------------+----------------------------------------------+----------------+
   | Stage Krylov restarts               | :c:func:`FIRKodeSetStageSolverMaxRestarts`   | 1              |
   +-------------------------------------+----------------------------------------------+----------------+
   | Stage Krylov tolerance factor       | :c:func:`FIRKodeSetStageEpsLin`              | 0.05           |
   +-------------------------------------+----------------------------------------------+----------------+

.. c:function:: int FIRKodeSetJacFn(void* firkode_mem, FIRKLsJacFn jac)

   The function :c:func:`FIRKodeSetJacFn` specifies the Jacobian
   approximation routine to be used for a matrix-based solver within the
   FIRKLS interface.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``jac`` -- user-defined Jacobian approximation function; pass ``NULL`` to use the internal difference-quotient approximation.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.
      - ``FIRKLS_ILL_INPUT`` -- The linear solver is not matrix-based, or the internal difference-quotient approximation is not available for the matrix type in use.

   **Notes:**
      By default, FIRKLS uses an internal difference quotient function for
      the :ref:`SUNMATRIX_DENSE <SUNMatrix.Dense>` and :ref:`SUNMATRIX_BAND
      <SUNMatrix.Band>` modules. If ``NULL`` is passed to ``jac``, this
      default function is used. An error will occur if no ``jac`` is supplied
      when using other matrix types. The function type :c:type:`FIRKLsJacFn`
      is described in :numref:`FIRKODE.Usage.user_supplied`.

.. c:function:: int FIRKodeSetLinSysFn(void* firkode_mem, FIRKLsLinSysFn linsys)

   The function :c:func:`FIRKodeSetLinSysFn` specifies the linear system
   approximation routine to be used for a matrix-based solver within the
   FIRKLS interface.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``linsys`` -- user-defined linear system approximation function; pass ``NULL`` to use the internal routine.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.
      - ``FIRKLS_ILL_INPUT`` -- The linear solver is not matrix-based.

   **Notes:**
      By default, FIRKLS uses an internal linear system function that
      leverages the ``SUNMatrix`` API to form the system :math:`M - \gamma J`.
      If ``NULL`` is passed to ``linsys``, this default function is used. The
      function type :c:type:`FIRKLsLinSysFn` is described in
      :numref:`FIRKODE.Usage.user_supplied`.

.. c:function:: int FIRKodeSetJacEvalFrequency(void* firkode_mem, long int msbj)

   The function :c:func:`FIRKodeSetJacEvalFrequency` specifies the number of
   steps after which the Jacobian information is considered out-of-date.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``msbj`` -- the Jacobian re-computation or preconditioner update frequency (:math:`\ge 0`).

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.
      - ``FIRKLS_ILL_INPUT`` -- The value of ``msbj`` is negative.

   **Notes:**
      The default value is 51; passing 0 restores the default. If the
      Jacobian or preconditioner data is out-of-date when a linear solver
      setup is performed, it will be recomputed; otherwise the saved data is
      reused with the current :math:`\gamma` (see
      :numref:`FIRKODE.Mathematics.linear`).

.. c:function:: int FIRKodeSetDeltaGammaMaxBadJac(void* firkode_mem, sunrealtype dgmax_jbad)

   The function :c:func:`FIRKodeSetDeltaGammaMaxBadJac` specifies the
   relative change in :math:`\gamma` below which a nonlinear solver failure
   is attributed to an outdated Jacobian (which is then re-evaluated) rather
   than to the change in :math:`\gamma`.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``dgmax_jbad`` -- the relative change threshold (:math:`> 0`).

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

   **Notes:**
      The default value is 0.2; a non-positive value restores the default.

.. c:function:: int FIRKodeSetEpsLin(void* firkode_mem, sunrealtype eplifac)

   The function :c:func:`FIRKodeSetEpsLin` specifies the factor by which the
   nonlinear convergence coefficient is multiplied to obtain the residual
   tolerance of an iterative linear solver for the block systems (see
   :numref:`SUNLinSol.FIRKODE.Iterative.Tolerance`).

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``eplifac`` -- linear convergence safety factor (:math:`\ge 0.0`).

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

   **Notes:**
      The default value is 0.05; passing a non-positive value restores the
      default.

.. c:function:: int FIRKodeSetLSNormFactor(void* firkode_mem, sunrealtype nrmfac)

   The function :c:func:`FIRKodeSetLSNormFactor` specifies the factor to use
   when converting from the integrator tolerance (WRMS norm) to the linear
   solver tolerance (L2 norm) for the block systems.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nrmfac`` -- the norm conversion factor. If ``nrmfac`` is :math:`> 0` then the provided value is used; if ``nrmfac`` is :math:`= 0` then the conversion factor is computed as the square root of the vector length (the default); if ``nrmfac`` is :math:`< 0` then the conversion factor is computed using the vector dot product :math:`\sqrt{v^T v}` for a vector of ones.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

.. c:function:: int FIRKodeSetPreconditioner(void* firkode_mem, FIRKLsPrecSetupFn pset, FIRKLsPrecSolveFn psolve)

   The function :c:func:`FIRKodeSetPreconditioner` specifies the
   preconditioner setup and solve functions for an iterative linear solver.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``pset`` -- user-defined preconditioner setup function. Pass ``NULL`` if no setup is necessary.
      - ``psolve`` -- user-defined preconditioner solve function.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.
      - ``FIRKLS_SUNLS_FAIL`` -- An error occurred when setting up preconditioning in the ``SUNLinearSolver`` object used by the FIRKLS interface.
      - ``FIRKLS_ILL_INPUT`` -- The ``SUNLinearSolver`` object does not support user-supplied preconditioning.

   **Notes:**
      The default is ``NULL`` for both arguments (i.e., no preconditioning).
      The function types :c:type:`FIRKLsPrecSetupFn` and
      :c:type:`FIRKLsPrecSolveFn` are described in
      :numref:`FIRKODE.Usage.user_supplied`. When the user's linear solver is
      matrix-free, the preconditioner solve routine is also applied
      block-wise as the preconditioner of the internal Krylov iteration on
      the coupled stage system.

.. c:function:: int FIRKodeSetJacTimes(void* firkode_mem, FIRKLsJacTimesSetupFn jtsetup, FIRKLsJacTimesVecFn jtimes)

   The function :c:func:`FIRKodeSetJacTimes` specifies the Jacobian-vector
   setup and product functions.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``jtsetup`` -- user-defined Jacobian-vector setup function. Pass ``NULL`` if no setup is necessary.
      - ``jtimes`` -- user-defined Jacobian-vector product function; pass ``NULL`` to use the internal difference-quotient approximation.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

   **Notes:**
      The Jacobian-vector product is used by matrix-free linear solvers and,
      for methods with more than one stage, by the internal Krylov iteration
      on the coupled stage system (unless the saved Jacobian matrix of a
      matrix-based linear solver provides :c:func:`SUNMatMatvec`). The
      function types :c:type:`FIRKLsJacTimesSetupFn` and
      :c:type:`FIRKLsJacTimesVecFn` are described in
      :numref:`FIRKODE.Usage.user_supplied`.

.. c:function:: int FIRKodeSetJacTimesRhsFn(void* firkode_mem, FIRKRhsFn jtimesRhsFn)

   The function :c:func:`FIRKodeSetJacTimesRhsFn` specifies an alternative
   ODE right-hand side function for use in the internal Jacobian-vector
   product difference quotient approximation.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``jtimesRhsFn`` -- the C function which computes the alternative ODE right-hand side function to use in Jacobian-vector product difference quotient approximations; pass ``NULL`` to use the right-hand side function of the ODE.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.
      - ``FIRKLS_ILL_INPUT`` -- The internal difference quotient approximation is disabled (a user-supplied ``jtimes`` function is in use).

   **Notes:**
      This function is intended for problems where the right-hand side
      contains non-differentiable terms that should be excluded from the
      Jacobian approximation.

.. c:function:: int FIRKodeSetStageSolverMaxl(void* firkode_mem, int maxl)

   The function :c:func:`FIRKodeSetStageSolverMaxl` specifies the maximum
   dimension of the Krylov subspace of the internal flexible GMRES iteration
   on the coupled stage system.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``maxl`` -- the maximum Krylov subspace dimension (:math:`> 0`); a non-positive value restores the default.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

   **Notes:**
      The default is :math:`\max(5, 3s)` for an :math:`s`-stage method. The
      internal solver stores :math:`2(\mathrm{maxl}+1)` vectors of the size
      of the coupled stage system, i.e., :math:`2 s (\mathrm{maxl}+1)`
      vectors of the size of the ODE system, so this value controls the
      memory footprint of FIRKODE. It has no effect for the one-stage method.
      This function must be called before the first call to
      :c:func:`FIRKodeEvolve` or a call to :c:func:`FIRKodeReInit`.

.. c:function:: int FIRKodeSetStageSolverMaxRestarts(void* firkode_mem, int maxrs)

   The function :c:func:`FIRKodeSetStageSolverMaxRestarts` specifies the
   maximum number of restarts of the internal flexible GMRES iteration on the
   coupled stage system.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``maxrs`` -- the maximum number of restarts (:math:`\ge 0`); a negative value restores the default.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

   **Notes:**
      The default value is 1.

.. c:function:: int FIRKodeSetStageEpsLin(void* firkode_mem, sunrealtype eplifac)

   The function :c:func:`FIRKodeSetStageEpsLin` specifies the factor by which
   the nonlinear convergence coefficient is multiplied to obtain the residual
   tolerance of the internal flexible GMRES iteration on the coupled stage
   system.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``eplifac`` -- linear convergence safety factor (:math:`\ge 0.0`).

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

   **Notes:**
      The default value is 0.05; passing a non-positive value restores the
      default.


.. _FIRKODE.Usage.optional_input.mass:

Mass matrix solver interface optional input functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All of these functions must be called *after* the mass matrix solver
interface has been initialized with :c:func:`FIRKodeSetMassLinearSolver`.

.. _FIRKODE.Usage.optional_input.mass_table:

.. table:: Optional inputs for the FIRKLS mass matrix solver interface

   +-------------------------------------+----------------------------------------------+----------------+
   |      **Optional input**             |              **Function name**               |  **Default**   |
   +=====================================+==============================================+================+
   | Mass matrix function                | :c:func:`FIRKodeSetMassFn`                   | none           |
   +-------------------------------------+----------------------------------------------+----------------+
   | Mass matrix-times-vector functions  | :c:func:`FIRKodeSetMassTimes`                | none           |
   +-------------------------------------+----------------------------------------------+----------------+
   | Mass matrix preconditioner          | :c:func:`FIRKodeSetMassPreconditioner`       | ``NULL``       |
   | functions                           |                                              |                |
   +-------------------------------------+----------------------------------------------+----------------+
   | Mass matrix linear solve tolerance  | :c:func:`FIRKodeSetMassEpsLin`               | 0.05           |
   | factor                              |                                              |                |
   +-------------------------------------+----------------------------------------------+----------------+
   | Mass matrix linear solve tolerance  | :c:func:`FIRKodeSetMassLSNormFactor`         | vector length  |
   | conversion factor                   |                                              |                |
   +-------------------------------------+----------------------------------------------+----------------+

.. c:function:: int FIRKodeSetMassFn(void* firkode_mem, FIRKLsMassFn mass)

   The function :c:func:`FIRKodeSetMassFn` specifies the routine that fills
   the mass matrix for a matrix-based mass linear solver. This function is
   required when a non-``NULL`` ``SUNMatrix`` was passed to
   :c:func:`FIRKodeSetMassLinearSolver`.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``mass`` -- user-defined mass matrix function of type :c:type:`FIRKLsMassFn`.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_MASSMEM_NULL`` -- The mass matrix solver interface has not been initialized.
      - ``FIRKLS_ILL_INPUT`` -- The function ``mass`` is ``NULL`` or the mass linear solver is matrix-free.

.. c:function:: int FIRKodeSetMassTimes(void* firkode_mem, FIRKLsMassTimesSetupFn mtsetup, FIRKLsMassTimesVecFn mtimes, void* mtimes_data)

   The function :c:func:`FIRKodeSetMassTimes` specifies the mass
   matrix-vector product setup and multiply functions. This function is
   required when a ``NULL`` ``SUNMatrix`` was passed to
   :c:func:`FIRKodeSetMassLinearSolver`; it is optional otherwise (in which
   case the product is computed with :c:func:`SUNMatMatvec` by default).

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``mtsetup`` -- user-defined mass matrix-vector product setup function of type :c:type:`FIRKLsMassTimesSetupFn`, or ``NULL`` if no setup is needed.
      - ``mtimes`` -- user-defined mass matrix-vector product function of type :c:type:`FIRKLsMassTimesVecFn`.
      - ``mtimes_data`` -- a pointer to user data passed to ``mtsetup`` and ``mtimes``.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_MASSMEM_NULL`` -- The mass matrix solver interface has not been initialized.
      - ``FIRKLS_ILL_INPUT`` -- The function ``mtimes`` is ``NULL``.

   **Notes:**
      Since the mass matrix is constant, the setup function is called at most
      once.

.. c:function:: int FIRKodeSetMassPreconditioner(void* firkode_mem, FIRKLsMassPrecSetupFn psetup, FIRKLsMassPrecSolveFn psolve)

   The function :c:func:`FIRKodeSetMassPreconditioner` specifies the
   preconditioner setup and solve functions for an iterative mass matrix
   linear solver.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``psetup`` -- user-defined preconditioner setup function of type :c:type:`FIRKLsMassPrecSetupFn`, or ``NULL`` if no setup is needed.
      - ``psolve`` -- user-defined preconditioner solve function of type :c:type:`FIRKLsMassPrecSolveFn`.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_MASSMEM_NULL`` -- The mass matrix solver interface has not been initialized.
      - ``FIRKLS_SUNLS_FAIL`` -- An error occurred when setting up preconditioning in the ``SUNLinearSolver`` object.
      - ``FIRKLS_ILL_INPUT`` -- The ``SUNLinearSolver`` object does not support user-supplied preconditioning.

   **Notes:**
      Since the mass matrix is constant, the setup function is called at most
      once.

.. c:function:: int FIRKodeSetMassEpsLin(void* firkode_mem, sunrealtype eplifac)

   The function :c:func:`FIRKodeSetMassEpsLin` specifies the factor by which
   the nonlinear convergence coefficient is multiplied to obtain the residual
   tolerance of an iterative mass matrix linear solver.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``eplifac`` -- linear convergence safety factor (:math:`\ge 0.0`).

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_MASSMEM_NULL`` -- The mass matrix solver interface has not been initialized.

   **Notes:**
      The default value is 0.05; passing a non-positive value restores the
      default.

.. c:function:: int FIRKodeSetMassLSNormFactor(void* firkode_mem, sunrealtype nrmfac)

   The function :c:func:`FIRKodeSetMassLSNormFactor` specifies the factor to
   use when converting from the integrator tolerance (WRMS norm) to the mass
   matrix linear solver tolerance (L2 norm).

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nrmfac`` -- the norm conversion factor, interpreted as in :c:func:`FIRKodeSetLSNormFactor`.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_MASSMEM_NULL`` -- The mass matrix solver interface has not been initialized.
      - ``FIRKLS_ILL_INPUT`` -- A negative ``nrmfac`` was given but the vector does not implement :c:func:`N_VDotProd`.


.. _FIRKODE.Usage.optional_input.root:

Rootfinding optional input functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following functions can be called to set optional inputs to control the
rootfinding algorithm.

.. _FIRKODE.Usage.optional_input.root_table:

.. table:: Optional rootfinding inputs for FIRKODE

   +-------------------------------------+----------------------------------------------+----------------+
   |      **Optional input**             |              **Function name**               |  **Default**   |
   +=====================================+==============================================+================+
   | Direction of zero-crossing          | :c:func:`FIRKodeSetRootDirection`            | both           |
   +-------------------------------------+----------------------------------------------+----------------+
   | Disable rootfinding warnings        | :c:func:`FIRKodeSetNoInactiveRootWarn`       | none           |
   +-------------------------------------+----------------------------------------------+----------------+

.. c:function:: int FIRKodeSetRootDirection(void* firkode_mem, int* rootdir)

   The function :c:func:`FIRKodeSetRootDirection` specifies the direction of
   zero-crossings to be located and returned.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``rootdir`` -- state array of length ``nrtfn``, the number of root functions :math:`g_i`, as specified in the call to the function :c:func:`FIRKodeRootInit`. A value of 0 for ``rootdir[i]`` indicates that crossing in either direction for :math:`g_i` should be reported. A value of +1 or -1 indicates that the solver should report only zero-crossings where :math:`g_i` is increasing or decreasing, respectively.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.
      - ``FIRK_ILL_INPUT`` -- Rootfinding has not been activated through a call to :c:func:`FIRKodeRootInit`.

   **Notes:**
      The default behavior is to monitor for both zero-crossing directions.

.. c:function:: int FIRKodeSetNoInactiveRootWarn(void* firkode_mem)

   The function :c:func:`FIRKodeSetNoInactiveRootWarn` disables issuing a
   warning if some root function appears to be identically zero at the
   beginning of the integration.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.
      - ``FIRK_ILL_INPUT`` -- Rootfinding has not been activated through a call to :c:func:`FIRKodeRootInit`.

   **Notes:**
      FIRKODE will not report the initial conditions as a possible
      zero-crossing (assuming that one or more components :math:`g_i` are
      zero at the initial time). However, if it appears that some
      :math:`g_i` is identically zero at the initial time (i.e., :math:`g_i`
      is zero at the initial time *and* after the first step), FIRKODE will
      issue a warning which can be disabled with this optional input
      function.


.. _FIRKODE.Usage.dky:

Interpolated output function
----------------------------

An optional function :c:func:`FIRKodeGetDky` is available to obtain
additional output values. This function should only be called after a
successful return from :c:func:`FIRKodeEvolve` as it provides interpolated
values either of :math:`y` or of its derivatives (up to the number of stages
of the method) interpolated to any value of :math:`t` in the last internal
step taken by FIRKODE (see :numref:`FIRKODE.Mathematics.dense`).

.. c:function:: int FIRKodeGetDky(void* firkode_mem, sunrealtype t, int k, N_Vector dky)

   The function :c:func:`FIRKodeGetDky` computes the ``k``-th derivative of
   the function :math:`y` at time ``t``, i.e.
   :math:`\dfrac{\mathrm{d}^{k}y}{\mathrm{d}t^{k}}(t)`, where
   :math:`t_n - h_u \leq t \leq t_n`, :math:`t_n` denotes the current internal
   time reached, and :math:`h_u` is the last internal step size successfully
   used by the solver. The user may request ``k`` :math:`= 0, 1, \ldots, s`,
   where :math:`s` is the number of stages of the method.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``t`` -- the value of the independent variable at which the derivative is to be evaluated.
      - ``k`` -- the derivative order requested.
      - ``dky`` -- vector containing the derivative. This vector must be allocated by the user.

   **Return value:**
      - ``FIRK_SUCCESS`` -- :c:func:`FIRKodeGetDky` succeeded.
      - ``FIRK_BAD_K`` -- ``k`` is not in the range :math:`0, 1, \ldots, s`.
      - ``FIRK_BAD_T`` -- ``t`` is not in the interval :math:`[t_n - h_u, t_n]`.
      - ``FIRK_BAD_DKY`` -- The ``dky`` argument was ``NULL``.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

   **Notes:**
      It is only legal to call the function :c:func:`FIRKodeGetDky` after a
      successful return from :c:func:`FIRKodeEvolve`. Before the first step
      only the current state (``k = 0`` at :math:`t = t_0`) is available. See
      :c:func:`FIRKodeGetCurrentTime` and :c:func:`FIRKodeGetLastStep` in the
      next section for access to :math:`t_n` and :math:`h_u`, respectively.
      The interpolant is the collocation polynomial of the step, whose
      accuracy at interior points is of order :math:`s+1`, lower than the
      order :math:`2s-1` of the method at the step end points.


.. _FIRKODE.Usage.optional_output:

Optional output functions
-------------------------

FIRKODE provides an extensive set of functions that can be used to obtain
solver performance information. The optional outputs are divided into the
following groups: main solver outputs, nonlinear solver outputs, FIRKLS
linear solver interface outputs, mass matrix solver interface outputs, and
rootfinding outputs. Each of these is described in detail in the remainder
of this section.

.. _FIRKODE.Usage.optional_output.main:

Main solver optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: int FIRKodeGetNumSteps(void* firkode_mem, long int* nsteps)

   Returns the cumulative number of internal steps taken by the solver (total
   so far).

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nsteps`` -- number of steps taken by FIRKODE.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

.. c:function:: int FIRKodeGetNumStepAttempts(void* firkode_mem, long int* step_attempts)

   Returns the cumulative number of steps attempted by the solver (total so
   far), including rejected steps.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``step_attempts`` -- number of steps attempted by FIRKODE.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

.. c:function:: int FIRKodeGetNumRhsEvals(void* firkode_mem, long int* nfevals)

   Returns the number of calls to the user's right-hand side function.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nfevals`` -- number of calls to the user's ``f`` function.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

   **Notes:**
      The ``nfevals`` value returned by :c:func:`FIRKodeGetNumRhsEvals` does
      not account for calls made to ``f`` by a linear solver or preconditioner
      module (see :c:func:`FIRKodeGetNumLinRhsEvals`).

.. c:function:: int FIRKodeGetNumLinSolvSetups(void* firkode_mem, long int* nlinsetups)

   Returns the number of calls made to the linear solver's setup function.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nlinsetups`` -- number of calls made to the linear solver setup function.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

.. c:function:: int FIRKodeGetNumErrTestFails(void* firkode_mem, long int* netfails)

   Returns the number of local error test failures that have occurred.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``netfails`` -- number of error test failures.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

.. c:function:: int FIRKodeGetNumStepSolveFails(void* firkode_mem, long int* nncfails)

   Returns the number of failed steps due to a nonlinear solver failure.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nncfails`` -- number of step failures.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

.. c:function:: int FIRKodeGetActualInitStep(void* firkode_mem, sunrealtype* hinused)

   Returns the value of the integration step size used on the first step.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``hinused`` -- actual value of initial step size.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

   **Notes:**
      Even if the value of the initial integration step size was specified by
      the user through a call to :c:func:`FIRKodeSetInitStep`, this value
      might have been changed by FIRKODE to ensure that the step size is
      within the prescribed bounds (:math:`h_{min} \le h_0 \le h_{max}`), or
      to satisfy the local error test condition.

.. c:function:: int FIRKodeGetLastStep(void* firkode_mem, sunrealtype* hlast)

   Returns the integration step size taken on the last internal step.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``hlast`` -- step size taken on the last internal step.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

.. c:function:: int FIRKodeGetCurrentStep(void* firkode_mem, sunrealtype* hcur)

   Returns the integration step size to be attempted on the next internal
   step.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``hcur`` -- step size to be attempted on the next internal step.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

.. c:function:: int FIRKodeGetCurrentTime(void* firkode_mem, sunrealtype* tcur)

   Returns the current internal time reached by the solver.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``tcur`` -- current internal time reached.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

.. c:function:: int FIRKodeGetCurrentState(void* firkode_mem, N_Vector* y)

   Returns a pointer to the current state vector, i.e., the solution at the
   current internal time.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``y`` -- pointer to the current state vector.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

   **Notes:**
      The vector is owned by FIRKODE and must not be modified or freed by the
      user. The pointer may change between steps.

.. c:function:: int FIRKodeGetCurrentGamma(void* firkode_mem, sunrealtype* gamma)

   Returns the current value of the scalar :math:`\gamma` in the matrix
   :math:`M - \gamma J`.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``gamma`` -- the current value of :math:`\gamma`.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

.. c:function:: int FIRKodeGetTolScaleFactor(void* firkode_mem, sunrealtype* tolsfac)

   Returns a suggested factor by which the user's tolerances should be
   scaled when too much accuracy has been requested for some internal step.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``tolsfac`` -- suggested scaling factor for user-supplied tolerances.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

.. c:function:: int FIRKodeGetErrWeights(void* firkode_mem, N_Vector eweight)

   Returns the solution error weights at the current time. These are the
   reciprocals of the :math:`W_i` given by :eq:`FIRKODE_errwt`.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``eweight`` -- solution error weights at the current time.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

   **Notes:**
      The user must allocate memory for ``eweight``.

.. c:function:: int FIRKodeGetEstLocalErrors(void* firkode_mem, N_Vector ele)

   Returns the vector of estimated local errors of the last step attempt.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``ele`` -- estimated local errors.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

   **Notes:**
      The user must allocate memory for ``ele``. The values returned in
      ``ele`` are valid only if :c:func:`FIRKodeEvolve` returned a
      non-negative value and a step was taken with time step adaptivity
      enabled. The ``ele`` vector, together with the ``eweight`` vector from
      :c:func:`FIRKodeGetErrWeights`, can be used to determine how the
      various components of the system contributed to the estimated local
      error test. Specifically, that error test uses the RMS norm of a vector
      whose components are the products of the components of these two
      vectors. Thus, for example, if there were recent error test failures,
      the components causing the failures are those with largest values for
      the products, denoted loosely as ``eweight[i]*ele[i]``.

.. c:function:: int FIRKodeGetIntegratorStats(void* firkode_mem, long int* nsteps, long int* nfevals, long int* nlinsetups, long int* netfails, sunrealtype* hinused, sunrealtype* hlast, sunrealtype* hcur, sunrealtype* tcur)

   Returns the main integrator statistics as a group.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nsteps`` -- number of steps taken by FIRKODE.
      - ``nfevals`` -- number of calls to the user's ``f`` function.
      - ``nlinsetups`` -- number of calls made to the linear solver setup function.
      - ``netfails`` -- number of error test failures.
      - ``hinused`` -- actual value of initial step size.
      - ``hlast`` -- step size taken on the last internal step.
      - ``hcur`` -- step size to be attempted on the next internal step.
      - ``tcur`` -- current internal time reached.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output values have been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

.. c:function:: int FIRKodeGetUserData(void* firkode_mem, void** user_data)

   Returns the user data pointer set with :c:func:`FIRKodeSetUserData`.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``user_data`` -- the user data pointer.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

.. c:function:: int FIRKodePrintAllStats(void* firkode_mem, FILE* outfile, SUNOutputFormat fmt)

   Outputs all of the integrator, nonlinear solver, linear solver, mass
   matrix solver, and rootfinding statistics.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``outfile`` -- pointer to output file.
      - ``fmt`` -- the output format: ``SUN_OUTPUTFORMAT_TABLE`` for a human-readable table or ``SUN_OUTPUTFORMAT_CSV`` for comma-separated values.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The output was successfully written.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.
      - ``FIRK_ILL_INPUT`` -- An invalid formatting option was provided.

.. c:function:: int FIRKodeWriteParameters(void* firkode_mem, FILE* fp)

   Outputs all FIRKODE solver parameters to the provided file pointer.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``fp`` -- pointer to the output file.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The parameters were successfully written.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

.. c:function:: char* FIRKodeGetReturnFlagName(long int flag)

   Returns the name of the FIRKODE constant corresponding to ``flag``.

   **Arguments:**
      - ``flag`` -- the flag returned by a call to a FIRKODE function.

   **Return value:**
      - A string containing the name of the corresponding constant. The
        string is allocated with ``malloc`` and must be freed by the user.


.. _FIRKODE.Usage.optional_output.nls:

Nonlinear solver optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: int FIRKodeGetNumNonlinSolvIters(void* firkode_mem, long int* nniters)

   Returns the number of Newton iterations performed.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nniters`` -- number of Newton iterations.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

.. c:function:: int FIRKodeGetNumNonlinSolvConvFails(void* firkode_mem, long int* nnfails)

   Returns the number of nonlinear solver convergence failures.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nnfails`` -- number of nonlinear convergence failures.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

.. c:function:: int FIRKodeGetNonlinSolvStats(void* firkode_mem, long int* nniters, long int* nnfails)

   Returns the nonlinear solver statistics as a group.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nniters`` -- number of Newton iterations.
      - ``nnfails`` -- number of nonlinear convergence failures.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output values have been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.


.. _FIRKODE.Usage.optional_output.ls:

Linear solver interface optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following optional outputs are available from the FIRKLS interface. In
the descriptions below, "block" linear solves refer to solves with the
matrices :math:`M - \gamma J` performed by the user's linear solver, and
"stage" linear iterations refer to the internal flexible GMRES iteration on
the coupled stage system.

.. c:function:: int FIRKodeGetJac(void* firkode_mem, SUNMatrix* J)

   Returns the internally stored copy of the Jacobian matrix of the ODE
   right-hand side function.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``J`` -- the Jacobian matrix.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

   **Notes:**
      This function is provided for debugging purposes and the values in the
      returned matrix should not be altered. The pointer is ``NULL`` for
      matrix-free linear solvers.

.. c:function:: int FIRKodeGetJacTime(void* firkode_mem, sunrealtype* t_J)

   Returns the time at which the internally stored copy of the Jacobian
   matrix was evaluated.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``t_J`` -- the time at which the Jacobian was evaluated.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

.. c:function:: int FIRKodeGetJacNumSteps(void* firkode_mem, long int* nst_J)

   Returns the value of the internal step counter at the time the internally
   stored copy of the Jacobian matrix was evaluated.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nst_J`` -- the value of the internal step counter.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

.. c:function:: int FIRKodeGetNumJacEvals(void* firkode_mem, long int* njevals)

   Returns the number of calls made to the Jacobian approximation function.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``njevals`` -- the number of calls to the Jacobian function.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

.. c:function:: int FIRKodeGetNumLinRhsEvals(void* firkode_mem, long int* nfevalsLS)

   Returns the number of calls made to the user-supplied right-hand side
   function due to the finite difference Jacobian approximation or finite
   difference Jacobian-vector product approximation.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nfevalsLS`` -- the number of calls made to the user-supplied right-hand side function.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

   **Notes:**
      The value ``nfevalsLS`` is incremented only if the default internal
      difference quotient function is used.

.. c:function:: int FIRKodeGetNumPrecEvals(void* firkode_mem, long int* npevals)

   Returns the number of preconditioner evaluations, i.e., the number of
   calls made to ``pset`` with ``jok = SUNFALSE``.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``npevals`` -- the current number of calls to ``pset``.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

.. c:function:: int FIRKodeGetNumPrecSolves(void* firkode_mem, long int* npsolves)

   Returns the cumulative number of calls made to the preconditioner solve
   function, ``psolve``.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``npsolves`` -- the current number of calls to ``psolve``.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

.. c:function:: int FIRKodeGetNumLinIters(void* firkode_mem, long int* nliters)

   Returns the cumulative number of linear iterations performed by the user's
   iterative linear solver in block solves.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nliters`` -- the current number of linear iterations.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

.. c:function:: int FIRKodeGetNumLinConvFails(void* firkode_mem, long int* nlcfails)

   Returns the cumulative number of convergence failures of the user's
   iterative linear solver in block solves.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nlcfails`` -- the current number of linear convergence failures.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

.. c:function:: int FIRKodeGetNumJTSetupEvals(void* firkode_mem, long int* njtsetups)

   Returns the cumulative number of calls made to the Jacobian-vector setup
   function, ``jtsetup``.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``njtsetups`` -- the current number of calls to ``jtsetup``.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

.. c:function:: int FIRKodeGetNumJtimesEvals(void* firkode_mem, long int* njvevals)

   Returns the cumulative number of calls made to the Jacobian-vector
   function, ``jtimes`` (user-supplied or difference quotient).

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``njvevals`` -- the current number of calls to ``jtimes``.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

.. c:function:: int FIRKodeGetNumBlockSolves(void* firkode_mem, long int* nbsolves)

   Returns the cumulative number of block linear solves, i.e., solves with
   the matrices :math:`M - \gamma J` performed by the user's linear solver
   (as preconditioner applications and for the error estimate).

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nbsolves`` -- the current number of block linear solves.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

.. c:function:: int FIRKodeGetNumStageLinIters(void* firkode_mem, long int* nsliters)

   Returns the cumulative number of iterations of the internal flexible GMRES
   solver on the coupled stage system.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nsliters`` -- the current number of stage linear iterations.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

.. c:function:: int FIRKodeGetNumStageLinConvFails(void* firkode_mem, long int* nslcfails)

   Returns the cumulative number of convergence failures of the internal
   flexible GMRES solver on the coupled stage system.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nslcfails`` -- the current number of stage linear convergence failures.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

.. c:function:: int FIRKodeGetLinSolveStats(void* firkode_mem, long int* njevals, long int* nfevalsLS, long int* nliters, long int* nlcfails, long int* npevals, long int* npsolves, long int* njtsetups, long int* njtimes)

   Returns the linear solver statistics as a group.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``njevals`` -- the number of calls to the Jacobian function.
      - ``nfevalsLS`` -- the number of calls to the right-hand side function for difference quotients.
      - ``nliters`` -- the number of block linear iterations.
      - ``nlcfails`` -- the number of block linear convergence failures.
      - ``npevals`` -- the number of calls to ``pset``.
      - ``npsolves`` -- the number of calls to ``psolve``.
      - ``njtsetups`` -- the number of calls to ``jtsetup``.
      - ``njtimes`` -- the number of calls to ``jtimes``.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output values have been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

.. c:function:: int FIRKodeGetLastLinFlag(void* firkode_mem, long int* flag)

   Returns the last return value from a FIRKLS routine.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``flag`` -- the value of the last return flag from a FIRKLS function.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

   **Notes:**
      If the FIRKLS setup function failed when using the
      :ref:`SUNLINSOL_DENSE <SUNLinSol_Dense>` or :ref:`SUNLINSOL_BAND
      <SUNLinSol_Band>` modules, then the value of ``flag`` is equal to the
      column index (numbered from one) at which a zero diagonal element was
      encountered during the LU factorization of the (dense or banded)
      Jacobian matrix. For all other failures, ``flag`` is negative.

.. c:function:: char* FIRKodeGetLinReturnFlagName(long int flag)

   Returns the name of the FIRKLS constant corresponding to ``flag``.

   **Arguments:**
      - ``flag`` -- the flag returned by a call to a FIRKLS function.

   **Return value:**
      - A string containing the name of the corresponding constant. If
        ``1 <= flag <= N`` (LU factorization failed), this function returns
        "NONE". The string is allocated with ``malloc`` and must be freed by
        the user.


.. _FIRKODE.Usage.optional_output.mass:

Mass matrix solver interface optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All of the following functions return ``FIRKLS_MASSMEM_NULL`` if the mass
matrix solver interface has not been initialized with
:c:func:`FIRKodeSetMassLinearSolver`.

.. c:function:: int FIRKodeGetCurrentMassMatrix(void* firkode_mem, SUNMatrix* M)

   Returns the mass matrix object supplied to :c:func:`FIRKodeSetMassLinearSolver`.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``M`` -- the mass matrix (``NULL`` for a matrix-free mass linear solver).

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_MASSMEM_NULL`` -- The mass matrix solver interface has not been initialized.

.. c:function:: int FIRKodeGetNumMassSetups(void* firkode_mem, long int* nmsetups)

   Returns the number of calls made to the mass matrix linear solver setup
   function.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nmsetups`` -- number of calls to the mass matrix solver setup function.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_MASSMEM_NULL`` -- The mass matrix solver interface has not been initialized.

.. c:function:: int FIRKodeGetNumMassMultSetups(void* firkode_mem, long int* nmvsetups)

   Returns the number of calls made to the mass matrix-vector product setup
   routines (the user-supplied ``mtsetup`` routine and the
   :c:func:`SUNMatMatvecSetup` routine of the mass matrix).

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nmvsetups`` -- number of calls to the mass matrix-vector product setup routines.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_MASSMEM_NULL`` -- The mass matrix solver interface has not been initialized.

.. c:function:: int FIRKodeGetNumMassMult(void* firkode_mem, long int* nmvevals)

   Returns the number of mass matrix-vector products computed (with the
   user-supplied ``mtimes`` routine or with :c:func:`SUNMatMatvec`).

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nmvevals`` -- number of mass matrix-vector products.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_MASSMEM_NULL`` -- The mass matrix solver interface has not been initialized.

.. c:function:: int FIRKodeGetNumMassSolves(void* firkode_mem, long int* nmsolves)

   Returns the number of linear solves with the mass matrix.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nmsolves`` -- number of mass matrix linear solves.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_MASSMEM_NULL`` -- The mass matrix solver interface has not been initialized.

.. c:function:: int FIRKodeGetNumMassPrecEvals(void* firkode_mem, long int* nmpevals)

   Returns the number of calls made to the mass matrix preconditioner setup
   function.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nmpevals`` -- number of calls to the mass matrix preconditioner setup function.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_MASSMEM_NULL`` -- The mass matrix solver interface has not been initialized.

.. c:function:: int FIRKodeGetNumMassPrecSolves(void* firkode_mem, long int* nmpsolves)

   Returns the number of calls made to the mass matrix preconditioner solve
   function.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nmpsolves`` -- number of calls to the mass matrix preconditioner solve function.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_MASSMEM_NULL`` -- The mass matrix solver interface has not been initialized.

.. c:function:: int FIRKodeGetNumMassIters(void* firkode_mem, long int* nmiters)

   Returns the cumulative number of linear iterations performed by an
   iterative mass matrix linear solver.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nmiters`` -- number of mass matrix linear iterations.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_MASSMEM_NULL`` -- The mass matrix solver interface has not been initialized.

.. c:function:: int FIRKodeGetNumMassConvFails(void* firkode_mem, long int* nmcfails)

   Returns the cumulative number of convergence failures of the mass matrix
   linear solver.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nmcfails`` -- number of mass matrix linear solver convergence failures.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_MASSMEM_NULL`` -- The mass matrix solver interface has not been initialized.

.. c:function:: int FIRKodeGetLastMassFlag(void* firkode_mem, long int* flag)

   Returns the last return value from a mass matrix solver interface routine.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``flag`` -- the value of the last return flag.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_MASSMEM_NULL`` -- The mass matrix solver interface has not been initialized.


.. _FIRKODE.Usage.optional_output.root:

Rootfinding optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: int FIRKodeGetRootInfo(void* firkode_mem, int* rootsfound)

   Returns an array showing which functions were found to have a root.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``rootsfound`` -- array of length ``nrtfn`` with the indices of the user functions :math:`g_i` found to have a root. For :math:`i = 0, \ldots,` ``nrtfn`` :math:`-1`, ``rootsfound[i]`` is nonzero if :math:`g_i` has a root, and 0 if not.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.
      - ``FIRK_ILL_INPUT`` -- Rootfinding has not been activated through a call to :c:func:`FIRKodeRootInit`.

   **Notes:**
      Note that, for the components :math:`g_i` for which a root was found,
      the sign of ``rootsfound[i]`` indicates the direction of
      zero-crossing. A value of +1 indicates that :math:`g_i` is increasing,
      while a value of -1 indicates a decreasing :math:`g_i`. The user must
      allocate memory for the vector ``rootsfound``.

.. c:function:: int FIRKodeGetNumGEvals(void* firkode_mem, long int* ngevals)

   Returns the cumulative number of calls made to the user's root function
   :math:`g`.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``ngevals`` -- number of calls made to the user's function :math:`g` thus far (zero if rootfinding is not active).

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.


.. _FIRKODE.Usage.reinit:

FIRKODE reinitialization and reset functions
--------------------------------------------

The function :c:func:`FIRKodeReInit` reinitializes the main FIRKODE solver
for the solution of a new problem, where a prior call to :c:func:`FIRKodeInit`
has been made. The new problem must have the same size as the previous one.
:c:func:`FIRKodeReInit` performs the same input checking and initializations
that :c:func:`FIRKodeInit` does, but does no memory allocation as it assumes
that the existing internal memory is sufficient for the new problem. A call
to :c:func:`FIRKodeReInit` resets all integrator statistics and, if a
different number of stages or a different Krylov subspace dimension has been
selected, reallocates the internal stage and Krylov storage at the next call
to :c:func:`FIRKodeEvolve`.

The use of :c:func:`FIRKodeReInit` requires that the number of Runge--Kutta
stages, the ``SUNLinearSolver`` and ``SUNMatrix`` objects, and any
user-supplied functions remain valid for the new problem; these may be
changed with the corresponding ``Set`` functions before the next call to
:c:func:`FIRKodeEvolve`. If the tolerances are to be changed, the appropriate
tolerance function must be called after :c:func:`FIRKodeReInit`.

.. c:function:: int FIRKodeReInit(void* firkode_mem, sunrealtype t0, N_Vector y0)

   The function :c:func:`FIRKodeReInit` provides required problem
   specifications and reinitializes FIRKODE.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``t0`` -- is the initial value of :math:`t`.
      - ``y0`` -- is the initial value of :math:`y`.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The call was successful.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.
      - ``FIRK_NO_MALLOC`` -- Memory space for the FIRKODE memory block was not allocated through a previous call to :c:func:`FIRKodeInit`.
      - ``FIRK_ILL_INPUT`` -- An input argument to :c:func:`FIRKodeReInit` has an illegal value.

   **Notes:**
      If an error occurred, :c:func:`FIRKodeReInit` also sends an error
      message to the error handler function.

.. c:function:: int FIRKodeReset(void* firkode_mem, sunrealtype tR, N_Vector yR)

   The function :c:func:`FIRKodeReset` resets the integrator state to the
   given independent variable value and dependent variable vector, while
   retaining all optional inputs and integrator statistics.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``tR`` -- the value of the independent variable :math:`t`.
      - ``yR`` -- the value of the dependent variable vector :math:`y(t_R)`.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The call was successful.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.
      - ``FIRK_NO_MALLOC`` -- Memory space for the FIRKODE memory block was not allocated through a previous call to :c:func:`FIRKodeInit`.
      - ``FIRK_ILL_INPUT`` -- An input argument to :c:func:`FIRKodeReset` has an illegal value.

   **Notes:**
      By default the next call to :c:func:`FIRKodeEvolve` will use the
      initial step size estimation heuristic; a step size may be specified
      with :c:func:`FIRKodeSetInitStep`. All counters are retained and the
      history used by the stage predictor is discarded.
