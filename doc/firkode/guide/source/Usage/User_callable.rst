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

      A second call to :c:func:`FIRKodeInit` on the same memory block frees
      the internal vectors of the first and discards the tolerances, which
      must be set again before :c:func:`FIRKodeEvolve`; to continue a problem
      of the same size use :c:func:`FIRKodeReInit` or :c:func:`FIRKodeReset`.

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
      - ``FIRK_MEM_FAIL`` -- The copy of the absolute tolerance vector could not be allocated.

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
by a factor derived from the relative tolerance as in RADAU5, or set with
:c:func:`FIRKodeSetNonlinConvCoef`. Since the
tolerances control local errors, the global error at the end of a long
integration is typically larger than the tolerances suggest; the actual
error can be assessed by comparing runs with different tolerances.


.. _FIRKODE.Usage.mode:

Mode selection functions
------------------------

FIRKODE has two modes. In BEGINNER mode, the default, only the essential
functions are available: creation, initialization and deallocation, the
tolerances, the user data, the linear and mass-matrix solver interfaces with
their user-supplied functions, the stop time, the step-size bounds and the
initial and fixed step size, the maximum number of steps, rootfinding, the
solver and dense-output functions, every getter and statistic, the output
functions, the warning limit (:c:func:`FIRKodeSetMaxHnilWarns`), the
tolerances and scaling of the user's linear solver
(:c:func:`FIRKodeSetEpsLin`, :c:func:`FIRKodeSetLSNormFactor`,
:c:func:`FIRKodeSetLinearSolutionScaling` and their mass-matrix
counterparts), :c:func:`FIRKodeSetOptions` for the accepted keys, and
:c:func:`FIRKodeSetMode`. Every other optional input changes how the method,
the Newton iteration, the step-size control or the stage solver works, and is
accepted only in EXPERT mode; in BEGINNER mode it returns ``FIRK_ILL_INPUT``,
so that a BEGINNER run always uses the validated defaults. The mode column of
the tables below marks each optional input B or E.

.. c:function:: int FIRKodeSetMode(void* firkode_mem, int mode)

   Selects the mode of FIRKODE.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``mode`` -- ``FIRK_BEGINNER`` (the default) or ``FIRK_EXPERT``.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The call was successful.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.
      - ``FIRK_ILL_INPUT`` -- The mode is not one of the two values, or a step
        has been attempted since :c:func:`FIRKodeCreate` or
        :c:func:`FIRKodeReInit`.
      - ``FIRK_LINIT_FAIL`` -- The options of the stacked stage solver could
        not be restored on the return to BEGINNER mode.
      - ``FIRK_MEM_FAIL`` or ``FIRK_CONTROLLER_ERR`` -- The default step size
        controller could not be created on the return to BEGINNER mode.

   **Notes:**
      The mode can be set at any time before the first step, and again after
      :c:func:`FIRKodeInit` or :c:func:`FIRKodeReInit`, but not after
      :c:func:`FIRKodeReset`, which keeps the history of the integration. A return to BEGINNER mode
      restores the default of every EXPERT option; if that fails, the mode
      is unchanged. :c:func:`FIRKodeWriteParameters` prints the mode first.

.. c:function:: int FIRKodeGetMode(void* firkode_mem, int* mode_ptr)

   Returns the mode of FIRKODE.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``mode_ptr`` -- the mode, ``FIRK_BEGINNER`` or ``FIRK_EXPERT``.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The call was successful.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.


.. _FIRKODE.Usage.method:

Method selection functions
--------------------------

FIRKODE integrates with the three-stage Radau IIA method of order five, the
only method of this version. The following functions, which select the method
by its number of stages or its order, accept only this method; they must be
called before the first call to :c:func:`FIRKodeEvolve` or directly after
:c:func:`FIRKodeReInit`.

.. c:function:: int FIRKodeSetNumStages(void* firkode_mem, int s)

   The function :c:func:`FIRKodeSetNumStages` selects the :math:`s`-stage
   Radau IIA method, of order :math:`2s-1`. This version of FIRKODE
   provides only the three-stage method of order 5, the method of RADAU5;
   the function exists so that its signature is stable when later versions
   add the other stage counts.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``s`` -- the number of stages, which must be 3.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The call was successful.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.
      - ``FIRK_ILL_INPUT`` -- ``s`` is not 3.
      - ``FIRK_ILL_INPUT`` -- The mode is BEGINNER; this is an EXPERT option (see :c:func:`FIRKodeSetMode`).

   **Notes:**
      The default is :math:`s = 3`. The tables of the Radau IIA methods with
      :math:`1 \le s \le 9` stages are nevertheless available through
      :c:func:`FIRKodeTable_RadauIIA` (see :numref:`FIRKODE.Usage.Tables`).

.. c:function:: int FIRKodeSetOrder(void* firkode_mem, int ord)

   The function :c:func:`FIRKodeSetOrder` selects the Radau IIA method with
   the smallest number of stages whose order is at least ``ord``, i.e.,
   :math:`s = \lceil (\mathrm{ord}+1)/2 \rceil`. As this version of FIRKODE
   provides only :math:`s = 3`, it accepts only the orders 4 and 5.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``ord`` -- the requested order of accuracy, 4 or 5; a non-positive value selects the default method.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The call was successful.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.
      - ``FIRK_ILL_INPUT`` -- ``ord`` selects a number of stages other than 3.
      - ``FIRK_ILL_INPUT`` -- The mode is BEGINNER; this is an EXPERT option (see :c:func:`FIRKodeSetMode`).

.. c:function:: int FIRKodeGetNumStages(void* firkode_mem, int* s_ptr)

   The function :c:func:`FIRKodeGetNumStages` returns the number of stages
   of the current method.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``s_ptr`` -- the number of stages.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The call was successful.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

.. c:function:: int FIRKodeGetCurrentTable(void* firkode_mem, FIRKodeTable* tbl_ptr)

   The function :c:func:`FIRKodeGetCurrentTable` returns a pointer to the
   coefficient table of the current method (see :numref:`FIRKODE.Usage.Tables`).

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``tbl_ptr`` -- the coefficient table.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The call was successful.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

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
      implementation provides this operation; otherwise with a user-supplied
      Jacobian-times-vector routine (see :c:func:`FIRKodeSetJacTimes`) or,
      without one, by a difference quotient, with a warning when the linear
      solver is initialized.

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
      :c:func:`FIRKodeReInit`, and reused for the rest of the integration,
      also after :c:func:`FIRKodeReset`. Attaching a mass routine again
      (:c:func:`FIRKodeSetMassFn`, :c:func:`FIRKodeSetMassTimes` or
      :c:func:`FIRKodeSetMassPreconditioner`) has it set up anew.

      Besides the initial step size estimate and rootfinding, the mass
      linear solver is used to convert the right-hand sides of the iterative
      linear solves to the units of :math:`y` (see
      :numref:`FIRKODE.Mathematics.mass`): for a method with :math:`s > 1`
      stages it performs :math:`s` solves in every Newton iteration, so it
      should be efficient.


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

.. c:function:: int FIRKodeEvolve(void* firkode_mem, sunrealtype tout, N_Vector yout, sunrealtype* tret_ptr, int itask)

   The function :c:func:`FIRKodeEvolve` integrates the ODE over an interval in :math:`t`.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``tout`` -- the next time at which a computed solution is desired.
      - ``yout`` -- the computed solution vector.
      - ``tret_ptr`` -- the time reached by the solver (output).
      - ``itask`` -- a flag indicating the job of the solver for the next user step. The ``FIRK_NORMAL`` option causes the solver to take internal steps until it has reached or just passed the user-specified ``tout`` parameter. The solver then interpolates in order to return an approximate value of :math:`y(t_{out})`; the interpolant is the collocation polynomial of degree :math:`s`, whose accuracy (order :math:`s + 1`) is below that of the method at the step ends (order :math:`2s - 1`), so for a full-order value at a final time make it a stop time with :c:func:`FIRKodeSetStopTime`. The ``FIRK_ONE_STEP`` option tells the solver to take just one internal step and then return the solution at the point reached by that step.

   **Return value:**
      - ``FIRK_SUCCESS`` -- :c:func:`FIRKodeEvolve` succeeded and no roots were found.
      - ``FIRK_TSTOP_RETURN`` -- :c:func:`FIRKodeEvolve` succeeded by reaching the stopping point specified through the optional input function :c:func:`FIRKodeSetStopTime`.
      - ``FIRK_ROOT_RETURN`` -- :c:func:`FIRKodeEvolve` succeeded and found one or more roots. In this case, ``tret_ptr`` is the location of the root. If ``nrtfn`` :math:`> 1`, call :c:func:`FIRKodeGetRootInfo` to see which :math:`g_i` were found to have a root.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized through a previous call to :c:func:`FIRKodeCreate`.
      - ``FIRK_NO_MALLOC`` -- The FIRKODE memory was not allocated by a call to :c:func:`FIRKodeInit`.
      - ``FIRK_ILL_INPUT`` -- One of the inputs to :c:func:`FIRKodeEvolve` was illegal, or some other input to the solver was illegal or missing. The latter category includes the following situations: the tolerances have not been set; no linear solver has been attached; a component of the error weight vector became zero during internal time-stepping; the direction of a fixed step size disagrees with ``tout``; a root of one of the root functions was found both at a point :math:`t` and also very near :math:`t`.
      - ``FIRK_TOO_CLOSE`` -- The initial time :math:`t_0` and the output time :math:`t_{out}` are too close to each other and the user did not specify an initial step size.
      - ``FIRK_TOO_MUCH_WORK`` -- The solver took ``mxstep`` internal steps but still could not reach ``tout``. The default value for ``mxstep`` is 500.
      - ``FIRK_TOO_MUCH_ACC`` -- The solver could not satisfy the accuracy demanded by the user for some internal step.
      - ``FIRK_ERR_FAILURE`` -- Either error test failures occurred too many times (default 7) during one internal time step, or with :math:`|h| = h_{min}`.
      - ``FIRK_CONV_FAILURE`` -- Either convergence test failures occurred too many times (default 10) during one internal time step, or with :math:`|h| = h_{min}`, or, in fixed-step mode, the Newton iteration failed with current Jacobian data (a failure with older data is retried once with a new Jacobian).
      - ``FIRK_LINIT_FAIL`` -- The linear solver interface's initialization function failed.
      - ``FIRK_LSETUP_FAIL`` -- The linear solver interface's setup function failed in an unrecoverable manner.
      - ``FIRK_LSOLVE_FAIL`` -- The linear solver interface's solve function failed in an unrecoverable manner.
      - ``FIRK_RHSFUNC_FAIL`` -- The right-hand side function failed in an unrecoverable manner.
      - ``FIRK_FIRST_RHSFUNC_ERR`` -- The right-hand side function had a recoverable error at the first call.
      - ``FIRK_REPTD_RHSFUNC_ERR`` -- Convergence test failures occurred too many times due to repeated recoverable errors in the right-hand side function. This flag will also be returned if the right-hand side function had repeated recoverable errors during the estimation of an initial step size.
      - ``FIRK_UNREC_RHSFUNC_ERR`` -- The right-hand function had a recoverable error, but no recovery was possible. This occurs when the right-hand side function fails recoverably when evaluated at the beginning of a step.
      - ``FIRK_RTFUNC_FAIL`` -- The rootfinding function failed.
      - ``FIRK_MASSINIT_FAIL`` -- The mass matrix solver interface's initialization function failed.
      - ``FIRK_MASSSETUP_FAIL`` -- The mass matrix solver interface's setup function failed.
      - ``FIRK_MASSSOLVE_FAIL`` -- The mass matrix solve failed in an unrecoverable manner.
      - ``FIRK_MASSMULT_FAIL`` -- The mass matrix-vector product failed.
      - ``FIRK_TABLE_FAIL`` -- The coefficient table of the selected method could not be constructed.
      - ``FIRK_CONTROLLER_ERR`` -- A call to the time step adaptivity controller failed.
      - ``FIRK_VECTOROP_ERR`` -- A vector operation failed.
      - ``FIRK_MEM_FAIL`` -- A memory allocation of the initial setup failed.
      - ``FIRK_UNRECOGNIZED_ERR`` -- An internal failure with no code of its own; please report it.

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
      :c:func:`FIRKodeEvolve`, the returned values of ``tret_ptr`` and ``yout``
      correspond to the farthest point reached in the integration. On all
      other error returns, ``tret_ptr`` and ``yout`` are left unchanged from the
      previous :c:func:`FIRKodeEvolve` return.


.. _FIRKODE.Usage.optional_input:

Optional input functions
------------------------

There are numerous optional input parameters that control the behavior of
the FIRKODE solver. FIRKODE provides functions that can be used to change
these optional input parameters from their default values. The **Mode**
column of the tables below says whether a function is available in BEGINNER
mode (B) or only in EXPERT mode (E); in BEGINNER mode an E function returns
``FIRK_ILL_INPUT`` (see :numref:`FIRKODE.Usage.mode`). The main inputs are
divided into the following categories:

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

   +-------------------------------+----------------------------------------------+----------------+----------+
   |      **Optional input**       |              **Function name**               |  **Default**   | **Mode** |
   +===============================+==============================================+================+==========+
   | User data                     | :c:func:`FIRKodeSetUserData`                 | ``NULL``       | B        |
   +-------------------------------+----------------------------------------------+----------------+----------+
   | Maximum no. of internal steps | :c:func:`FIRKodeSetMaxNumSteps`              | 500            | B        |
   | before :math:`t_{out}`        |                                              |                |          |
   +-------------------------------+----------------------------------------------+----------------+----------+
   | Maximum no. of warnings for   | :c:func:`FIRKodeSetMaxHnilWarns`             | 10             | B        |
   | :math:`t_n+h=t_n`             |                                              |                |          |
   +-------------------------------+----------------------------------------------+----------------+----------+
   | Initial step size             | :c:func:`FIRKodeSetInitStep`                 | estimated      | B        |
   +-------------------------------+----------------------------------------------+----------------+----------+
   | Minimum absolute step size    | :c:func:`FIRKodeSetMinStep`                  | 0.0            | B        |
   +-------------------------------+----------------------------------------------+----------------+----------+
   | Maximum absolute step size    | :c:func:`FIRKodeSetMaxStep`                  | :math:`\infty` | B        |
   +-------------------------------+----------------------------------------------+----------------+----------+
   | Fixed step size               | :c:func:`FIRKodeSetFixedStep`                | disabled       | B        |
   +-------------------------------+----------------------------------------------+----------------+----------+
   | Value of :math:`t_{stop}`     | :c:func:`FIRKodeSetStopTime`                 | undefined      | B        |
   +-------------------------------+----------------------------------------------+----------------+----------+
   | Interpolate at                | :c:func:`FIRKodeSetInterpolateStopTime`      | ``SUNFALSE``   | B        |
   | :math:`t_{stop}`              |                                              |                |          |
   +-------------------------------+----------------------------------------------+----------------+----------+
   | Disable the stop time         | :c:func:`FIRKodeClearStopTime`               | N/A            | B        |
   +-------------------------------+----------------------------------------------+----------------+----------+
   | Maximum no. of error test     | :c:func:`FIRKodeSetMaxErrTestFails`          | 7              | E        |
   | failures                      |                                              |                |          |
   +-------------------------------+----------------------------------------------+----------------+----------+
   | Maximum no. of convergence    | :c:func:`FIRKodeSetMaxConvFails`             | 10             | E        |
   | failures                      |                                              |                |          |
   +-------------------------------+----------------------------------------------+----------------+----------+
   | Jacobian and setup reuse      | :c:func:`FIRKodeSetReusePolicy`              | RADAU5         | E        |
   | policy                        |                                              |                |          |
   +-------------------------------+----------------------------------------------+----------------+----------+
   | Linear solver setup frequency | :c:func:`FIRKodeSetLSetupFrequency`          | 20             | E        |
   +-------------------------------+----------------------------------------------+----------------+----------+
   | Max change in :math:`\gamma`  | :c:func:`FIRKodeSetDeltaGammaMaxLSetup`      | 0.3            | E        |
   | without a linear solver setup |                                              |                |          |
   +-------------------------------+----------------------------------------------+----------------+----------+
   | Stage predictor               | :c:func:`FIRKodeSetPredictorMethod`          | extrapolation  | E        |
   +-------------------------------+----------------------------------------------+----------------+----------+
   | Second filtering of the error | :c:func:`FIRKodeSetErrorRefilter`            | ``SUNTRUE``    | E        |
   | estimate                      |                                              |                |          |
   +-------------------------------+----------------------------------------------+----------------+----------+

.. c:function:: int FIRKodeSetOptions(void* firkode_mem, const char* firkid, const char* file_name, int argc, char* argv[])

   The function :c:func:`FIRKodeSetOptions` sets FIRKODE options from an array
   of strings, typically the command-line arguments of the program.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``firkid`` -- the prefix of the options to read; ``NULL`` or ``""`` selects ``"firkode"``.
      - ``file_name`` -- must be ``NULL`` or ``""``: file-based options are not supported.
      - ``argc`` -- the number of strings in ``argv``.
      - ``argv`` -- the strings; ``argv[0]`` is skipped, as in the arguments of ``main``.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The options were set.
      - ``FIRK_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRK_ILL_INPUT`` -- A file name was given, a key lacks its value, or a value is not one of the names accepted.
      - ``FIRK_MEM_FAIL`` -- A memory allocation failed.
      - Any other value is the return value of the setter of a key, e.g. ``FIRK_ILL_INPUT`` for an EXPERT option in BEGINNER mode or ``FIRKLS_LMEM_NULL`` for a linear solver option before :c:func:`FIRKodeSetLinearSolver`.

   **Notes:**
      Each option is a key, ``<firkid>.<key>``, followed by its value(s) as
      separate strings. There is one key per scalar setter, named after it
      in snake case, e.g. ``firkode.max_num_steps 10000`` for
      :c:func:`FIRKodeSetMaxNumSteps`, ``firkode.scalar_tolerances 1e-6
      1e-10`` for :c:func:`FIRKodeSStolerances`, ``firkode.eps_lin 0.01``
      for :c:func:`FIRKodeSetEpsLin`. The keys are:

      - one integer: ``max_conv_fails``, ``max_err_test_fails``,
        ``max_hnil_warns``, ``max_nonlin_iters``, ``num_stages``, ``order``,
        ``predictor_method``, ``small_num_efails``, ``error_refilter``,
        ``interpolate_stop_time``, ``newton_count_safety``,
        ``linear_solution_scaling`` (the flags take 0 or 1),
        ``stage_solver_maxl``, ``stage_solver_max_restarts``;
      - one integer of type ``long int``: ``max_num_steps``,
        ``jac_eval_frequency``, ``lsetup_frequency``;
      - one real: ``init_step``, ``min_step``, ``max_step``, ``fixed_step``,
        ``stop_time``, ``error_bias``, ``safety_factor``, ``max_growth``,
        ``max_first_growth``, ``min_reduction``, ``max_efail_growth``,
        ``max_cfail_growth``, ``nonlin_conv_coef``,
        ``nonlin_divergence_rate``, ``nonlin_rate_exponent``,
        ``jac_reuse_conv_rate``, ``jac_bad_conv_rate``,
        ``delta_gamma_max_lsetup``, ``delta_gamma_max_bad_jac``, ``eps_lin``,
        ``ls_norm_factor``, ``stage_eps_lin``, ``mass_eps_lin``,
        ``mass_ls_norm_factor``;
      - two reals: ``scalar_tolerances``, ``fixed_step_bounds``;
      - a name: ``mode`` (``beginner`` or ``expert``, in either case),
        ``reuse_policy`` (``radau5`` or ``periodic``, in either case),
        ``stage_solver_gs_type`` (``modified`` or ``classical``),
        ``adapt_controller_by_name`` (a controller name, see
        :c:func:`FIRKodeSetAdaptControllerByName`);
      - no value: ``clear_stop_time``, ``no_inactive_root_warn``.

      The keys ``small_num_efails``, ``lsetup_frequency``,
      ``max_efail_growth``, ``max_cfail_growth`` and
      ``delta_gamma_max_lsetup`` are spelled as the corresponding ARKODE and
      CVODE keys rather than as their setters.

      A key of an EXPERT option is rejected in BEGINNER mode, as its setter
      is. The key ``mode`` is applied before every other key of the same
      call, wherever it stands, so that ``firkode.mode expert`` enables the
      EXPERT keys given with it. An unknown key is reported as a warning
      and skipped; strings without the prefix are ignored, and ``argv`` is
      not changed.

      .. code-block:: console

         $ ./a.out firkode.mode expert firkode.reuse_policy periodic firkode.max_num_steps 10000

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
      ``NULL`` pointer is passed. The pointer reaches every such function
      whenever this function is called, before or after the routines are
      attached and between calls to :c:func:`FIRKodeEvolve`: the right-hand
      side, error-weight, root, Jacobian, linear-system,
      Jacobian-times-vector, preconditioner, mass-matrix and mass
      preconditioner functions. The mass matrix-times-vector routine
      receives the data given to :c:func:`FIRKodeSetMassTimes` instead.

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
      As in RADAU5, an adaptive step is never shorter than
      :math:`10 u |t_n|`, :math:`u` the unit roundoff, whatever ``hmin``: a
      step that fails at that size ends the integration with
      ``FIRK_ERR_FAILURE`` or ``FIRK_CONV_FAILURE``.

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
      - ``FIRK_ILL_INPUT`` -- The mode is BEGINNER; this is an EXPERT option (see :c:func:`FIRKodeSetMode`).

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
      - ``FIRK_ILL_INPUT`` -- The mode is BEGINNER; this is an EXPERT option (see :c:func:`FIRKodeSetMode`).

   **Notes:**
      The default value is 10. A non-positive value restores the default.

.. c:function:: int FIRKodeSetReusePolicy(void* firkode_mem, int policy)

   The function :c:func:`FIRKodeSetReusePolicy` selects when the Jacobian is
   re-evaluated and the linear solver set up again (the matrix
   :math:`M - \gamma J` refactored or the preconditioner recomputed).

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``policy`` -- ``FIRK_POLICY_RADAU5`` or ``FIRK_POLICY_PERIODIC``.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.
      - ``FIRK_ILL_INPUT`` -- ``policy`` is not one of the two values.
      - ``FIRK_ILL_INPUT`` -- The mode is BEGINNER; this is an EXPERT option (see :c:func:`FIRKodeSetMode`).

   **Notes:**
      ``FIRK_POLICY_RADAU5``, the default, is the strategy of RADAU5: after an
      accepted step the Jacobian is re-evaluated unless the Newton iteration
      converged with a rate of at most :math:`\theta_{jac}` (see
      :c:func:`FIRKodeSetJacReuseConvRate`), after a failed attempt only if
      it was not evaluated at the current point, and the linear solver is set
      up whenever the Jacobian or :math:`\gamma` changed; a step size kept in
      the dead band therefore reuses the setup. This suits problems whose
      Jacobian and factorization are cheap relative to a Newton iteration.

      ``FIRK_POLICY_PERIODIC`` is the strategy of CVODE: the linear solver is
      set up on the first step, after a failed attempt, every ``msbp`` steps
      (:c:func:`FIRKodeSetLSetupFrequency`) and when :math:`\gamma` changed by
      more than ``dgmax_lsetup`` (:c:func:`FIRKodeSetDeltaGammaMaxLSetup`);
      the Jacobian on the first step, every ``msbj`` steps
      (:c:func:`FIRKodeSetJacEvalFrequency`) and after a convergence failure
      with a Jacobian that was not current (see
      :c:func:`FIRKodeSetDeltaGammaMaxBadJac`). This suits expensive setups,
      e.g., algebraic multigrid preconditioners. Those four options have no
      effect under ``FIRK_POLICY_RADAU5``.

      Under both policies a Jacobian is never evaluated twice at the same
      point, and a Newton rate above the one set with
      :c:func:`FIRKodeSetJacBadConvRate` forces a new Jacobian and setup.

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
      - ``FIRK_ILL_INPUT`` -- The mode is BEGINNER; this is an EXPERT option (see :c:func:`FIRKodeSetMode`).

   **Notes:**
      The default value is 20; passing 0 restores the default. Under
      ``FIRK_POLICY_PERIODIC`` (see :c:func:`FIRKodeSetReusePolicy`) the
      linear solver setup is additionally called on the first step, after a
      nonlinear solver or error test failure, and whenever :math:`\gamma`
      has changed by more than the fraction set with
      :c:func:`FIRKodeSetDeltaGammaMaxLSetup` (see
      :numref:`FIRKODE.Mathematics.linear`); the option has no effect under
      ``FIRK_POLICY_RADAU5``.

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
      - ``FIRK_ILL_INPUT`` -- The mode is BEGINNER; this is an EXPERT option (see :c:func:`FIRKodeSetMode`).

   **Notes:**
      The default value is 0.3. A negative value restores the default. The
      option is used by ``FIRK_POLICY_PERIODIC`` only (see
      :c:func:`FIRKodeSetReusePolicy`).

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
      - ``FIRK_ILL_INPUT`` -- The mode is BEGINNER; this is an EXPERT option (see :c:func:`FIRKodeSetMode`).

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
      - ``FIRK_ILL_INPUT`` -- The mode is BEGINNER; this is an EXPERT option (see :c:func:`FIRKodeSetMode`).


.. _FIRKODE.Usage.optional_input.adapt:

Time step adaptivity optional input functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following functions control the selection of the step size after each
step attempt (see :numref:`FIRKODE.Mathematics.adaptivity`). None of them
has any effect in fixed-step mode.

.. _FIRKODE.Usage.optional_input.adapt_table:

.. table:: Optional time step adaptivity inputs for FIRKODE

   +-------------------------------+----------------------------------------------+----------------+----------+
   |      **Optional input**       |              **Function name**               |  **Default**   | **Mode** |
   +===============================+==============================================+================+==========+
   | Time step controller object   | :c:func:`FIRKodeSetAdaptController`          | ImpGus         | E        |
   +-------------------------------+----------------------------------------------+----------------+----------+
   | Time step controller by name  | :c:func:`FIRKodeSetAdaptControllerByName`    | ImpGus         | E        |
   +-------------------------------+----------------------------------------------+----------------+----------+
   | Safety factor                 | :c:func:`FIRKodeSetSafetyFactor`             | 0.9            | E        |
   +-------------------------------+----------------------------------------------+----------------+----------+
   | Error bias                    | :c:func:`FIRKodeSetErrorBias`                | controller     | E        |
   +-------------------------------+----------------------------------------------+----------------+----------+
   | Maximum step growth           | :c:func:`FIRKodeSetMaxGrowth`                | 8.0            | E        |
   +-------------------------------+----------------------------------------------+----------------+----------+
   | Minimum step reduction        | :c:func:`FIRKodeSetMinReduction`             | 0.2            | E        |
   +-------------------------------+----------------------------------------------+----------------+----------+
   | Step size dead band           | :c:func:`FIRKodeSetFixedStepBounds`          | [1.0, 1.2]     | E        |
   +-------------------------------+----------------------------------------------+----------------+----------+
   | Maximum first step growth     | :c:func:`FIRKodeSetMaxFirstGrowth`           | 10000.0        | E        |
   +-------------------------------+----------------------------------------------+----------------+----------+
   | Maximum growth after error    | :c:func:`FIRKodeSetMaxEFailGrowth`           | 0.3            | E        |
   | test failures                 |                                              |                |          |
   +-------------------------------+----------------------------------------------+----------------+----------+
   | No. of error test failures    | :c:func:`FIRKodeSetSmallNumEFails`           | 2              | E        |
   | before limiting growth        |                                              |                |          |
   +-------------------------------+----------------------------------------------+----------------+----------+
   | Step reduction after a        | :c:func:`FIRKodeSetMaxCFailGrowth`           | 0.25           | E        |
   | convergence failure           |                                              |                |          |
   +-------------------------------+----------------------------------------------+----------------+----------+
   | Newton iteration count        | :c:func:`FIRKodeSetNewtonCountSafety`        | ``SUNTRUE``    | E        |
   | safety adjustment             |                                              |                |          |
   +-------------------------------+----------------------------------------------+----------------+----------+

.. c:function:: int FIRKodeSetAdaptController(void* firkode_mem, SUNAdaptController C)

   Sets a user-supplied time step controller object.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``C`` -- user-supplied time adaptivity controller. If ``C`` is ``NULL`` then the default controller is created.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.
      - ``FIRK_MEM_FAIL`` -- ``C`` was ``NULL`` and the default controller could not be allocated, or the previous controller could not be destroyed.
      - ``FIRK_CONTROLLER_ERR`` -- ``C`` was ``NULL`` and the parameters of the default controller could not be set.
      - ``FIRK_ILL_INPUT`` -- The controller is not of type ``SUN_ADAPTCONTROLLER_H``.
      - ``FIRK_ILL_INPUT`` -- The mode is BEGINNER; this is an EXPERT option (see :c:func:`FIRKodeSetMode`).

   **Notes:**
      The default controller is the implicit Gustafsson controller
      :c:func:`SUNAdaptController_ImpGus` (see
      :numref:`SUNAdaptController.Soderlind`) with both parameters equal to 1
      (see :c:func:`SUNAdaptController_SetParams_ImpGus`), which is the
      predictive controller of RADAU5; the I controller of
      :numref:`SUNAdaptController.Soderlind` reproduces the classical
      (non-predictive) RADAU5 step size selection. The controller is
      supplied with the order :math:`p = s` of the error estimate. The
      user-supplied controller is not freed by FIRKODE.

.. c:function:: int FIRKodeSetAdaptControllerByName(void* firkode_mem, const char* cname)

   Creates a time step controller of the Soderlind family (see
   :numref:`SUNAdaptController.Soderlind`) by its name and makes it the
   controller of FIRKODE, which frees it with the integrator or when it is
   replaced.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``cname`` -- the name of the controller: ``"Soderlind"``, ``"PID"``, ``"PI"``, ``"I"``, ``"ExpGus"``, ``"ImpGus"``, ``"H0211"``, ``"H0321"``, ``"H211"`` or ``"H312"``, which selects the constructor ``SUNAdaptController_<cname>``.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The controller was created and set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.
      - ``FIRK_MEM_FAIL`` -- The controller could not be allocated.
      - ``FIRK_ILL_INPUT`` -- ``cname`` is ``NULL`` or not one of the names above; the current controller is kept.
      - ``FIRK_ILL_INPUT`` -- The mode is BEGINNER; this is an EXPERT option (see :c:func:`FIRKodeSetMode`).

   **Notes:**
      The names are those of ``ARKodeSetAdaptControllerByName`` except
      ``"ImExGus"``, a controller for IMEX methods. The controller is created
      with the default parameters of its constructor; ``"ImpGus"`` thus
      differs from FIRKODE's default controller, whose parameters are 1.

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
      - ``FIRK_ILL_INPUT`` -- The mode is BEGINNER; this is an EXPERT option (see :c:func:`FIRKodeSetMode`).

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
      - ``FIRK_ILL_INPUT`` -- The mode is BEGINNER; this is an EXPERT option (see :c:func:`FIRKodeSetMode`).

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
      - ``FIRK_ILL_INPUT`` -- The mode is BEGINNER; this is an EXPERT option (see :c:func:`FIRKodeSetMode`).

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
      - ``FIRK_ILL_INPUT`` -- The mode is BEGINNER; this is an EXPERT option (see :c:func:`FIRKodeSetMode`).

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
      - ``FIRK_ILL_INPUT`` -- The mode is BEGINNER; this is an EXPERT option (see :c:func:`FIRKodeSetMode`).

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
      - ``FIRK_ILL_INPUT`` -- The mode is BEGINNER; this is an EXPERT option (see :c:func:`FIRKodeSetMode`).

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
      - ``FIRK_ILL_INPUT`` -- The mode is BEGINNER; this is an EXPERT option (see :c:func:`FIRKodeSetMode`).

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
      - ``FIRK_ILL_INPUT`` -- The mode is BEGINNER; this is an EXPERT option (see :c:func:`FIRKodeSetMode`).

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
      - ``FIRK_ILL_INPUT`` -- The mode is BEGINNER; this is an EXPERT option (see :c:func:`FIRKodeSetMode`).

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
      - ``FIRK_ILL_INPUT`` -- The mode is BEGINNER; this is an EXPERT option (see :c:func:`FIRKodeSetMode`).


.. _FIRKODE.Usage.optional_input.nls:

Nonlinear solver optional input functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following functions control the simplified Newton iteration used to
solve the stage equations (see :numref:`FIRKODE.Mathematics.Newton`).

.. _FIRKODE.Usage.optional_input.nls_table:

.. table:: Optional nonlinear solver inputs for FIRKODE

   +-------------------------------+----------------------------------------------+----------------+----------+
   |      **Optional input**       |              **Function name**               |  **Default**   | **Mode** |
   +===============================+==============================================+================+==========+
   | Maximum no. of Newton         | :c:func:`FIRKodeSetMaxNonlinIters`           | 7              | E        |
   | iterations                    |                                              |                |          |
   +-------------------------------+----------------------------------------------+----------------+----------+
   | Nonlinear convergence         | :c:func:`FIRKodeSetNonlinConvCoef`           | automatic      | E        |
   | coefficient                   |                                              |                |          |
   +-------------------------------+----------------------------------------------+----------------+----------+
   | Newton divergence rate        | :c:func:`FIRKodeSetNonlinDivergenceRate`     | 0.99           | E        |
   +-------------------------------+----------------------------------------------+----------------+----------+
   | Convergence rate triggering a | :c:func:`FIRKodeSetJacBadConvRate`           | 1.0 (disabled) | E        |
   | Jacobian update               |                                              |                |          |
   +-------------------------------+----------------------------------------------+----------------+----------+
   | Exponent of the first         | :c:func:`FIRKodeSetNonlinRateExponent`       | 0.8            | E        |
   | iteration's rate factor       |                                              |                |          |
   +-------------------------------+----------------------------------------------+----------------+----------+
   | Initial rate of each Newton   | :c:func:`FIRKodeSetJacReuseConvRate`         | 0.001          | E        |
   | solve                         |                                              |                |          |
   +-------------------------------+----------------------------------------------+----------------+----------+

.. c:function:: int FIRKodeSetMaxNonlinIters(void* firkode_mem, int maxcor)

   The function :c:func:`FIRKodeSetMaxNonlinIters` specifies the maximum
   number of Newton iterations permitted per step.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``maxcor`` -- maximum number of Newton iterations per step (:math:`> 0`).

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.
      - ``FIRK_ILL_INPUT`` -- The mode is BEGINNER; this is an EXPERT option (see :c:func:`FIRKodeSetMode`).

   **Notes:**
      The default value is 7; a non-positive value restores the default.

.. c:function:: int FIRKodeSetNonlinConvCoef(void* firkode_mem, sunrealtype nlscoef)

   The function :c:func:`FIRKodeSetNonlinConvCoef` specifies the safety
   factor :math:`\epsilon_{nls}` used in the nonlinear convergence test
   :eq:`FIRKODE_nlstest`.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nlscoef`` -- coefficient in the nonlinear convergence test (:math:`> 0`), or a non-positive value for the automatic default.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.
      - ``FIRK_ILL_INPUT`` -- The mode is BEGINNER; this is an EXPERT option (see :c:func:`FIRKodeSetMode`).

   **Notes:**
      By default :math:`\epsilon_{nls}` follows the relative tolerance as in
      RADAU5, :math:`\max(10 U / \mathrm{rtol}, \min(0.03, \sqrt{\mathrm{rtol}}))`
      with the unit roundoff :math:`U`, recomputed whenever the tolerances are
      set, and is 0.03 before the tolerances are set, for a zero relative
      tolerance and with a user weight function;
      :c:func:`FIRKodeWriteParameters` marks this automatic value. A positive
      ``nlscoef`` fixes the coefficient and a non-positive value restores the
      automatic default. The value also scales the tolerance of the stacked
      stage solve (see :c:func:`FIRKodeSetStageEpsLin`) and of an iterative
      block solver's Newton solves (see :c:func:`FIRKodeSetEpsLin`), not that
      of the error estimate.

.. c:function:: int FIRKodeSetNonlinDivergenceRate(void* firkode_mem, sunrealtype theta_max)

   The function :c:func:`FIRKodeSetNonlinDivergenceRate` specifies the
   contraction rate above which the Newton iteration is declared divergent.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``theta_max`` -- the divergence threshold (:math:`> 0`).

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.
      - ``FIRK_ILL_INPUT`` -- The mode is BEGINNER; this is an EXPERT option (see :c:func:`FIRKodeSetMode`).

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
      - ``FIRK_ILL_INPUT`` -- The mode is BEGINNER; this is an EXPERT option (see :c:func:`FIRKodeSetMode`).

   **Notes:**
      The test is disabled by default (the value 1.0); a non-positive value
      restores the default. A threshold of, e.g., 0.1 mimics the strategy of
      RADAU5 for problems with rapidly varying Jacobians.

.. c:function:: int FIRKodeSetNonlinRateExponent(void* firkode_mem, sunrealtype rate_exp)

   The function :c:func:`FIRKodeSetNonlinRateExponent` specifies the exponent
   :math:`p` of the factor :math:`\eta_1 = \max(\eta, U)^p` that replaces
   :math:`\theta/(1-\theta)` in the convergence test of the first Newton
   iteration, where :math:`\eta` is the last factor of the previous solve
   (see :numref:`FIRKODE.Mathematics.Newton`).

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``rate_exp`` -- the exponent :math:`p`; a non-positive value restores the default 0.8 of RADAU5.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.
      - ``FIRK_ILL_INPUT`` -- The mode is BEGINNER; this is an EXPERT option (see :c:func:`FIRKodeSetMode`).

   **Notes:**
      A smaller exponent makes the first iteration's test stricter.

.. c:function:: int FIRKodeSetJacReuseConvRate(void* firkode_mem, sunrealtype theta_jac)

   The function :c:func:`FIRKodeSetJacReuseConvRate` specifies the rate
   :math:`\theta_{jac}` of the Newton iteration below which the Jacobian
   counts as reusable, RADAU5's ``THET``; each Newton solve starts from this
   rate, which is the one reported by
   :c:func:`FIRKodeGetLastNonlinConvRate` after a solve that converged in one
   iteration.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``theta_jac`` -- the rate, less than 1; a non-positive value restores the default 0.001.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.
      - ``FIRK_ILL_INPUT`` -- ``theta_jac`` is 1 or larger.
      - ``FIRK_ILL_INPUT`` -- The mode is BEGINNER; this is an EXPERT option (see :c:func:`FIRKodeSetMode`).


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

   +-------------------------------------+----------------------------------------------+----------------+----------+
   |      **Optional input**             |              **Function name**               |  **Default**   | **Mode** |
   +=====================================+==============================================+================+==========+
   | Jacobian function                   | :c:func:`FIRKodeSetJacFn`                    | DQ             | B        |
   +-------------------------------------+----------------------------------------------+----------------+----------+
   | Linear system function              | :c:func:`FIRKodeSetLinSysFn`                 | internal       | B        |
   +-------------------------------------+----------------------------------------------+----------------+----------+
   | Jacobian evaluation frequency       | :c:func:`FIRKodeSetJacEvalFrequency`         | 51             | E        |
   +-------------------------------------+----------------------------------------------+----------------+----------+
   | Max change in :math:`\gamma`        | :c:func:`FIRKodeSetDeltaGammaMaxBadJac`      | 0.2            | E        |
   | attributed to a bad Jacobian        |                                              |                |          |
   +-------------------------------------+----------------------------------------------+----------------+----------+
   | Block linear solve tolerance factor | :c:func:`FIRKodeSetEpsLin`                   | 0.05           | B        |
   +-------------------------------------+----------------------------------------------+----------------+----------+
   | Newton linear solve tolerance       | :c:func:`FIRKodeSetLSNormFactor`             | sqrt(length)   | B        |
   | conversion factor                   |                                              |                |          |
   +-------------------------------------+----------------------------------------------+----------------+----------+
   | Linear solution scaling             | :c:func:`FIRKodeSetLinearSolutionScaling`    | on (matrix)    | B        |
   +-------------------------------------+----------------------------------------------+----------------+----------+
   | Preconditioner functions            | :c:func:`FIRKodeSetPreconditioner`           | ``NULL``       | B        |
   +-------------------------------------+----------------------------------------------+----------------+----------+
   | Jacobian-times-vector functions     | :c:func:`FIRKodeSetJacTimes`                 | DQ             | B        |
   +-------------------------------------+----------------------------------------------+----------------+----------+
   | Jacobian-times-vector DQ RHS        | :c:func:`FIRKodeSetJacTimesRhsFn`            | ``NULL``       | B        |
   | function                            |                                              |                |          |
   +-------------------------------------+----------------------------------------------+----------------+----------+
   | Stage Krylov subspace dimension     | :c:func:`FIRKodeSetStageSolverMaxl`          | max(5, 3s)     | E        |
   +-------------------------------------+----------------------------------------------+----------------+----------+
   | Stage Krylov restarts               | :c:func:`FIRKodeSetStageSolverMaxRestarts`   | 1              | E        |
   +-------------------------------------+----------------------------------------------+----------------+----------+
   | Stage Krylov orthogonalization      | :c:func:`FIRKodeSetStageSolverGSType`        | modified GS    | E        |
   +-------------------------------------+----------------------------------------------+----------------+----------+
   | Stage Krylov tolerance factor       | :c:func:`FIRKodeSetStageEpsLin`              | 0.05           | E        |
   +-------------------------------------+----------------------------------------------+----------------+----------+

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
      - ``FIRKLS_ILL_INPUT`` -- A Jacobian routine was supplied for a linear solver without a ``SUNMatrix``. (Without a Jacobian routine, a matrix type for which no difference-quotient approximation exists is reported by :c:func:`FIRKodeEvolve` as ``FIRK_LINIT_FAIL``.)

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
      - ``FIRKLS_MEM_FAIL`` -- The copy of the Jacobian, saved again when the default is restored, could not be allocated.

   **Notes:**
      By default, FIRKLS uses an internal linear system function that
      leverages the ``SUNMatrix`` API to form the system :math:`M - \gamma J`.
      If ``NULL`` is passed to ``linsys``, this default function is used. The
      function type :c:type:`FIRKLsLinSysFn` is described in
      :numref:`FIRKODE.Usage.user_supplied`.

      Only the default function saves the Jacobian, which the internal
      Krylov iteration on the coupled stage system uses for its
      matrix-vector products; with a user function these come from a
      Jacobian-times-vector routine or difference quotients. This function
      may be called between calls to :c:func:`FIRKodeEvolve`: setting a user
      function drops the saved Jacobian, and setting ``NULL`` saves one again
      from the next Jacobian evaluation.

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
      - ``FIRK_ILL_INPUT`` -- The mode is BEGINNER; this is an EXPERT option (see :c:func:`FIRKodeSetMode`).

   **Notes:**
      The default value is 51; passing 0 restores the default. If the
      Jacobian or preconditioner data is out-of-date when a linear solver
      setup is performed, it will be recomputed; otherwise the saved data is
      reused with the current :math:`\gamma` (see
      :numref:`FIRKODE.Mathematics.linear`). The option is used by
      ``FIRK_POLICY_PERIODIC`` only (see :c:func:`FIRKodeSetReusePolicy`).

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
      - ``FIRK_ILL_INPUT`` -- The mode is BEGINNER; this is an EXPERT option (see :c:func:`FIRKodeSetMode`).

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
      - ``FIRKLS_ILL_INPUT`` -- ``eplifac`` is negative.

   **Notes:**
      The default value is 0.05; passing a zero value restores the
      default.

      With a mass matrix the tolerance is applied to the residual converted
      to the units of :math:`y`, see :numref:`FIRKODE.Mathematics.mass`.

.. c:function:: int FIRKodeSetLSNormFactor(void* firkode_mem, sunrealtype nrmfac)

   The function :c:func:`FIRKodeSetLSNormFactor` specifies the factor to use
   when converting from the integrator tolerance (WRMS norm) to the linear
   solver tolerance (L2 norm) for the block systems. The tolerance of the
   internal Krylov iteration on the coupled stage system of :math:`s` stages
   uses :math:`\sqrt{s}` times this factor.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nrmfac`` -- the norm conversion factor. If ``nrmfac`` is :math:`> 0` then the provided value is used; if ``nrmfac`` is :math:`= 0` then the conversion factor is computed as the square root of the vector length (the default); if ``nrmfac`` is :math:`< 0` then the conversion factor is computed using the vector dot product :math:`\sqrt{v^T v}` for a vector of ones.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

.. c:function:: int FIRKodeSetLinearSolutionScaling(void* firkode_mem, sunbooleantype onoff)

   The function :c:func:`FIRKodeSetLinearSolutionScaling` enables or disables
   scaling the Newton correction computed by a matrix-based linear solver by
   :math:`2/(1 + \gamma/\gamma_p)`, where :math:`\gamma_p` is the value of
   :math:`\gamma` at the last setup of the linear system, as in CVODE. The
   scaling applies only when the block solver solves the Newton system
   itself, i.e., for the one-stage method, which this version does not offer;
   with :math:`s = 3` stages the block solver preconditions the coupled stage
   system and its solutions are never scaled.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``onoff`` -- flag to enable (``SUNTRUE``) or disable (``SUNFALSE``) the scaling.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The flag value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.
      - ``FIRKLS_ILL_INPUT`` -- The attached linear solver is not matrix-based.

   **Notes:**
      The scaling is enabled by default when a matrix-based linear solver is
      attached. This function must be called after
      :c:func:`FIRKodeSetLinearSolver`.

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
      setup function ``jtsetup``, if supplied, is called once per Newton
      iteration of that stage system, and once per error estimate with a
      matrix-free linear solver, also when the products come from the
      saved matrix and ``jtimes`` is not called. The function types
      :c:type:`FIRKLsJacTimesSetupFn` and :c:type:`FIRKLsJacTimesVecFn` are
      described in :numref:`FIRKODE.Usage.user_supplied`.

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
      - ``FIRK_ILL_INPUT`` -- The mode is BEGINNER; this is an EXPERT option (see :c:func:`FIRKodeSetMode`).

   **Notes:**
      The default is :math:`\max(5, 3s)` for an :math:`s`-stage method. The
      internal solver stores :math:`2(\mathrm{maxl}+1)` vectors of the size
      of the coupled stage system, i.e., :math:`2 s (\mathrm{maxl}+1)`
      vectors of the size of the ODE system, so this value controls the
      memory footprint of FIRKODE. It has no effect for the one-stage method.
      This function may be called between calls to :c:func:`FIRKodeEvolve`;
      the internal solver is created again with the new dimension at its
      next use.

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
      - ``FIRKLS_SUNLS_FAIL`` -- The stacked solver rejected the value.
      - ``FIRK_ILL_INPUT`` -- The mode is BEGINNER; this is an EXPERT option (see :c:func:`FIRKodeSetMode`).

   **Notes:**
      The default value is 1.

.. c:function:: int FIRKodeSetStageSolverGSType(void* firkode_mem, int gstype)

   The function :c:func:`FIRKodeSetStageSolverGSType` specifies the
   Gram-Schmidt orthogonalization of the internal flexible GMRES iteration on
   the coupled stage system.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``gstype`` -- ``SUN_MODIFIED_GS`` (the default) or ``SUN_CLASSICAL_GS``.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.
      - ``FIRKLS_ILL_INPUT`` -- ``gstype`` is not one of the two values.
      - ``FIRKLS_SUNLS_FAIL`` -- The internal solver rejected the value.
      - ``FIRK_ILL_INPUT`` -- The mode is BEGINNER; this is an EXPERT option (see :c:func:`FIRKodeSetMode`).

   **Notes:**
      Classical Gram-Schmidt needs fewer global reductions per iteration,
      which can matter at scale, at the price of a less stable
      orthogonalization. It has no effect for the one-stage method. This
      function may be called between calls to :c:func:`FIRKodeEvolve`.

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
      - ``FIRKLS_ILL_INPUT`` -- ``eplifac`` is negative.
      - ``FIRK_ILL_INPUT`` -- The mode is BEGINNER; this is an EXPERT option (see :c:func:`FIRKodeSetMode`).

   **Notes:**
      The default value is 0.05; passing a zero value restores the
      default.

      With a mass matrix the tolerance is applied to the residual converted
      to the units of :math:`y`, see :numref:`FIRKODE.Mathematics.mass`.


.. _FIRKODE.Usage.optional_input.mass:

Mass matrix solver interface optional input functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All of these functions must be called *after* the mass matrix solver
interface has been initialized with :c:func:`FIRKodeSetMassLinearSolver`.

.. _FIRKODE.Usage.optional_input.mass_table:

.. table:: Optional inputs for the FIRKLS mass matrix solver interface

   +-------------------------------------+----------------------------------------------+----------------+----------+
   |      **Optional input**             |              **Function name**               |  **Default**   | **Mode** |
   +=====================================+==============================================+================+==========+
   | Mass matrix function                | :c:func:`FIRKodeSetMassFn`                   | none           | B        |
   +-------------------------------------+----------------------------------------------+----------------+----------+
   | Mass matrix-times-vector functions  | :c:func:`FIRKodeSetMassTimes`                | none           | B        |
   +-------------------------------------+----------------------------------------------+----------------+----------+
   | Mass matrix preconditioner          | :c:func:`FIRKodeSetMassPreconditioner`       | ``NULL``       | B        |
   | functions                           |                                              |                |          |
   +-------------------------------------+----------------------------------------------+----------------+----------+
   | Mass matrix linear solve tolerance  | :c:func:`FIRKodeSetMassEpsLin`               | 0.05           | B        |
   | factor                              |                                              |                |          |
   +-------------------------------------+----------------------------------------------+----------------+----------+
   | Mass matrix linear solve tolerance  | :c:func:`FIRKodeSetMassLSNormFactor`         | sqrt(length)   | B        |
   | conversion factor                   |                                              |                |          |
   +-------------------------------------+----------------------------------------------+----------------+----------+

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

.. c:function:: int FIRKodeSetMassTimes(void* firkode_mem, FIRKLsMassTimesSetupFn msetup, FIRKLsMassTimesVecFn mtimes, void* mtimes_data)

   The function :c:func:`FIRKodeSetMassTimes` specifies the mass
   matrix-vector product setup and multiply functions. This function is
   required when a ``NULL`` ``SUNMatrix`` was passed to
   :c:func:`FIRKodeSetMassLinearSolver`; it is optional otherwise (in which
   case the product is computed with :c:func:`SUNMatMatvec` by default).

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``msetup`` -- user-defined mass matrix-vector product setup function of type :c:type:`FIRKLsMassTimesSetupFn`, or ``NULL`` if no setup is needed.
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
   the WRMS norm of the right-hand side of a mass matrix linear system is
   multiplied to obtain the residual tolerance of an iterative mass matrix
   linear solver.

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

      The tolerance is relative to the right-hand side, so that the mass
      solves have the same relative accuracy whatever the scaling of
      :math:`M` and :math:`f`, both in the solves that convert the
      right-hand sides of the iterative linear solves to the units of
      :math:`y` and in those for :math:`\dot{y} = M^{-1} f(t, y)` in the
      initial step size estimate and in rootfinding (see
      :numref:`FIRKODE.Mathematics.mass`). Since :math:`\dot{y}` only feeds
      these heuristics, an approximate solution whose residual the solver
      reduced without reaching the tolerance is accepted for it.

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

   +-------------------------------------+----------------------------------------------+----------------+----------+
   |      **Optional input**             |              **Function name**               |  **Default**   | **Mode** |
   +=====================================+==============================================+================+==========+
   | Direction of zero-crossing          | :c:func:`FIRKodeSetRootDirection`            | both           | B        |
   +-------------------------------------+----------------------------------------------+----------------+----------+
   | Disable rootfinding warnings        | :c:func:`FIRKodeSetNoInactiveRootWarn`       | none           | B        |
   +-------------------------------------+----------------------------------------------+----------------+----------+

.. c:function:: int FIRKodeSetRootDirection(void* firkode_mem, int* rootdir_1d)

   The function :c:func:`FIRKodeSetRootDirection` specifies the direction of
   zero-crossings to be located and returned.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``rootdir_1d`` -- state array of length ``nrtfn``, the number of root functions :math:`g_i`, as specified in the call to the function :c:func:`FIRKodeRootInit`. A value of 0 for ``rootdir_1d[i]`` indicates that crossing in either direction for :math:`g_i` should be reported. A value of +1 or -1 indicates that the solver should report only zero-crossings where :math:`g_i` is increasing or decreasing, respectively.

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

.. c:function:: int FIRKodeGetNumSteps(void* firkode_mem, long int* nsteps_ptr)

   Returns the cumulative number of internal steps taken by the solver (total
   so far).

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nsteps_ptr`` -- number of steps taken by FIRKODE.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

.. c:function:: int FIRKodeGetNumStepAttempts(void* firkode_mem, long int* step_attempts_ptr)

   Returns the cumulative number of steps attempted by the solver (total so
   far), including rejected steps.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``step_attempts_ptr`` -- number of steps attempted by FIRKODE.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

.. c:function:: int FIRKodeGetNumRhsEvals(void* firkode_mem, long int* nfevals_ptr)

   Returns the number of calls to the user's right-hand side function.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nfevals_ptr`` -- number of calls to the user's ``f`` function.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

   **Notes:**
      The ``nfevals_ptr`` value returned by :c:func:`FIRKodeGetNumRhsEvals` does
      not account for calls made to ``f`` by a linear solver or preconditioner
      module (see :c:func:`FIRKodeGetNumLinRhsEvals`).

.. c:function:: int FIRKodeGetNumLinSolvSetups(void* firkode_mem, long int* nlinsetups_ptr)

   Returns the number of calls made to the linear solver's setup function.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nlinsetups_ptr`` -- number of calls made to the linear solver setup function.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

.. c:function:: int FIRKodeGetNumErrTestFails(void* firkode_mem, long int* netfails_ptr)

   Returns the number of local error test failures that have occurred.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``netfails_ptr`` -- number of error test failures.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

.. c:function:: int FIRKodeGetNumStepSolveFails(void* firkode_mem, long int* nncfails_ptr)

   Returns the number of failed steps due to a nonlinear solver failure.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nncfails_ptr`` -- number of step failures.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

.. c:function:: int FIRKodeGetActualInitStep(void* firkode_mem, sunrealtype* hinused_ptr)

   Returns the value of the integration step size used on the first step.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``hinused_ptr`` -- actual value of initial step size.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

   **Notes:**
      Even if the value of the initial integration step size was specified by
      the user through a call to :c:func:`FIRKodeSetInitStep`, this value
      might have been changed by FIRKODE to ensure that the step size is
      within the prescribed bounds (:math:`h_{min} \le h_0 \le h_{max}`), or
      to satisfy the local error test condition.

.. c:function:: int FIRKodeGetLastStep(void* firkode_mem, sunrealtype* hlast_ptr)

   Returns the integration step size taken on the last internal step.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``hlast_ptr`` -- step size taken on the last internal step.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

.. c:function:: int FIRKodeGetCurrentStep(void* firkode_mem, sunrealtype* hcur_ptr)

   Returns the integration step size to be attempted on the next internal
   step.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``hcur_ptr`` -- step size to be attempted on the next internal step.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

.. c:function:: int FIRKodeGetCurrentTime(void* firkode_mem, sunrealtype* tcur_ptr)

   Returns the current internal time reached by the solver.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``tcur_ptr`` -- current internal time reached.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

.. c:function:: int FIRKodeGetCurrentState(void* firkode_mem, N_Vector* y_ptr)

   Returns a pointer to the current state vector, i.e., the solution at the
   current internal time.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``y_ptr`` -- pointer to the current state vector.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

   **Notes:**
      The vector is owned by FIRKODE and must not be modified or freed by the
      user. The pointer may change between steps.

.. c:function:: int FIRKodeGetCurrentGamma(void* firkode_mem, sunrealtype* gamma_ptr)

   Returns the current value of the scalar :math:`\gamma` in the matrix
   :math:`M - \gamma J`.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``gamma_ptr`` -- the current value of :math:`\gamma`.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

.. c:function:: int FIRKodeGetTolScaleFactor(void* firkode_mem, sunrealtype* tolsfac_ptr)

   Returns a suggested factor by which the user's tolerances should be
   scaled when too much accuracy has been requested for some internal step.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``tolsfac_ptr`` -- suggested scaling factor for user-supplied tolerances.

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

.. c:function:: int FIRKodeGetIntegratorStats(void* firkode_mem, long int* nsteps_ptr, long int* nfevals_ptr, long int* nlinsetups_ptr, long int* netfails_ptr, sunrealtype* hinused_ptr, sunrealtype* hlast_ptr, sunrealtype* hcur_ptr, sunrealtype* tcur_ptr)

   Returns the main integrator statistics as a group.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nsteps_ptr`` -- number of steps taken by FIRKODE.
      - ``nfevals_ptr`` -- number of calls to the user's ``f`` function.
      - ``nlinsetups_ptr`` -- number of calls made to the linear solver setup function.
      - ``netfails_ptr`` -- number of error test failures.
      - ``hinused_ptr`` -- actual value of initial step size.
      - ``hlast_ptr`` -- step size taken on the last internal step.
      - ``hcur_ptr`` -- step size to be attempted on the next internal step.
      - ``tcur_ptr`` -- current internal time reached.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output values have been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

.. c:function:: int FIRKodeGetUserData(void* firkode_mem, void** user_data_ptr)

   Returns the user data pointer set with :c:func:`FIRKodeSetUserData`.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``user_data_ptr`` -- the user data pointer.

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

.. c:function:: int FIRKodeGetNumNonlinSolvIters(void* firkode_mem, long int* nniters_ptr)

   Returns the number of Newton iterations performed.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nniters_ptr`` -- number of Newton iterations.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

.. c:function:: int FIRKodeGetNumNonlinSolvConvFails(void* firkode_mem, long int* nnfails_ptr)

   Returns the number of nonlinear solver convergence failures.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nnfails_ptr`` -- number of nonlinear convergence failures.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

.. c:function:: int FIRKodeGetNonlinSolvStats(void* firkode_mem, long int* nniters_ptr, long int* nnfails_ptr)

   Returns the nonlinear solver statistics as a group.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nniters_ptr`` -- number of Newton iterations.
      - ``nnfails_ptr`` -- number of nonlinear convergence failures.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output values have been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.

.. c:function:: int FIRKodeGetLastNonlinConvRate(void* firkode_mem, sunrealtype* theta_ptr)

   Returns the rate of convergence :math:`\theta` of the Newton iteration at
   the end of the last nonlinear solve (see
   :numref:`FIRKODE.Mathematics.Newton`): the estimate of the last iteration
   that computed one, or the initial rate set with
   :c:func:`FIRKodeSetJacReuseConvRate` if the solve ended before the second
   iteration.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``theta_ptr`` -- the rate of convergence.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.


.. _FIRKODE.Usage.optional_output.ls:

Linear solver interface optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following optional outputs are available from the FIRKLS interface. In
the descriptions below, "block" linear solves refer to solves with the
matrices :math:`M - \gamma J` performed by the user's linear solver, and
"stage" linear iterations refer to the internal flexible GMRES iteration on
the coupled stage system.

.. c:function:: int FIRKodeGetJac(void* firkode_mem, SUNMatrix* J_ptr)

   Returns the internally stored copy of the Jacobian matrix of the ODE
   right-hand side function.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``J_ptr`` -- the Jacobian matrix.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

   **Notes:**
      This function is provided for debugging purposes and the values in the
      returned matrix should not be altered. The pointer is ``NULL`` for
      matrix-free linear solvers.

.. c:function:: int FIRKodeGetJacTime(void* firkode_mem, sunrealtype* t_J_ptr)

   Returns the time at which the internally stored copy of the Jacobian
   matrix was evaluated.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``t_J_ptr`` -- the time at which the Jacobian was evaluated.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

.. c:function:: int FIRKodeGetJacNumSteps(void* firkode_mem, long int* nst_J_ptr)

   Returns the value of the internal step counter at the time the internally
   stored copy of the Jacobian matrix was evaluated.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nst_J_ptr`` -- the value of the internal step counter.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

.. c:function:: int FIRKodeGetNumJacEvals(void* firkode_mem, long int* njevals_ptr)

   Returns the number of calls made to the Jacobian approximation function.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``njevals_ptr`` -- the number of calls to the Jacobian function.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

.. c:function:: int FIRKodeGetNumLinRhsEvals(void* firkode_mem, long int* nfevalsLS_ptr)

   Returns the number of calls made to the user-supplied right-hand side
   function due to the finite difference Jacobian approximation or finite
   difference Jacobian-vector product approximation.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nfevalsLS_ptr`` -- the number of calls made to the user-supplied right-hand side function.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

   **Notes:**
      The value ``nfevalsLS_ptr`` is incremented only if the default internal
      difference quotient function is used.

.. c:function:: int FIRKodeGetNumPrecEvals(void* firkode_mem, long int* npevals_ptr)

   Returns the number of preconditioner evaluations, i.e., the number of
   calls made to ``pset``, whatever their ``jok`` argument.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``npevals_ptr`` -- the current number of calls to ``pset``.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

.. c:function:: int FIRKodeGetNumPrecSolves(void* firkode_mem, long int* npsolves_ptr)

   Returns the cumulative number of calls made to the preconditioner solve
   function, ``psolve``.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``npsolves_ptr`` -- the current number of calls to ``psolve``.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

.. c:function:: int FIRKodeGetNumLinIters(void* firkode_mem, long int* nliters_ptr)

   Returns the cumulative number of linear iterations performed by the user's
   iterative linear solver in block solves.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nliters_ptr`` -- the current number of linear iterations.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

.. c:function:: int FIRKodeGetNumLinConvFails(void* firkode_mem, long int* nlcfails_ptr)

   Returns the cumulative number of convergence failures of the user's
   iterative linear solver in block solves.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nlcfails_ptr`` -- the current number of linear convergence failures.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

.. c:function:: int FIRKodeGetNumJTSetupEvals(void* firkode_mem, long int* njtsetups_ptr)

   Returns the cumulative number of calls made to the Jacobian-vector setup
   function, ``jtsetup``.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``njtsetups_ptr`` -- the current number of calls to ``jtsetup``.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

.. c:function:: int FIRKodeGetNumJtimesEvals(void* firkode_mem, long int* njvevals_ptr)

   Returns the cumulative number of calls made to the Jacobian-vector
   function, ``jtimes`` (user-supplied or difference quotient).

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``njvevals_ptr`` -- the current number of calls to ``jtimes``.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

.. c:function:: int FIRKodeGetNumBlockSolves(void* firkode_mem, long int* nbsolves_ptr)

   Returns the cumulative number of block linear solves, i.e., solves with
   the matrices :math:`M - \gamma J` performed by the user's linear solver
   (as preconditioner applications and for the error estimate).

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nbsolves_ptr`` -- the current number of block linear solves.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

.. c:function:: int FIRKodeGetNumStageLinIters(void* firkode_mem, long int* nsliters_ptr)

   Returns the cumulative number of iterations of the internal flexible GMRES
   solver on the coupled stage system.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nsliters_ptr`` -- the current number of stage linear iterations.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

.. c:function:: int FIRKodeGetNumStageLinConvFails(void* firkode_mem, long int* nslcfails_ptr)

   Returns the cumulative number of convergence failures of the internal
   flexible GMRES solver on the coupled stage system.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nslcfails_ptr`` -- the current number of stage linear convergence failures.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

.. c:function:: int FIRKodeGetLastLinFlag(void* firkode_mem, long int* flag_ptr)

   Returns the last return value from a FIRKLS routine.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``flag_ptr`` -- the value of the last return flag from a FIRKLS function.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_LMEM_NULL`` -- The FIRKLS linear solver interface has not been initialized.

   **Notes:**
      If the FIRKLS setup function failed when using the
      :ref:`SUNLINSOL_DENSE <SUNLinSol_Dense>` or :ref:`SUNLINSOL_BAND
      <SUNLinSol_Band>` modules, then the value of ``flag_ptr`` is equal to the
      column index (numbered from one) at which a zero diagonal element was
      encountered during the LU factorization of the (dense or banded)
      Jacobian matrix. For all other failures, ``flag_ptr`` is negative.

.. c:function:: char* FIRKodeGetLinReturnFlagName(long int flag)

   Returns the name of the FIRKLS constant corresponding to ``flag``, or of
   the FIRKODE constant for a negative ``flag`` that is not a FIRKLS
   constant, such as ``FIRK_ILL_INPUT``, which the setters of EXPERT options
   return in BEGINNER mode.

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

.. c:function:: int FIRKodeGetCurrentMassMatrix(void* firkode_mem, SUNMatrix* M_ptr)

   Returns the mass matrix object supplied to :c:func:`FIRKodeSetMassLinearSolver`.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``M_ptr`` -- the mass matrix (``NULL`` for a matrix-free mass linear solver).

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_MASSMEM_NULL`` -- The mass matrix solver interface has not been initialized.

.. c:function:: int FIRKodeGetNumMassSetups(void* firkode_mem, long int* nmsetups_ptr)

   Returns the number of calls made to the mass matrix linear solver setup
   function.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nmsetups_ptr`` -- number of calls to the mass matrix solver setup function.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_MASSMEM_NULL`` -- The mass matrix solver interface has not been initialized.

.. c:function:: int FIRKodeGetNumMassMultSetups(void* firkode_mem, long int* nmvsetups_ptr)

   Returns the number of calls made to the mass matrix-vector product setup
   routines (the user-supplied ``mtsetup`` routine and the
   :c:func:`SUNMatMatvecSetup` routine of the mass matrix).

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nmvsetups_ptr`` -- number of calls to the mass matrix-vector product setup routines.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_MASSMEM_NULL`` -- The mass matrix solver interface has not been initialized.

.. c:function:: int FIRKodeGetNumMassMult(void* firkode_mem, long int* nmvevals_ptr)

   Returns the number of mass matrix-vector products computed (with the
   user-supplied ``mtimes`` routine or with :c:func:`SUNMatMatvec`).

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nmvevals_ptr`` -- number of mass matrix-vector products.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_MASSMEM_NULL`` -- The mass matrix solver interface has not been initialized.

.. c:function:: int FIRKodeGetNumMassSolves(void* firkode_mem, long int* nmsolves_ptr)

   Returns the number of linear solves with the mass matrix.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nmsolves_ptr`` -- number of mass matrix linear solves.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_MASSMEM_NULL`` -- The mass matrix solver interface has not been initialized.

.. c:function:: int FIRKodeGetNumMassPrecEvals(void* firkode_mem, long int* nmpevals_ptr)

   Returns the number of calls made to the mass matrix preconditioner setup
   function.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nmpevals_ptr`` -- number of calls to the mass matrix preconditioner setup function.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_MASSMEM_NULL`` -- The mass matrix solver interface has not been initialized.

.. c:function:: int FIRKodeGetNumMassPrecSolves(void* firkode_mem, long int* nmpsolves_ptr)

   Returns the number of calls made to the mass matrix preconditioner solve
   function.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nmpsolves_ptr`` -- number of calls to the mass matrix preconditioner solve function.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_MASSMEM_NULL`` -- The mass matrix solver interface has not been initialized.

.. c:function:: int FIRKodeGetNumMassIters(void* firkode_mem, long int* nmiters_ptr)

   Returns the cumulative number of linear iterations performed by an
   iterative mass matrix linear solver.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nmiters_ptr`` -- number of mass matrix linear iterations.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_MASSMEM_NULL`` -- The mass matrix solver interface has not been initialized.

.. c:function:: int FIRKodeGetNumMassConvFails(void* firkode_mem, long int* nmcfails_ptr)

   Returns the cumulative number of convergence failures of the mass matrix
   linear solver.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``nmcfails_ptr`` -- number of mass matrix linear solver convergence failures.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_MASSMEM_NULL`` -- The mass matrix solver interface has not been initialized.

.. c:function:: int FIRKodeGetLastMassFlag(void* firkode_mem, long int* flag_ptr)

   Returns the last return value from a mass matrix solver interface routine.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``flag_ptr`` -- the value of the last return flag.

   **Return value:**
      - ``FIRKLS_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRKLS_MEM_NULL`` -- The ``firkode_mem`` pointer is ``NULL``.
      - ``FIRKLS_MASSMEM_NULL`` -- The mass matrix solver interface has not been initialized.


.. _FIRKODE.Usage.optional_output.root:

Rootfinding optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: int FIRKodeGetRootInfo(void* firkode_mem, int* rootsfound_1d)

   Returns an array showing which functions were found to have a root.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``rootsfound_1d`` -- array of length ``nrtfn`` with the indices of the user functions :math:`g_i` found to have a root. For :math:`i = 0, \ldots,` ``nrtfn`` :math:`-1`, ``rootsfound_1d[i]`` is nonzero if :math:`g_i` has a root, and 0 if not.

   **Return value:**
      - ``FIRK_SUCCESS`` -- The optional output value has been successfully set.
      - ``FIRK_MEM_NULL`` -- The FIRKODE memory block was not initialized.
      - ``FIRK_ILL_INPUT`` -- Rootfinding has not been activated through a call to :c:func:`FIRKodeRootInit`.

   **Notes:**
      Note that, for the components :math:`g_i` for which a root was found,
      the sign of ``rootsfound[i]`` indicates the direction of
      zero-crossing. A value of +1 indicates that :math:`g_i` is increasing,
      while a value of -1 indicates a decreasing :math:`g_i`. The user must
      allocate memory for the vector ``rootsfound_1d``.

.. c:function:: int FIRKodeGetNumGEvals(void* firkode_mem, long int* ngevals_ptr)

   Returns the cumulative number of calls made to the user's root function
   :math:`g`.

   **Arguments:**
      - ``firkode_mem`` -- pointer to the FIRKODE memory block.
      - ``ngevals_ptr`` -- number of calls made to the user's function :math:`g` thus far (zero if rootfinding is not active).

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
      All counters, of the integrator and of the linear, stage and mass
      solvers and the rootfinding, are reset to zero, and the mass matrix is
      set up anew at the next call to :c:func:`FIRKodeEvolve`.

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
      with :c:func:`FIRKodeSetInitStep`. All counters, of the integrator and
      of the linear, stage and mass solvers and the rootfinding, are
      retained, the mass matrix is not set up again, and the history used by
      the stage predictor is discarded.
