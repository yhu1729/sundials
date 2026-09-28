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

.. _FIRKODE.Usage.user_supplied:

User-supplied functions
=======================

The user-supplied functions consist of one function defining the ODE
right-hand side, (optionally) a function that provides the error weight
vector, (optionally) a function that provides the root functions,
(optionally) one or more functions that provide Jacobian-related information
for the linear solver, (optionally) one or two functions that define the
preconditioner for use in any of the Krylov iterative algorithms, and
(optionally) functions that define the mass matrix.

.. _FIRKODE.Usage.user_supplied.rhs:

ODE right-hand side
-------------------

The user must provide a function of type defined as follows:

.. c:type:: int (*FIRKRhsFn)(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)

   This function computes the ODE right-hand side :math:`f(t,y)` of the
   system :math:`M \dot{y} = f(t,y)` for a given value of the independent
   variable :math:`t` and state vector :math:`y`.

   **Arguments:**
      - ``t`` -- is the current value of the independent variable.
      - ``y`` -- is the current value of the dependent variable vector, :math:`y(t)`.
      - ``ydot`` -- is the output vector :math:`f(t,y)`.
      - ``user_data`` -- is the ``user_data`` pointer passed to :c:func:`FIRKodeSetUserData`.

   **Return value:**
      A :c:type:`FIRKRhsFn` should return 0 if successful, a positive value if
      a recoverable error occurred (in which case FIRKODE will attempt to
      correct), or a negative value if it failed unrecoverably (in which case
      the integration is halted and ``FIRK_RHSFUNC_FAIL`` is returned).

   **Notes:**
      Allocation of memory for ``ydot`` is handled within FIRKODE.

      A recoverable failure error return from the :c:type:`FIRKRhsFn` is
      typically used to flag a value of the dependent variable :math:`y` that
      is "illegal" in some way (e.g., negative where only a non-negative value
      is physically meaningful). If such a return is made, FIRKODE will
      attempt to recover (possibly repeating the nonlinear solve, or reducing
      the step size) in order to avoid this recoverable error return.

      The right-hand side is evaluated at the stage values :math:`y_n + Z_i`
      of each Newton iteration and, at the beginning of each step, at the
      accepted solution :math:`(t_n, y_n)`. There are two situations in which
      recovery is not possible even if the right-hand side function returns a
      recoverable error flag: at the very first call to the
      :c:type:`FIRKRhsFn` (in which case FIRKODE returns
      ``FIRK_FIRST_RHSFUNC_ERR``), and when the evaluation at the beginning of
      a step fails recoverably (in which case FIRKODE returns
      ``FIRK_UNREC_RHSFUNC_ERR``).

      When a mass matrix is supplied, the function must return :math:`f(t,y)`
      and not :math:`M^{-1} f(t,y)`.

.. _FIRKODE.Usage.user_supplied.ewt:

Error weight function
---------------------

As an alternative to providing the relative and absolute tolerances, the
user may provide a function of type :c:type:`FIRKEwtFn` to compute a vector
``ewt`` containing the weights in the WRMS norm
:math:`\|v\|_{\text{WRMS}} = \left( \frac{1}{N} \sum_{i=1}^N (W_i\, v_i)^2
\right)^{1/2}`. These weights will be used in place of those defined by
:eq:`FIRKODE_errwt`. The function type is defined as follows:

.. c:type:: int (*FIRKEwtFn)(N_Vector y, N_Vector ewt, void* user_data)

   This function computes the WRMS error weights for the vector :math:`y`.

   **Arguments:**
      - ``y`` -- the value of the dependent variable vector at which the weight vector is to be computed.
      - ``ewt`` -- the output vector containing the error weights.
      - ``user_data`` -- a pointer to user data, the same as the ``user_data`` parameter passed to :c:func:`FIRKodeSetUserData`.

   **Return value:**
      Should return 0 if successful, or -1 if unsuccessful.

   **Notes:**
      Allocation of memory for ``ewt`` is handled within FIRKODE.

   .. warning::

      The error weight vector must have all components positive. It is the
      user's responsibility to perform this test and return -1 if it is not
      satisfied.

.. _FIRKODE.Usage.user_supplied.root:

Rootfinding function
--------------------

If a rootfinding problem is to be solved during the integration of the ODE
system, the user must supply a C function of type :c:type:`FIRKRootFn`,
defined as follows:

.. c:type:: int (*FIRKRootFn)(sunrealtype t, N_Vector y, sunrealtype* gout, void* user_data)

   This function implements a vector-valued function :math:`g(t,y)` such that
   the roots of the ``nrtfn`` components :math:`g_i(t,y)` are sought.

   **Arguments:**
      - ``t`` -- the current value of the independent variable.
      - ``y`` -- the current value of the dependent variable vector, :math:`y(t)`.
      - ``gout`` -- the output array of length ``nrtfn`` with components :math:`g_i(t,y)`.
      - ``user_data`` -- a pointer to user data, the same as the ``user_data`` parameter passed to :c:func:`FIRKodeSetUserData`.

   **Return value:**
      A :c:type:`FIRKRootFn` should return 0 if successful or a non-zero
      value if an error occurred (in which case the integration is halted and
      :c:func:`FIRKodeEvolve` returns ``FIRK_RTFUNC_FAIL``).

   **Notes:**
      Allocation of memory for ``gout`` is automatically handled within
      FIRKODE.

.. _FIRKODE.Usage.user_supplied.jac:

Jacobian construction (matrix-based linear solvers)
---------------------------------------------------

If a matrix-based linear solver module is used (i.e., a non-``NULL``
``SUNMatrix`` object was supplied to :c:func:`FIRKodeSetLinearSolver`), the
user may optionally provide a function of type :c:type:`FIRKLsJacFn` for
evaluating the Jacobian of the ODE right-hand side function (or an
approximation of it). :c:type:`FIRKLsJacFn` is defined as follows:

.. c:type:: int (*FIRKLsJacFn)(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix Jac, void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)

   This function computes the Jacobian matrix :math:`J = \dfrac{\partial f}{\partial y}` (or an approximation to it).

   **Arguments:**
      - ``t`` -- the current value of the independent variable.
      - ``y`` -- the current value of the dependent variable vector, namely the solution :math:`y_n` at the beginning of the step.
      - ``fy`` -- the current value of the vector :math:`f(t,y)`.
      - ``Jac`` -- the output Jacobian matrix.
      - ``user_data`` -- a pointer to user data, the same as the ``user_data`` parameter passed to :c:func:`FIRKodeSetUserData`.
      - ``tmp1, tmp2, tmp3`` -- are pointers to memory allocated for variables of type ``N_Vector`` which can be used by a :c:type:`FIRKLsJacFn` function as temporary storage or work space.

   **Return value:**
      Should return 0 if successful, a positive value if a recoverable error
      occurred (in which case FIRKODE will attempt to correct, while FIRKLS
      sets ``last_flag`` to ``FIRKLS_JACFUNC_RECVR``), or a negative value if
      it failed unrecoverably (in which case the integration is halted,
      :c:func:`FIRKodeEvolve` returns ``FIRK_LSETUP_FAIL`` and FIRKLS sets
      ``last_flag`` to ``FIRKLS_JACFUNC_UNRECVR``).

   **Notes:**
      Information regarding the structure of the specific ``SUNMatrix``
      structure (e.g. number of rows, upper/lower bandwidth, sparsity type)
      may be obtained through using the implementation-specific ``SUNMatrix``
      interface functions (see :numref:`SUNMatrix` for details).

      With direct linear solvers (i.e., linear solvers with type
      ``SUNLINEARSOLVER_DIRECT``), the Jacobian matrix :math:`J(t,y)` is
      zeroed out prior to calling the user-supplied Jacobian function so only
      nonzero elements need to be loaded into ``Jac``.

      Unlike in CVODE and ARKODE, the Jacobian is always evaluated at the
      accepted solution :math:`(t_n, y_n)` at the beginning of a step, never
      at a predicted or stage value. Each call to the user's
      :c:type:`FIRKLsJacFn` function is preceded by a call to the
      :c:type:`FIRKRhsFn` user function with the same ``(t,y)`` arguments.
      Thus, the Jacobian function can use any auxiliary data that is computed
      and saved during the evaluation of the ODE right-hand side.

      If the user's :c:type:`FIRKLsJacFn` function uses difference quotient
      approximations, then it may need to access quantities not in the
      argument list. These include the current step size, the error weights,
      etc. To obtain these, the user will need to add a pointer to
      ``firkode_mem`` in ``user_data`` and then use the ``FIRKodeGet*``
      functions described in :numref:`FIRKODE.Usage.optional_output`. The
      unit roundoff can be accessed as ``SUN_UNIT_ROUNDOFF`` defined in
      ``sundials_types.h``.

      **Dense**: A user-supplied dense Jacobian function must load the
      :math:`N` by :math:`N` dense matrix ``Jac`` with an approximation to the
      Jacobian matrix :math:`J(t,y)` at the point :math:`(t, y)`. The accessor
      macros ``SM_ELEMENT_D`` and ``SM_COLUMN_D`` allow the user to read and
      write dense matrix elements without making explicit references to the
      underlying representation of the SUNMATRIX_DENSE type (see
      :numref:`SUNMatrix.Dense`).

      **Banded**: A user-supplied banded Jacobian function must load the
      :math:`N` by :math:`N` banded matrix ``Jac`` with the elements of the
      Jacobian :math:`J(t,y)` at the point :math:`(t,y)`. The accessor macros
      ``SM_ELEMENT_B``, ``SM_COLUMN_B``, and ``SM_COLUMN_ELEMENT_B`` allow the
      user to read and write band matrix elements without making specific
      references to the underlying representation of the SUNMATRIX_BAND type
      (see :numref:`SUNMatrix.Band`).

      **Sparse**: A user-supplied sparse Jacobian function must load the
      :math:`N` by :math:`N` compressed-sparse-column or compressed-sparse-row
      matrix ``Jac`` with an approximation to the Jacobian matrix
      :math:`J(t,y)` at the point :math:`(t,y)`. Storage for ``Jac`` already
      exists on entry to this function, although the user should ensure that
      sufficient space is allocated in ``Jac`` to hold the nonzero values to
      be set; if the existing space is insufficient the user may reallocate
      the data and index arrays as needed. The amount of allocated space in a
      SUNMATRIX_SPARSE object may be accessed using the macro ``SM_NNZ_S`` or
      the routine :c:func:`SUNSparseMatrix_NNZ` (see :numref:`SUNMatrix.Sparse`).

.. _FIRKODE.Usage.user_supplied.linsys:

Linear system construction (matrix-based linear solvers)
--------------------------------------------------------

With matrix-based linear solver modules, as an alternative to optionally
supplying a function for evaluating the Jacobian of the ODE right-hand side
function, the user may optionally supply a function of type
:c:type:`FIRKLsLinSysFn` for evaluating the linear system :math:`\mathcal{A}
= M - \gamma J` (or an approximation of it). :c:type:`FIRKLsLinSysFn` is
defined as follows:

.. c:type:: int (*FIRKLsLinSysFn)(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix A, SUNMatrix M, sunbooleantype jok, sunbooleantype* jcur, sunrealtype gamma, void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)

   This function computes the linear system matrix :math:`\mathcal{A} = M -
   \gamma J` (or an approximation to it).

   **Arguments:**
      - ``t`` -- the current value of the independent variable.
      - ``y`` -- the current value of the dependent variable vector, namely the solution :math:`y_n` at the beginning of the step.
      - ``fy`` -- the current value of the vector :math:`f(t,y)`.
      - ``A`` -- the output linear system matrix.
      - ``M`` -- the mass matrix supplied to :c:func:`FIRKodeSetMassLinearSolver`, or ``NULL`` if the mass matrix is the identity.
      - ``jok`` -- an input flag indicating whether the Jacobian-related data needs to be updated. The ``jok`` argument provides for the reuse of Jacobian data. When ``jok = SUNFALSE``, the Jacobian-related data should be recomputed from scratch. When ``jok = SUNTRUE`` the Jacobian data, if saved from the previous call to this function, can be reused (with the current value of ``gamma``). A call with ``jok = SUNTRUE`` can only occur after a call with ``jok = SUNFALSE``.
      - ``jcur`` -- a pointer to a flag which should be set to ``SUNTRUE`` if Jacobian data was recomputed, or set to ``SUNFALSE`` if Jacobian data was not recomputed, but saved data was still reused.
      - ``gamma`` -- the scalar :math:`\gamma` appearing in the matrix :math:`\mathcal{A} = M - \gamma J`.
      - ``user_data`` -- a pointer to user data, the same as the ``user_data`` parameter passed to :c:func:`FIRKodeSetUserData`.
      - ``tmp1, tmp2, tmp3`` -- are pointers to memory allocated for variables of type ``N_Vector`` which can be used by a :c:type:`FIRKLsLinSysFn` function as temporary storage or work space.

   **Return value:**
      Should return 0 if successful, a positive value if a recoverable error
      occurred (in which case FIRKODE will attempt to correct, while FIRKLS
      sets ``last_flag`` to ``FIRKLS_JACFUNC_RECVR``), or a negative value if
      it failed unrecoverably (in which case the integration is halted,
      :c:func:`FIRKodeEvolve` returns ``FIRK_LSETUP_FAIL`` and FIRKLS sets
      ``last_flag`` to ``FIRKLS_JACFUNC_UNRECVR``).

   **Notes:**
      Since the same matrix :math:`\mathcal{A}` is used as the block
      preconditioner of the Newton iteration and in the error estimate, a
      user-supplied linear system function must produce a matrix that is
      usable for both purposes. When the method has more than one stage, the
      Jacobian-vector products required by the Newton iteration are computed
      independently of this function (see :c:func:`FIRKodeSetJacTimes`).

.. _FIRKODE.Usage.user_supplied.jtimes:

Jacobian-vector product (matrix-free linear solvers)
----------------------------------------------------

If a matrix-free linear solver is to be used (i.e., a ``NULL``-valued
``SUNMatrix`` was supplied to :c:func:`FIRKodeSetLinearSolver`), or if the
method has more than one stage and the user wishes to supply the products
required by the coupled stage system, the user may provide a function of
type :c:type:`FIRKLsJacTimesVecFn` in the following form, to compute
matrix-vector products :math:`Jv`. If such a function is not supplied, the
default is a difference quotient approximation to these products (or, for
matrix-based solvers, the product with the saved Jacobian matrix).

.. c:type:: int (*FIRKLsJacTimesVecFn)(N_Vector v, N_Vector Jv, sunrealtype t, N_Vector y, N_Vector fy, void* user_data, N_Vector tmp)

   This function computes the product :math:`Jv = \dfrac{\partial f(t,y)}{\partial y} v` (or an approximation to it).

   **Arguments:**
      - ``v`` -- the vector by which the Jacobian must be multiplied.
      - ``Jv`` -- the output vector computed.
      - ``t`` -- the current value of the independent variable.
      - ``y`` -- the current value of the dependent variable vector, namely the solution :math:`y_n` at the beginning of the step.
      - ``fy`` -- the current value of the vector :math:`f(t,y)`.
      - ``user_data`` -- a pointer to user data, the same as the ``user_data`` parameter passed to :c:func:`FIRKodeSetUserData`.
      - ``tmp`` -- a pointer to memory allocated for a variable of type ``N_Vector`` which can be used for work space.

   **Return value:**
      The value returned by the Jacobian-vector product function should be 0
      if successful. Any other return value will result in an unrecoverable
      error of the Krylov solver, in which case the integration is halted.

   **Notes:**
      This function must return a value of :math:`Jv` that uses the
      *current* value of :math:`J`, i.e., as evaluated at the current
      :math:`(t,y)`, which is the beginning of the step.

      If the user's :c:type:`FIRKLsJacTimesVecFn` function uses difference
      quotient approximations, it may need to access quantities not in the
      argument list. These include the current step size, the error weights,
      etc. To obtain these, the user will need to add a pointer to
      ``firkode_mem`` to ``user_data`` and then use the ``FIRKodeGet*``
      functions described in :numref:`FIRKODE.Usage.optional_output`. The
      unit roundoff can be accessed as ``SUN_UNIT_ROUNDOFF`` defined in
      ``sundials_types.h``.

.. _FIRKODE.Usage.user_supplied.jtsetup:

Jacobian-vector product setup (matrix-free linear solvers)
----------------------------------------------------------

If the user's Jacobian-times-vector routine requires that any
Jacobian-related data be preprocessed or evaluated, then this needs to be
done in a user-supplied function of type :c:type:`FIRKLsJacTimesSetupFn`,
defined as follows:

.. c:type:: int (*FIRKLsJacTimesSetupFn)(sunrealtype t, N_Vector y, N_Vector fy, void* user_data)

   This function preprocesses and/or evaluates Jacobian data needed by the
   Jacobian-times-vector routine.

   **Arguments:**
      - ``t`` -- the current value of the independent variable.
      - ``y`` -- the current value of the dependent variable vector.
      - ``fy`` -- the current value of the vector :math:`f(t,y)`.
      - ``user_data`` -- a pointer to user data, the same as the ``user_data`` parameter passed to :c:func:`FIRKodeSetUserData`.

   **Return value:**
      The value returned by the Jacobian-vector setup function should be 0 if
      successful, positive for a recoverable error (in which case the step
      will be retried), or negative for an unrecoverable error (in which case
      the integration is halted).

   **Notes:**
      Each call to the Jacobian-vector setup function is preceded by a call
      to the :c:type:`FIRKRhsFn` user function with the same :math:`(t,y)`
      arguments. Thus, the setup function can use any auxiliary data that is
      computed and saved during the evaluation of the ODE right-hand side.
      The setup function is called once per Newton linear solve, before the
      Jacobian-vector products of that solve.

.. _FIRKODE.Usage.user_supplied.psolve:

Preconditioner solve (iterative linear solvers)
-----------------------------------------------

If a user-supplied preconditioner is to be used with a ``SUNLinearSolver``
module, then the user must provide a function to solve the linear system
:math:`Pz = r`, where :math:`P` may be either a left or right preconditioner
matrix. Here :math:`P` should approximate (at least crudely) the matrix
:math:`\mathcal{A} = M - \gamma J`, where :math:`J = \dfrac{\partial
f}{\partial y}`. If preconditioning is done on both sides, the product of the
two preconditioner matrices should approximate :math:`\mathcal{A}`. This
function must be of type :c:type:`FIRKLsPrecSolveFn`, defined as follows:

.. c:type:: int (*FIRKLsPrecSolveFn)(sunrealtype t, N_Vector y, N_Vector fy, N_Vector r, N_Vector z, sunrealtype gamma, sunrealtype delta, int lr, void* user_data)

   This function solves the preconditioned system :math:`Pz = r`.

   **Arguments:**
      - ``t`` -- the current value of the independent variable.
      - ``y`` -- the current value of the dependent variable vector.
      - ``fy`` -- the current value of the vector :math:`f(t,y)`.
      - ``r`` -- the right-hand side vector of the linear system.
      - ``z`` -- the computed output vector.
      - ``gamma`` -- the scalar :math:`\gamma` in the matrix given by :math:`\mathcal{A} = M - \gamma J`.
      - ``delta`` -- an input tolerance to be used if an iterative method is employed in the solution. In that case, the residual vector :math:`Res = r - Pz` of the system should be made less than ``delta`` in the weighted :math:`l_2` norm, i.e., :math:`\sqrt{\sum_i (Res_i \cdot ewt_i)^2 } <` ``delta``. To obtain the ``N_Vector`` ``ewt``, call :c:func:`FIRKodeGetErrWeights`.
      - ``lr`` -- an input flag indicating whether the preconditioner solve function is to use the left preconditioner (``lr = 1``) or the right preconditioner (``lr = 2``).
      - ``user_data`` -- a pointer to user data, the same as the ``user_data`` parameter passed to :c:func:`FIRKodeSetUserData`.

   **Return value:**
      The value returned by the preconditioner solve function is a flag
      indicating whether it was successful. This value should be 0 if
      successful, positive for a recoverable error (in which case the step
      will be retried), or negative for an unrecoverable error (in which case
      the integration is halted).

   **Notes:**
      When the user's linear solver is matrix-free and the method has more
      than one stage, this function is also called directly by FIRKODE, once
      per stage with ``lr = 2``, as the block preconditioner of the internal
      Krylov iteration on the coupled stage system.

.. _FIRKODE.Usage.user_supplied.psetup:

Preconditioner setup (iterative linear solvers)
-----------------------------------------------

If the user's preconditioner requires that any Jacobian-related data be
preprocessed or evaluated, then this needs to be done in a user-supplied
function of type :c:type:`FIRKLsPrecSetupFn`, defined as follows:

.. c:type:: int (*FIRKLsPrecSetupFn)(sunrealtype t, N_Vector y, N_Vector fy, sunbooleantype jok, sunbooleantype* jcurPtr, sunrealtype gamma, void* user_data)

   This function preprocesses and/or evaluates Jacobian-related data needed by the preconditioner.

   **Arguments:**
      - ``t`` -- the current value of the independent variable.
      - ``y`` -- the current value of the dependent variable vector, namely the solution :math:`y_n` at the beginning of the step.
      - ``fy`` -- the current value of the vector :math:`f(t,y)`.
      - ``jok`` -- an input flag indicating whether the Jacobian-related data needs to be updated. The ``jok`` argument provides for the reuse of Jacobian data in the preconditioner solve function. ``jok = SUNFALSE`` means that the Jacobian-related data must be recomputed from scratch. ``jok = SUNTRUE`` means that the Jacobian data, if saved from the previous call to this function, can be reused (with the current value of :math:`\gamma`). A call with ``jok = SUNTRUE`` can only occur after a call with ``jok = SUNFALSE``.
      - ``jcurPtr`` -- a pointer to a flag which should be set to ``SUNTRUE`` if Jacobian data was recomputed, or set to ``SUNFALSE`` if Jacobian data was not recomputed, but saved data was still reused.
      - ``gamma`` -- the scalar :math:`\gamma` appearing in the matrix :math:`\mathcal{A} = M - \gamma J`.
      - ``user_data`` -- a pointer to user data, the same as the ``user_data`` parameter passed to :c:func:`FIRKodeSetUserData`.

   **Return value:**
      The value returned by the preconditioner setup function is a flag
      indicating whether it was successful. This value should be 0 if
      successful, positive for a recoverable error (in which case the step
      will be retried), or negative for an unrecoverable error (in which case
      the integration is halted).

   **Notes:**
      The operations performed by this function might include forming a
      crude approximate Jacobian and performing an LU factorization of the
      resulting approximation to :math:`\mathcal{A} = M - \gamma J`.

      Each call to the preconditioner setup function is preceded by a call
      to the :c:type:`FIRKRhsFn` user function with the same :math:`(t,y)`
      arguments. Thus, the preconditioner setup function can use any
      auxiliary data that is computed and saved during the evaluation of the
      ODE right-hand side.

      This function is not called in advance of every call to the
      preconditioner solve function, but rather is called only as often as
      needed to achieve convergence in the Newton iteration (see
      :numref:`FIRKODE.Mathematics.linear`).

.. _FIRKODE.Usage.user_supplied.mass:

Mass matrix construction (matrix-based mass linear solvers)
-----------------------------------------------------------

If a matrix-based mass linear solver is used (i.e., a non-``NULL``
``SUNMatrix`` object was supplied to :c:func:`FIRKodeSetMassLinearSolver`),
the user must provide a function of type :c:type:`FIRKLsMassFn` to fill the
mass matrix.

.. c:type:: int (*FIRKLsMassFn)(sunrealtype t, SUNMatrix M, void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)

   This function computes the mass matrix :math:`M`.

   **Arguments:**
      - ``t`` -- the current value of the independent variable (the initial time, since the mass matrix is constant).
      - ``M`` -- the output mass matrix.
      - ``user_data`` -- a pointer to user data, the same as the ``user_data`` parameter passed to :c:func:`FIRKodeSetUserData`.
      - ``tmp1, tmp2, tmp3`` -- pointers to memory allocated to variables of type ``N_Vector`` which can be used by a :c:type:`FIRKLsMassFn` as temporary storage or work space.

   **Return value:**
      A :c:type:`FIRKLsMassFn` function should return 0 if successful, a
      positive value if a recoverable error occurred, or a negative value if
      it failed unrecoverably. Any nonzero return value causes
      :c:func:`FIRKodeEvolve` to return ``FIRK_MASSSETUP_FAIL`` (and FIRKLS
      sets ``last_flag`` to ``FIRKLS_MASSFUNC_RECVR`` or
      ``FIRKLS_MASSFUNC_UNRECVR``), since the mass matrix is constructed only
      once.

   **Notes:**
      Information regarding the structure of the specific ``SUNMatrix``
      structure (e.g. number of rows, upper/lower bandwidth, sparsity type)
      may be obtained through using the implementation-specific ``SUNMatrix``
      interface functions (see :numref:`SUNMatrix` for details).

      With direct linear solvers, the mass matrix is zeroed out prior to
      calling the user-supplied mass matrix function, so only nonzero
      elements need to be loaded into ``M``. The mass matrix must have the
      same ``SUNMatrix`` type as the Jacobian template supplied to
      :c:func:`FIRKodeSetLinearSolver`, and its sparsity pattern (for sparse
      matrices) must be compatible with that of the Jacobian, since the
      matrix :math:`M - \gamma J` is formed in the storage of the Jacobian
      template.

.. _FIRKODE.Usage.user_supplied.mtimes:

Mass matrix-vector product (matrix-free mass linear solvers)
------------------------------------------------------------

If a matrix-free mass linear solver is used (i.e., a ``NULL`` ``SUNMatrix``
was supplied to :c:func:`FIRKodeSetMassLinearSolver`), the user *must*
provide a function of type :c:type:`FIRKLsMassTimesVecFn` to compute mass
matrix-vector products :math:`Mv`. Such a function may also be supplied for
a matrix-based mass linear solver, in which case it replaces the default
product computed with :c:func:`SUNMatMatvec`.

.. c:type:: int (*FIRKLsMassTimesVecFn)(N_Vector v, N_Vector Mv, sunrealtype t, void* mtimes_data)

   This function computes the product :math:`Mv`.

   **Arguments:**
      - ``v`` -- the vector to multiply.
      - ``Mv`` -- the output vector computed.
      - ``t`` -- the current value of the independent variable.
      - ``mtimes_data`` -- a pointer to user data, the same as the ``mtimes_data`` parameter passed to :c:func:`FIRKodeSetMassTimes`.

   **Return value:**
      The value returned by the mass matrix-vector product function should be
      0 if successful. Any other return value will result in an unrecoverable
      error, in which case the integration is halted.

   **Notes:**
      Mass matrix-vector products are required in every Newton residual
      evaluation, in every iteration of the internal Krylov solver on the
      coupled stage system, and in the error estimate, so this function
      should be efficient.

.. _FIRKODE.Usage.user_supplied.mtsetup:

Mass matrix-vector product setup (matrix-free mass linear solvers)
------------------------------------------------------------------

If the user's mass matrix-vector product routine requires that any
mass-matrix-related data be preprocessed or evaluated, then this needs to be
done in a user-supplied function of type :c:type:`FIRKLsMassTimesSetupFn`,
defined as follows:

.. c:type:: int (*FIRKLsMassTimesSetupFn)(sunrealtype t, void* mtimes_data)

   This function preprocesses and/or evaluates mass-matrix-related data
   needed by the mass matrix-vector product routine.

   **Arguments:**
      - ``t`` -- the current value of the independent variable.
      - ``mtimes_data`` -- a pointer to user data, the same as the ``mtimes_data`` parameter passed to :c:func:`FIRKodeSetMassTimes`.

   **Return value:**
      The value returned by the mass matrix-vector setup function should be 0
      if successful. Any other return value causes :c:func:`FIRKodeEvolve` to
      return ``FIRK_MASSSETUP_FAIL``.

   **Notes:**
      Since the mass matrix is constant, this function is called at most
      once, at the beginning of the integration.

.. _FIRKODE.Usage.user_supplied.mpsolve:

Mass matrix preconditioner solve (iterative mass linear solvers)
----------------------------------------------------------------

If a user-supplied preconditioner is to be used with an iterative mass
linear solver, then the user must provide a function of type
:c:type:`FIRKLsMassPrecSolveFn` to solve the linear system :math:`Pz = r`,
where :math:`P` may be either a left or right preconditioner matrix that
approximates the mass matrix :math:`M`.

.. c:type:: int (*FIRKLsMassPrecSolveFn)(sunrealtype t, N_Vector r, N_Vector z, sunrealtype delta, int lr, void* user_data)

   This function solves the preconditioned system :math:`Pz = r`.

   **Arguments:**
      - ``t`` -- the current value of the independent variable.
      - ``r`` -- the right-hand side vector of the linear system.
      - ``z`` -- the computed output vector.
      - ``delta`` -- an input tolerance to be used if an iterative method is employed in the solution. In that case, the residual vector :math:`Res = r - Pz` of the system should be made less than ``delta`` in the weighted :math:`l_2` norm, i.e., :math:`\sqrt{\sum_i (Res_i \cdot ewt_i)^2 } <` ``delta``. To obtain the ``N_Vector`` ``ewt``, call :c:func:`FIRKodeGetErrWeights`.
      - ``lr`` -- an input flag indicating whether the preconditioner solve function is to use the left preconditioner (``lr = 1``) or the right preconditioner (``lr = 2``).
      - ``user_data`` -- a pointer to user data, the same as the ``user_data`` parameter passed to :c:func:`FIRKodeSetUserData`.

   **Return value:**
      The value returned by the preconditioner solve function is a flag
      indicating whether it was successful. This value should be 0 if
      successful, positive for a recoverable error, or negative for an
      unrecoverable error (in which case the integration is halted).

.. _FIRKODE.Usage.user_supplied.mpsetup:

Mass matrix preconditioner setup (iterative mass linear solvers)
----------------------------------------------------------------

If the user's mass matrix preconditioner requires that any data be
preprocessed or evaluated, then this needs to be done in a user-supplied
function of type :c:type:`FIRKLsMassPrecSetupFn`, defined as follows:

.. c:type:: int (*FIRKLsMassPrecSetupFn)(sunrealtype t, void* user_data)

   This function preprocesses and/or evaluates data needed by the mass
   matrix preconditioner.

   **Arguments:**
      - ``t`` -- the current value of the independent variable.
      - ``user_data`` -- a pointer to user data, the same as the ``user_data`` parameter passed to :c:func:`FIRKodeSetUserData`.

   **Return value:**
      The value returned by the preconditioner setup function is a flag
      indicating whether it was successful. This value should be 0 if
      successful, positive for a recoverable error, or negative for an
      unrecoverable error (in which case the integration is halted).

   **Notes:**
      Since the mass matrix is constant, this function is called at most
      once.
