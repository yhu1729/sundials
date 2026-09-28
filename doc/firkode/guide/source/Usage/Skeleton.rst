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

.. _FIRKODE.Usage.Skeleton:

A skeleton of the user's main program
=====================================

The following is a skeleton of the user's main program (or calling program)
for the integration of an ODE IVP. Most of the steps are independent of the
``N_Vector``, ``SUNMatrix``, and ``SUNLinearSolver`` implementations used.
For the steps that are not, refer to :numref:`NVectors`, :numref:`SUNMatrix`,
and :numref:`SUNLinSol` for the specific name of the function to be called or
macro to be referenced.

  #. **Initialize parallel or multi-threaded environment, if appropriate**
     For example, call ``MPI_Init`` to initialize MPI if used, or set the number
     of threads to use within the threaded vector functions if used.

  #. **Create the SUNDIALS context object**
     Call :c:func:`SUNContext_Create` to allocate the ``SUNContext`` object.

  #. **Set problem dimensions etc.**
     This generally includes the problem size ``N``, and may include the local
     vector length ``Nlocal``.

     Note: The variables ``N`` and ``Nlocal`` should be of type ``sunindextype``.

  #. **Set vector of initial values**
     To set the vector of initial values, use the appropriate functions
     defined by the particular ``N_Vector`` implementation.

     For native SUNDIALS vector implementations, use a call of the form
     ``y0 = N_VMake_***(..., ydata)`` if the array containing the initial
     values of :math:`y` already exists. Otherwise, create a new vector by
     making a call of the form ``N_VNew_***(...)``, and then set its elements
     by accessing the underlying data with a call of the form
     ``ydata = N_VGetArrayPointer(y0)``. See :numref:`NVectors` for details.

  #. **Create FIRKODE object**
     Call :c:func:`FIRKodeCreate` to create the FIRKODE memory block.
     :c:func:`FIRKodeCreate` returns a pointer to the FIRKODE memory structure.

  #. **Initialize FIRKODE solver**
     Call :c:func:`FIRKodeInit` to provide required problem specifications,
     allocate internal memory for FIRKODE, and initialize FIRKODE.
     :c:func:`FIRKodeInit` returns a flag, the value of which indicates either
     success or an illegal argument value.

  #. **Specify integration tolerances**
     Call :c:func:`FIRKodeSStolerances` or :c:func:`FIRKodeSVtolerances` to
     specify either a scalar relative tolerance and scalar absolute tolerance,
     or a scalar relative tolerance and a vector of absolute tolerances,
     respectively. Alternatively, call :c:func:`FIRKodeWFtolerances` to specify
     a function which sets directly the weights used in evaluating WRMS vector
     norms.

  #. **Select the method** (*optional*)
     Call :c:func:`FIRKodeSetNumStages` or :c:func:`FIRKodeSetOrder` to
     select the number of stages of the Radau IIA method. The default is the
     three-stage, fifth-order method.

  #. **Create matrix object**
     If the linear solver will be a matrix-based linear solver, then a
     template Jacobian matrix must be created by calling the appropriate
     constructor function defined by the particular ``SUNMatrix``
     implementation.

     For the native SUNDIALS ``SUNMatrix`` implementations, the matrix object
     may be created using a call of the form ``SUN***Matrix(...)`` where
     ``***`` is the name of the matrix (see :numref:`SUNMatrix` for details).

  #. **Create linear solver object**
     The desired linear solver object must be created by calling the
     appropriate constructor function defined by the particular
     ``SUNLinearSolver`` implementation.

     For any of the SUNDIALS-supplied ``SUNLinearSolver`` implementations, the
     linear solver object may be created using a call of the form
     ``SUNLinearSolver LS = SUNLinSol_*(...);`` where ``*`` can be replaced
     with "Dense", "SPGMR", or other options, as discussed in
     :numref:`FIRKODE.Usage.linear_solver` and :numref:`SUNLinSol`.

  #. **Set linear solver optional inputs**
     Call functions from the selected linear solver module to change optional
     inputs specific to that linear solver. See the documentation for each
     ``SUNLinearSolver`` module in :numref:`SUNLinSol` for details.

  #. **Attach linear solver module**
     Initialize the FIRKLS linear solver interface by attaching the linear
     solver object (and matrix object, if applicable) with a call
     ``ier = FIRKodeSetLinearSolver(firkode_mem, LS, A)`` (for details see
     :numref:`FIRKODE.Usage.linear_solver`). This step is required.

  #. **Attach mass matrix solver module** (*optional*)
     For problems with a non-identity mass matrix, create a second matrix
     object (if applicable) and linear solver object for the mass matrix and
     attach them with a call ``ier = FIRKodeSetMassLinearSolver(firkode_mem,
     MLS, M, SUNFALSE)`` (for details see :numref:`FIRKODE.Usage.mass_solver`).
     Then supply the mass matrix routine with :c:func:`FIRKodeSetMassFn` or
     the mass matrix-vector product routine with :c:func:`FIRKodeSetMassTimes`.

  #. **Set optional inputs**
     Call ``FIRKodeSet***`` functions to change any optional inputs that
     control the behavior of FIRKODE from their default values. See
     :numref:`FIRKODE.Usage.optional_input` for details.

  #. **Specify rootfinding problem** (*optional*)
     Call :c:func:`FIRKodeRootInit` to initialize a rootfinding problem to be
     solved during the integration of the ODE system. See
     :numref:`FIRKODE.Usage.rootinit`, and see
     :numref:`FIRKODE.Usage.optional_input.root` for relevant optional input
     calls.

  #. **Advance solution in time**
     For each point at which output is desired, call
     ``ier = FIRKodeEvolve(firkode_mem, tout, yout, &tret, itask)``. Here
     ``itask`` specifies the return mode. The vector ``yout`` (which can be
     the same as the vector ``y0`` above) will contain :math:`y(t)`. See
     :c:func:`FIRKodeEvolve` for details.

  #. **Get optional outputs**
     Call ``FIRKodeGet***`` functions to obtain optional output. See
     :numref:`FIRKODE.Usage.optional_output` for details.

  #. **Deallocate memory for solution vector**
     Upon completion of the integration, deallocate memory for the vector
     ``y`` (or ``yout``) by calling the appropriate destructor function
     defined by the ``N_Vector`` implementation.

  #. **Free solver memory**
     Call :c:func:`FIRKodeFree` to free the memory allocated by FIRKODE.

  #. **Free linear solver and matrix memory**
     Call :c:func:`SUNLinSolFree` and :c:func:`SUNMatDestroy` to free any
     memory allocated for the linear solver and matrix objects created above.

  #. **Free the SUNContext object**
     Call :c:func:`SUNContext_Free` to free the memory allocated for the
     ``SUNContext`` object.

  #. **Finalize MPI, if used**
     Call ``MPI_Finalize`` to terminate MPI.
