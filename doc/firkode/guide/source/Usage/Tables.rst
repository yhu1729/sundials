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

.. _FIRKODE.Usage.Tables:

Radau IIA coefficient tables
============================

The coefficient tables of the Radau IIA methods, together with the derived
quantities used by the integrator (see :numref:`FIRKODE.Mathematics`), are
stored in objects of type :c:type:`FIRKodeTable`, declared in the header file
``firkode/firkode_tables.h`` (which is included by ``firkode/firkode.h``).
FIRKODE constructs the table of the selected method internally; the utility
routines below are provided for users who wish to inspect the coefficients,
e.g., through :c:func:`FIRKodeGetCurrentTable`, or to verify the tables.

.. c:type:: FIRKodeTableMem* FIRKodeTable

   Pointer to the structure holding an :math:`s`-stage fully implicit
   Runge--Kutta table.

.. c:struct:: FIRKodeTableMem

   Structure holding the coefficients of an :math:`s`-stage fully implicit
   Runge--Kutta method and the derived quantities used by FIRKODE. Its
   members are:

   - ``int s`` -- the number of stages.
   - ``int q`` -- the order of accuracy of the method (:math:`2s-1` for Radau IIA).
   - ``int p`` -- the order of the error estimate (:math:`s` for Radau IIA).
   - ``sunrealtype* c`` -- the stage nodes, an array of length ``s``.
   - ``sunrealtype** A`` -- the Butcher table coefficients, an ``s`` by ``s`` array.
   - ``sunrealtype* b`` -- the solution weights, an array of length ``s``.
   - ``sunrealtype** Ainv`` -- the inverse :math:`A^{-1}`.
   - ``sunrealtype* d`` -- the update weights :math:`d = b^T A^{-1}` of :eq:`FIRKODE_update`.
   - ``sunrealtype* e`` -- the error estimate weights :math:`e_i` of :eq:`FIRKODE_errest`.
   - ``sunrealtype gamma0`` -- the real shift :math:`\gamma_0` of :eq:`FIRKODE_precond`.
   - ``sunrealtype** P`` -- the dense output coefficients; ``P[k][i]`` is the coefficient of :math:`\theta^{k+1}` in the Lagrange basis polynomial :math:`L_i(\theta)` of :eq:`FIRKODE_dense`, for :math:`k = 0, \ldots, s-1`.
   - ``sunbooleantype stiffly_accurate`` -- ``SUNTRUE`` if :math:`b` equals the last row of :math:`A`.

.. c:macro:: FIRK_MAX_STAGES

   The maximum number of stages for which tables can be constructed (9).

.. c:function:: FIRKodeTable FIRKodeTable_RadauIIA(int s)

   Constructs the table of the :math:`s`-stage Radau IIA method, including
   all derived quantities.

   **Arguments:**
      - ``s`` -- the number of stages, :math:`1 \le s \le` :c:macro:`FIRK_MAX_STAGES`.

   **Return value:**
      - A pointer to the newly allocated table, or ``NULL`` if ``s`` is out
        of range, if a memory allocation failed, or if the construction of
        the table failed.

   **Notes:**
      The tables for :math:`s = 1, 2, 3` are given in closed form. For
      :math:`s \ge 4` the nodes are computed as the roots of the Radau
      polynomial and the coefficients from the collocation conditions, in
      extended precision where available, and rounded to ``sunrealtype``. The
      table is owned by the caller and must be freed with
      :c:func:`FIRKodeTable_Free`.

.. c:function:: FIRKodeTable FIRKodeTable_Alloc(int s)

   Allocates an empty table with ``s`` stages, with all coefficient arrays
   set to zero.

   **Arguments:**
      - ``s`` -- the number of stages, :math:`1 \le s \le` :c:macro:`FIRK_MAX_STAGES`.

   **Return value:**
      - A pointer to the newly allocated table, or ``NULL`` on failure.

.. c:function:: FIRKodeTable FIRKodeTable_Copy(FIRKodeTable T)

   Creates a copy of the table ``T``.

   **Arguments:**
      - ``T`` -- the table to copy.

   **Return value:**
      - A pointer to the newly allocated copy, or ``NULL`` on failure.

.. c:function:: void FIRKodeTable_Free(FIRKodeTable T)

   Frees the table ``T`` and all of its arrays.

   **Arguments:**
      - ``T`` -- the table to free (a ``NULL`` argument is ignored).

.. c:function:: void FIRKodeTable_Write(FIRKodeTable T, FILE* outfile)

   Writes the coefficients of the table ``T`` and its derived quantities to
   the file ``outfile``.

   **Arguments:**
      - ``T`` -- the table to write.
      - ``outfile`` -- the output file pointer.

.. c:function:: int FIRKodeTable_CheckOrder(FIRKodeTable T, int* q, int* p, FILE* outfile)

   Determines the order of accuracy of the method defined by the table ``T``
   by checking the quadrature, collocation, and simplifying order conditions
   numerically, and compares the result with the orders stored in the table.

   **Arguments:**
      - ``T`` -- the table to check.
      - ``q`` -- on output, the order of accuracy determined from the order conditions.
      - ``p`` -- on output, the order of the error estimate determined from the order conditions.
      - ``outfile`` -- if non-``NULL``, a file pointer to which the results of the individual condition checks are written.

   **Return value:**
      - 0 if the determined orders match those stored in the table, 1 if the
        determined orders are higher than the stored ones (a warning), and -1
        if the determined orders are lower than the stored ones or the table
        is invalid.

   **Notes:**
      The order conditions are checked with a tolerance that scales with the
      number of stages and the unit roundoff; for large numbers of stages the
      check of high-order conditions is limited by roundoff in the table
      coefficients.
