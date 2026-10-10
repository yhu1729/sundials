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
The tables of the Radau IIA methods with :math:`1 \le s \le`
:c:macro:`FIRK_MAX_STAGES` stages are built in: they are decimal literals with
40 significant digits, generated in 50-digit arithmetic by the script
``scripts/firkode_radau_tables.py``, which checks the order conditions and the
identities below before writing them, so that the tables are identical on every
platform and nothing is computed at run time. The utility routines below are
provided for users who wish to inspect the coefficients, e.g., through
:c:func:`FIRKodeGetCurrentTable`, or to verify the tables.

.. c:type:: FIRKodeTableMem* FIRKodeTable

   Pointer to the structure holding an :math:`s`-stage fully implicit
   Runge--Kutta table.

   .. versionadded:: x.y.z


.. c:struct:: FIRKodeTableMem

   Structure holding the coefficients of an :math:`s`-stage fully implicit
   Runge--Kutta method and the derived quantities used by FIRKODE. Its
   members are:

   - ``int s`` -- the number of stages.
   - ``int q`` -- the order of accuracy of the method (:math:`2s-1` for Radau IIA).
   - ``int p`` -- the order of the error estimate (:math:`s` for Radau IIA).
   - ``sunrealtype c[FIRK_MAX_STAGES]`` -- the stage nodes, in the first ``s`` entries.
   - ``sunrealtype A[FIRK_MAX_STAGES * FIRK_MAX_STAGES]`` -- the Butcher table coefficients, stored row-major: :math:`a_{ij}` is ``A[i * s + j]``.
   - ``sunrealtype b[FIRK_MAX_STAGES]`` -- the solution weights.
   - ``sunrealtype Ainv[FIRK_MAX_STAGES * FIRK_MAX_STAGES]`` -- the inverse :math:`A^{-1}`, stored like ``A``.
   - ``sunrealtype d[FIRK_MAX_STAGES]`` -- the update weights :math:`d = b^T A^{-1}` of :eq:`FIRKODE_update`.
   - ``sunrealtype e[FIRK_MAX_STAGES]`` -- the error estimate weights :math:`e_i` of :eq:`FIRKODE_errest`.
   - ``sunrealtype gamma0`` -- the real shift :math:`\gamma_0` of :eq:`FIRKODE_precond`.
   - ``sunrealtype P[FIRK_MAX_STAGES * FIRK_MAX_STAGES]`` -- the dense output coefficients; ``P[k * s + i]`` is the coefficient of :math:`\theta^{k+1}` in the Lagrange basis polynomial :math:`L_i(\theta)` of :eq:`FIRKODE_dense`, for :math:`k = 0, \ldots, s-1`.
   - ``sunbooleantype stiffly_accurate`` -- ``SUNTRUE`` if :math:`b` equals the last row of :math:`A`.

   Every table returned by :c:func:`FIRKodeTable_RadauIIA`,
   :c:func:`FIRKodeTable_Alloc` or :c:func:`FIRKodeTable_Copy` is owned by
   the caller, who frees it with :c:func:`FIRKodeTable_Free`.

.. c:macro:: FIRK_MAX_STAGES

   The maximum number of stages of a table (9).

.. c:function:: FIRKodeTable FIRKodeTable_RadauIIA(int s)

   Returns a copy of the built-in table of the :math:`s`-stage Radau IIA
   method, including all derived quantities.

   :param s: the number of stages, :math:`1 \le s \le` :c:macro:`FIRK_MAX_STAGES`.

   :returns: A pointer to a newly allocated table, or ``NULL`` if ``s`` is out of
      range or the allocation failed.

   .. versionadded:: x.y.z

   **Notes:**
      Every call allocates a new table that the caller owns, may modify, and
      frees with :c:func:`FIRKodeTable_Free`, as for the ARKODE Butcher
      tables. For the stiffly accurate Radau IIA methods :math:`b` equals the
      last row of :math:`A` and :math:`d` is the last unit vector exactly.

.. c:function:: FIRKodeTable FIRKodeTable_Alloc(int s)

   Allocates an empty table with ``s`` stages, with all coefficient arrays
   set to zero.

   :param s: the number of stages, :math:`1 \le s \le` :c:macro:`FIRK_MAX_STAGES`.

   :returns: A pointer to the newly allocated table, or ``NULL`` on failure.

   .. versionadded:: x.y.z


.. c:function:: FIRKodeTable FIRKodeTable_Copy(FIRKodeTable T)

   Creates a copy of the table ``T``, which the caller owns and may modify.

   :param T: the table to copy, for example a built-in table.

   :returns: A pointer to the newly allocated copy, or ``NULL`` on failure.

   .. versionadded:: x.y.z


.. c:function:: void FIRKodeTable_Free(FIRKodeTable T)

   Frees a table returned by :c:func:`FIRKodeTable_RadauIIA`,
   :c:func:`FIRKodeTable_Alloc` or :c:func:`FIRKodeTable_Copy`.

   :param T: the table to free; a ``NULL`` argument is ignored. The table
      returned by :c:func:`FIRKodeGetCurrentTable` belongs to the
      integrator and must not be passed here.

   .. versionadded:: x.y.z


.. c:function:: void FIRKodeTable_Write(FIRKodeTable T, FILE* outfile)

   Writes every member of the table ``T`` to the file ``outfile``.

   :param T: the table to write.
   :param outfile: the output file pointer.

   .. versionadded:: x.y.z


.. c:function:: int FIRKodeTable_CheckOrder(FIRKodeTable T, int* q_ptr, int* p_ptr, int* stage_order_ptr, FILE* outfile)

   Determines the orders of the method defined by the table ``T`` by checking
   the order conditions numerically, and compares them with the orders stored
   in the table.

   :param T: the table to check.
   :param q_ptr: on output, the order of the method,
      :math:`\min(q_B, 2\eta + 2, \eta + \zeta + 1)` from the simplifying
      assumptions :math:`B(q_B)`, :math:`C(\eta)` and :math:`D(\zeta)`.
   :param p_ptr: on output, the order of the error estimate: the largest
      :math:`k` such that the embedded quadrature, with the weights
      :math:`b^T + e^T A` on the nodes :math:`c` and :math:`\gamma_0` on the
      node 0, integrates polynomials of degree :math:`k - 1` exactly.
   :param stage_order_ptr: on output, the stage order :math:`\min(\eta, q_B)`.
   :param outfile: if non-``NULL``, a file pointer to which the results of
      the individual condition checks are written.

   :returns: 0 if the determined orders ``q_ptr`` and ``p_ptr`` match the orders stored in
      the table, 1 if either is higher than the stored order (a warning), and
      -1 if either is lower than the stored order or an argument is invalid.

   .. versionadded:: x.y.z

   **Notes:**
      A condition is taken to hold if its residual is at most
      :math:`1000 s^2 u`, with :math:`u` the unit roundoff. In double and
      extended precision this separates satisfied from violated conditions for
      every built-in table; in single precision only for :math:`s \le 3`.
