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

.. _FIRKODE.Mathematics:

***************************
Mathematical Considerations
***************************

FIRKODE solves ODE initial value problems (IVPs) in real :math:`N`-space,
which we write in the abstract form

.. math::
   M \dot{y} = f(t,y) \, , \qquad y(t_0) = y_0 \, ,
   :label: FIRKODE_ivp

where :math:`y \in \mathbb{R}^N` and :math:`M` is a constant, nonsingular
mass matrix (the identity unless the user supplies a mass matrix, see
:numref:`FIRKODE.Mathematics.mass`). Here, we use :math:`\dot{y}` to denote
:math:`\mathrm{d}y/\mathrm{d}t`. The independent variable :math:`t` may be
integrated forward or backward.

.. _FIRKODE.Mathematics.radau:

Radau IIA methods
=================

FIRKODE integrates :eq:`FIRKODE_ivp` with the fully implicit Radau IIA
Runge--Kutta methods :cite:p:`HaWa:91,Butcher:08`. The :math:`s`-stage method
is the collocation method with nodes :math:`0 < c_1 < \cdots < c_s = 1`
given by the roots of

.. math::
   P_s(2c - 1) - P_{s-1}(2c - 1) = 0 \, ,

where :math:`P_k` denotes the Legendre polynomial of degree :math:`k`, i.e.,
the nodes of the right Radau quadrature rule. The coefficients :math:`A =
(a_{ij})` are determined by the collocation conditions

.. math::
   \sum_{j=1}^{s} a_{ij}\, c_j^{k-1} = \frac{c_i^k}{k} \, , \qquad
   i, k = 1, \ldots, s \, ,

and the weights :math:`b_j = a_{sj}` coincide with the last row of
:math:`A`, so that the methods are *stiffly accurate*. The :math:`s`-stage
Radau IIA method has classical order :math:`2s - 1` and stage order
:math:`s`, and is A-stable and L-stable; in particular, the stability function
satisfies :math:`R(\infty) = 0`, so that stiff components are damped in a
single step. FIRKODE provides the tables for :math:`s = 1, 2, 3` (orders 1,
3, and 5; the one-stage method is the implicit Euler method and the
three-stage method is the method of RADAU5 :cite:p:`HaWa:99`) in closed form,
and constructs the tables for :math:`4 \le s \le 9` at run time from the
above characterization (see :numref:`FIRKODE.Usage.Tables`). The default
method has :math:`s = 3` stages.

Writing :math:`h = t_{n+1} - t_n` for the step size, a step of the method from
:math:`(t_n, y_n)` computes the stage increments :math:`Z_i = Y_i - y_n`,
:math:`i = 1, \ldots, s`, from the coupled nonlinear system

.. math::
   M Z_i = h \sum_{j=1}^{s} a_{ij}\, f(t_n + c_j h,\, y_n + Z_j) \, ,
   \qquad i = 1, \ldots, s \, ,
   :label: FIRKODE_stages

and then sets

.. math::
   y_{n+1} = y_n + \sum_{j=1}^{s} d_j Z_j \, , \qquad d = b^T A^{-1} \, .
   :label: FIRKODE_update

For the Radau IIA methods :math:`d = (0, \ldots, 0, 1)`, so that
:math:`y_{n+1} = y_n + Z_s = Y_s`. Formulating the step in terms of the
increments :math:`Z_i` rather than the stage values :math:`Y_i`, and updating
through :eq:`FIRKODE_update` rather than through the weights :math:`b`, avoids
the amplification of roundoff errors in the stiff regime :cite:p:`HaWa:91`.

.. _FIRKODE.Mathematics.Newton:

Nonlinear solution of the stage equations
=========================================

The system :eq:`FIRKODE_stages` of :math:`sN` equations is solved with a
*simplified Newton iteration* in which the Jacobian :math:`J = \partial f /
\partial y` is evaluated once at :math:`(t_n, y_n)` (and possibly reused
over several steps, see :numref:`FIRKODE.Mathematics.linear`). Writing
:math:`Z = (Z_1, \ldots, Z_s)` for the stacked vector of stage increments and
:math:`G(Z)` for the residual of :eq:`FIRKODE_stages`, the iteration reads

.. math::
   \left( I_s \otimes M - h A \otimes J \right) \delta^{(k)} = -G(Z^{(k)})
   \, , \qquad Z^{(k+1)} = Z^{(k)} + \delta^{(k)} \, .
   :label: FIRKODE_Newton

The stacked vectors :math:`Z`, :math:`\delta`, and :math:`G` are stored as
:ref:`NVECTOR_MANYVECTOR <NVectors.ManyVector>` objects whose :math:`s`
subvectors are clones of the user's vector, and the iteration is carried out
with the ``SUNNONLINSOL_NEWTON`` module.

The initial guess :math:`Z^{(0)}` is, by default, obtained by extrapolating
the collocation polynomial of the previous step (see
:numref:`FIRKODE.Mathematics.dense`) to the stage times of the new step; the
trivial predictor :math:`Z^{(0)} = 0` may be selected instead with
:c:func:`FIRKodeSetPredictorMethod`. On the first step, or after the
integrator has been reinitialized, the trivial predictor is used.

The convergence test follows RADAU5 :cite:p:`HaWa:91`. Let :math:`\|\cdot\|`
denote the weighted root-mean-square (WRMS) norm defined by the user's
tolerances (see :numref:`FIRKODE.Mathematics.error`), applied to all
:math:`sN` components of the stacked vectors. After the :math:`k`-th
iteration, with :math:`k \ge 2`, the contraction rate is estimated as
:math:`\theta_k = \|\delta^{(k)}\| / \|\delta^{(k-1)}\|`, and the iteration
is considered converged when

.. math::
   \frac{\theta_k}{1 - \theta_k}\, \|\delta^{(k)}\| \le \epsilon_{nls} \, ,
   :label: FIRKODE_nlstest

where the nonlinear convergence coefficient :math:`\epsilon_{nls}`
(default 0.1) may be set with :c:func:`FIRKodeSetNonlinConvCoef`. For the
first iteration, where no rate estimate is available, the factor
:math:`\theta_k / (1 - \theta_k)` is replaced by :math:`\max(\eta, U)^{0.8}`,
where :math:`\eta` is the last value of this factor from the previous
nonlinear solve and :math:`U` is the unit roundoff; a linear problem thus
typically converges in one or two iterations. The iteration is abandoned as
divergent if :math:`\theta_k \ge \theta_{\max}` (default 0.99, see
:c:func:`FIRKodeSetNonlinDivergenceRate`), and it is abandoned early when
convergence within the maximum number of iterations :math:`m` (default 7, see
:c:func:`FIRKodeSetMaxNonlinIters`) is predicted to be impossible, i.e., when

.. math::
   \frac{\theta_k}{1 - \theta_k}\, \theta_k^{\,m - 1 - k}\, \|\delta^{(k)}\|
   \ge \epsilon_{nls} \, .

In the latter case the step is retried with a step size reduced by the factor
:math:`0.8\, q^{-1/(4 + m - 1 - k)}`, where :math:`q` is the left-hand side
of the previous inequality divided by :math:`\epsilon_{nls}` (clipped to
:math:`[10^{-4}, 20]`); in all other failure cases the step size is reduced
by the factor set with :c:func:`FIRKodeSetMaxCFailGrowth` (default 0.25). A
nonlinear solver failure also forces the Jacobian and the linear solver
setup to be refreshed on the retried step.

.. _FIRKODE.Mathematics.linear:

Solution of the linear systems
==============================

The matrix of the Newton system :eq:`FIRKODE_Newton` couples all stages.
Rather than factoring the :math:`sN \times sN` matrix, FIRKODE solves each
Newton system with the right-preconditioned flexible GMRES iteration of the
:ref:`SUNLINSOL_SPFGMR <SUNLinSol.SPFGMR>` module, using the block-diagonal
preconditioner

.. math::
   P = I_s \otimes \left( M - \gamma J \right) \, , \qquad
   \gamma = h \gamma_0 \, ,
   :label: FIRKODE_precond

where :math:`\gamma_0` is the reciprocal of the real eigenvalue of
:math:`A^{-1}` for an odd number of stages, and the reciprocal of the real part
of the eigenvalue of :math:`A^{-1}` of smallest modulus for an even number of
stages. The preconditioner therefore requires only solves with a matrix of the
size of the ODE system, :math:`M - \gamma J`, exactly as in a diagonally
implicit method, and these solves are performed by the user's
:c:type:`SUNLinearSolver` object through the FIRKLS interface described in
:numref:`FIRKODE.Usage.linear_solver`. Every iteration of the flexible GMRES
method additionally requires one matrix-vector product with the Jacobian for
each stage, :math:`J v_j`, which is computed with the saved Jacobian matrix
(matrix-based linear solvers), with a user-supplied Jacobian-times-vector
routine, or with a difference quotient approximation (see
:c:func:`FIRKodeSetJacTimes`). For the one-stage method the preconditioner is
exact and the user's linear solver is applied directly, without the Krylov
iteration.

The Krylov iteration on the stacked system is stopped when the WRMS norm of
the preconditioned residual is below :math:`\epsilon_{stk}\, \epsilon_{nls}`,
where :math:`\epsilon_{stk}` (default 0.05) is set with
:c:func:`FIRKodeSetStageEpsLin`; its maximum Krylov subspace dimension and
number of restarts may be set with :c:func:`FIRKodeSetStageSolverMaxl` and
:c:func:`FIRKodeSetStageSolverMaxRestarts`. Because the block solves
performed by an iterative user linear solver may be inexact and may differ
from one iteration to the next, the *flexible* variant of GMRES is essential
here. When the user's linear solver is matrix-free, its preconditioner solve
routine is applied block-wise as the preconditioner :eq:`FIRKODE_precond`
directly, so that no Krylov iteration is nested inside another.

As in CVODE, the Jacobian :math:`J` and the linear solver setup (e.g., the
factorization of :math:`M - \gamma J`) are updated only periodically. The
linear solver setup is performed on the first step, after a nonlinear
solver failure or a failed error test, every :math:`msbp` steps (default 20,
see :c:func:`FIRKodeSetLSetupFrequency`), and whenever :math:`\gamma` has
changed by more than a fraction :math:`\Delta\gamma_{\max}` (default 0.2, see
:c:func:`FIRKodeSetDeltaGammaMaxLSetup`) since the last setup. Within a
setup, the Jacobian itself is re-evaluated on the first step, after a
nonlinear solver failure that was not attributed to a change in
:math:`\gamma`, every :math:`msbj` steps (default 51, see
:c:func:`FIRKodeSetJacEvalFrequency`), and optionally when the observed
Newton contraction rate exceeds a threshold (see
:c:func:`FIRKodeSetJacBadConvRate`); otherwise the saved Jacobian is reused
and only :math:`M - \gamma J` is rebuilt with the new :math:`\gamma`.

.. _FIRKODE.Mathematics.error:

Local error estimation and error test
=====================================

FIRKODE controls the local truncation error in the WRMS norm

.. math::
   \| v \|_{\text{WRMS}} = \left( \frac{1}{N} \sum_{i=1}^{N}
   \left( W_i\, v_i \right)^2 \right)^{1/2} \, , \qquad
   W_i = \frac{1}{\mathrm{rtol}\, |y_{n,i}| + \mathrm{atol}_i} \, ,
   :label: FIRKODE_errwt

where the weights are computed from the relative tolerance
:math:`\mathrm{rtol}` and the (scalar or vector) absolute tolerance
:math:`\mathrm{atol}` at the beginning of each step (or by a user-supplied
weight function, see :c:func:`FIRKodeWFtolerances`).

The local error estimate is the filtered embedded estimate of RADAU5
:cite:p:`HaWa:91,HaWa:99`, extended to an arbitrary number of stages.
Let :math:`u(t)` denote the collocation polynomial of the step, and let
:math:`\hat{b}` be the weights of the embedded quadrature rule of order
:math:`s` that includes the left endpoint of the step. The difference
:math:`(\hat{b} - b)^T A^{-1}` between the embedded and the collocation
solution can be written in terms of the stage increments, and the estimate is

.. math::
   \mathrm{err} = \left( M - \gamma J \right)^{-1}
   \left[ \gamma f(t_n, y_n) + \sum_{i=1}^{s} e_i\, M Z_i \right] \, ,
   \qquad \gamma = h \gamma_0 \, ,
   :label: FIRKODE_errest

with method constants :math:`e_i` that are stored with the coefficient table
(for :math:`s = 3` these are the constants of RADAU5). The estimate is of
order :math:`s + 1` and the linear solve applies the same matrix as the
preconditioner :eq:`FIRKODE_precond`, so that the estimate is bounded for
:math:`h \to \infty` for stiff problems. The estimate is scaled to the WRMS
norm, :math:`\mathrm{dsm} = \| \mathrm{err} \|_{\text{WRMS}}`, and the step is
accepted if :math:`\mathrm{dsm} \le 1`. As in RADAU5, when the test fails
on the first step or on a step immediately following a rejected step, the
estimate is filtered a second time by replacing :math:`f(t_n, y_n)` with
:math:`f(t_n, y_n + \mathrm{err})` in :eq:`FIRKODE_errest` before the test is
applied; this second filtering may be disabled with
:c:func:`FIRKodeSetErrorRefilter`. Rejected steps are retried with a reduced
step size, and the integration fails with the error code
``FIRK_ERR_FAILURE`` after a maximum number of consecutive error test
failures (default 7, see :c:func:`FIRKodeSetMaxErrTestFails`).

.. _FIRKODE.Mathematics.adaptivity:

Time step adaptivity
====================

After each step attempt the step size for the next step is proposed by a
:c:type:`SUNAdaptController` object (see :numref:`SUNAdaptController`),
supplied with the current step size, the value :math:`\mathrm{dsm}`, and the
order :math:`p = s` of the error estimate. The default controller is the
implicit Gustafsson controller (:numref:`SUNAdaptController.ImExGus`)
used by RADAU5; any other controller may be supplied with
:c:func:`FIRKodeSetAdaptController`.

The step size proposed by the controller is multiplied by a safety factor
(default 0.9, see :c:func:`FIRKodeSetSafetyFactor`) which, following RADAU5,
is further reduced according to the number of Newton iterations
:math:`k_{\text{nls}}` of the last nonlinear solve,

.. math::
   \mathrm{safety} \leftarrow \min\left( \mathrm{safety},\;
   \frac{2 m + 1}{2 m + k_{\text{nls}}} \right) \, ,

where :math:`m` is the maximum number of Newton iterations; this
adjustment may be disabled with :c:func:`FIRKodeSetNewtonCountSafety`. The
ratio :math:`\eta = h_{\text{new}} / h` is then limited to the interval
:math:`[\eta_{\min}, \eta_{\max}]`, where :math:`\eta_{\min}` (default 0.2)
is set with :c:func:`FIRKodeSetMinReduction` and :math:`\eta_{\max}` is the
maximum growth factor (default 8, see :c:func:`FIRKodeSetMaxGrowth`), except
on the first step where a much larger growth (default :math:`10^4`, see
:c:func:`FIRKodeSetMaxFirstGrowth`) is allowed, and after repeated error test
failures where the growth is limited to a factor :math:`\eta_{\max,f}`
(default 0.3, see :c:func:`FIRKodeSetMaxEFailGrowth` and
:c:func:`FIRKodeSetSmallNumEFails`). Finally, to avoid unnecessary
refactorizations of :math:`M - \gamma J`, a successful step keeps the current
step size whenever :math:`\eta` lies within a dead band
:math:`[\eta_{lb}, \eta_{ub}]` (default :math:`[1.0, 1.2]`, see
:c:func:`FIRKodeSetFixedStepBounds`). The step size is also limited by the
user-specified minimum and maximum step sizes (see :c:func:`FIRKodeSetMinStep`
and :c:func:`FIRKodeSetMaxStep`).

When a fixed step size is requested with :c:func:`FIRKodeSetFixedStep`, the
error estimate and the controller are bypassed and a nonlinear solver failure
immediately terminates the integration.

.. _FIRKODE.Mathematics.initial_step:

Initial step size
=================

Unless the user supplies an initial step size with
:c:func:`FIRKodeSetInitStep`, FIRKODE estimates the initial step size with
the heuristic used by CVODE. The step size :math:`h_0` is chosen so that
:math:`h_0^2 \, \|\ddot{y}(t_0)\|_{\text{WRMS}} / 2 \approx 1`, where the
second derivative is approximated by a difference quotient of
:math:`\dot{y} = M^{-1} f(t, y)` over a trial step; the trial is repeated a
few times and the result is limited to a fraction of the distance to the
first output time.

.. _FIRKODE.Mathematics.dense:

Dense output
============

Within each step, the collocation polynomial :math:`u(t)` of degree
:math:`s` that interpolates :math:`y_n` at :math:`t_n` and the stage values
:math:`Y_i` at :math:`t_n + c_i h` is available for output at any
:math:`t` in the step. Writing :math:`\theta = (t - t_n)/h`,

.. math::
   u(t_n + \theta h) = y_n + \sum_{j=1}^{s} L_j(\theta)\, Z_j \, ,
   :label: FIRKODE_dense

where :math:`L_j` is the Lagrange basis polynomial on the nodes
:math:`\{0, c_1, \ldots, c_s\}` with :math:`L_j(c_i) = \delta_{ij}` and
:math:`L_j(0) = 0`. The polynomial and its derivatives up to order
:math:`s` may be evaluated with :c:func:`FIRKodeGetDky` for
:math:`t` in the last step taken. The interpolation error of
:eq:`FIRKODE_dense` at interior points is of order :math:`s + 1` in
:math:`h`, which is lower than the classical order :math:`2s - 1` of the
method at the step end points. The same polynomial is used to predict the
stage increments of the next step (see :numref:`FIRKODE.Mathematics.Newton`)
and to locate roots (see :numref:`FIRKODE.Mathematics.rootfinding`).

.. _FIRKODE.Mathematics.mass:

Non-identity mass matrices
==========================

FIRKODE supports problems :eq:`FIRKODE_ivp` with a constant, nonsingular
mass matrix :math:`M \ne I`. The mass matrix enters the stage equations
:eq:`FIRKODE_stages`, the Newton matrix :eq:`FIRKODE_Newton`, the
preconditioner :eq:`FIRKODE_precond`, and the error estimate
:eq:`FIRKODE_errest` only through matrix-vector products :math:`M v` and
through the block matrices :math:`M - \gamma J`. The user supplies the mass
matrix either as a :c:type:`SUNMatrix` object of the same type as the
Jacobian matrix, together with a routine that fills it, or, for a
matrix-free linear solver, as a matrix-vector product routine (see
:c:func:`FIRKodeSetMassLinearSolver`). Linear solves with :math:`M` itself
are required only to compute :math:`\dot{y} = M^{-1} f(t, y)` in the initial
step size estimate and in the rootfinding module; these are performed with
the mass linear solver supplied by the user. The mass matrix is constructed
and factored once, at the beginning of the integration. Time-dependent mass
matrices are not supported.

.. _FIRKODE.Mathematics.rootfinding:

Rootfinding
===========

The FIRKODE solver has been augmented to include a rootfinding feature. This
means that, while integrating the Initial Value Problem :eq:`FIRKODE_ivp`,
FIRKODE can also find the roots of a set of user-defined functions
:math:`g_i(t,y)` that depend both on :math:`t` and on the solution vector
:math:`y = y(t)`. The number of these root functions is arbitrary, and if
more than one :math:`g_i` is found to have a root in any given interval, the
various root locations are found and reported in the order that they occur
on the :math:`t` axis, in the direction of integration.

Generally, this rootfinding feature finds only roots of odd multiplicity,
corresponding to changes in sign of :math:`g_i(t,y(t))`, denoted
:math:`g_i(t)` for short. If a user root function has a root of even
multiplicity (no sign change), it will probably be missed by FIRKODE. If such
a root is desired, the user should reformulate the root function so that it
changes sign at the desired root.

The basic scheme used is to check for sign changes of any :math:`g_i(t)`
over each time step taken, and then (when a sign change is found) to hone in
on the root(s) with a modified secant method :cite:p:`HeSh:80`. In addition,
each time :math:`g` is computed, FIRKODE checks to see if :math:`g_i(t) = 0`
exactly, and if so it reports this as a root. However, if an exact zero of
any :math:`g_i` is found at a point :math:`t`, FIRKODE computes :math:`g` at
:math:`t + \delta` for a small increment :math:`\delta`, slightly further in
the direction of integration, and if any :math:`g_i(t + \delta)=0` also,
FIRKODE stops and reports an error. This way, each time FIRKODE takes a time
step, it is guaranteed that the values of all :math:`g_i` are nonzero at some
past value of :math:`t`, beyond which a search for roots is to be done.

At any given time in the course of the time-stepping, after suitable checking
and adjusting has been done, FIRKODE has an interval :math:`(t_{lo},t_{hi}]`
in which roots of the :math:`g_i(t)` are to be sought, such that
:math:`t_{hi}` is further ahead in the direction of integration, and all
:math:`g_i(t_{lo}) \neq 0`. The endpoint :math:`t_{hi}` is either
:math:`t_n`, the end of the time step last taken, or the next requested
output time :math:`t_{\text{out}}` if this comes sooner. The endpoint
:math:`t_{lo}` is either :math:`t_{n-1}`, the last output time
:math:`t_{\text{out}}` (if this occurred within the last step), or the last
root location (if a root was just located within this step), possibly
adjusted slightly toward :math:`t_n` if an exact zero was found. The
algorithm checks :math:`g_i` at :math:`t_{hi}` for zeros and for sign changes
in :math:`(t_{lo},t_{hi})`. If no sign changes were found, then either a root
is reported (if some :math:`g_i(t_{hi}) = 0`) or we proceed to the next time
interval (starting at :math:`t_{hi}`). If one or more sign changes were
found, then a loop is entered to locate the root to within a rather tight
tolerance, given by

.. math:: \tau = 100 * U * (|t_n| + |h|)~~~ (U = \mbox{unit roundoff}) ~.

Whenever sign changes are seen in two or more root functions, the one deemed
most likely to have its root occur first is the one with the largest value of
:math:`|g_i(t_{hi})|/|g_i(t_{hi}) - g_i(t_{lo})|`, corresponding to the
closest to :math:`t_{lo}` of the secant method values. At each pass through
the loop, a new value :math:`t_{mid}` is set, strictly within the search
interval, and the values of :math:`g_i(t_{mid})` are checked. Then either
:math:`t_{lo}` or :math:`t_{hi}` is reset to :math:`t_{mid}` according to
which subinterval is found to include the sign change. If there is none in
:math:`(t_{lo},t_{mid})` but some :math:`g_i(t_{mid}) = 0`, then that root is
reported. The loop continues until :math:`|t_{hi}-t_{lo}| < \tau`, and then
the reported root location is :math:`t_{hi}`.

In the loop to locate the root of :math:`g_i(t)`, the formula for
:math:`t_{mid}` is

.. math::

   t_{mid} = t_{hi} - (t_{hi} - t_{lo})
                g_i(t_{hi}) / [g_i(t_{hi}) - \alpha g_i(t_{lo})] ~,

where :math:`\alpha` is a weight parameter. On the first two passes through
the loop, :math:`\alpha` is set to :math:`1`, making :math:`t_{mid}` the
secant method value. Thereafter, :math:`\alpha` is reset according to the
side of the subinterval (low vs. high, i.e., toward :math:`t_{lo}` vs. toward
:math:`t_{hi}`) in which the sign change was found in the previous two
passes. If the two sides were opposite, :math:`\alpha` is set to 1. If the
two sides were the same, :math:`\alpha` is halved (if on the low side) or
doubled (if on the high side). The value of :math:`t_{mid}` is closer to
:math:`t_{lo}` when :math:`\alpha < 1` and closer to :math:`t_{hi}` when
:math:`\alpha > 1`. If the above value of :math:`t_{mid}` is within
:math:`\tau/2` of :math:`t_{lo}` or :math:`t_{hi}`, it is adjusted inward,
such that its fractional distance from the endpoint (relative to the interval
size) is between .1 and .5 (.5 being the midpoint), and the actual distance
from the endpoint is at least :math:`\tau/2`.

The values of :math:`y(t)` needed to evaluate :math:`g` at interior points of
a step are obtained from the dense output :eq:`FIRKODE_dense`.
