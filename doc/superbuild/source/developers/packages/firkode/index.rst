..
   Author(s): Yifan Hu @ UMBC
   -----------------------------------------------------------------------------
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
   -----------------------------------------------------------------------------

.. _FIRKODE.Dev:

FIRKODE
=======

FIRKODE integrates :math:`M y' = f(t, y)` with the three-stage Radau IIA
method. Its design document, ``src/firkode/DESIGN.md``, is the specification:
the listings of its Chapter 8 are the control flow below, and its Table 11.3
is the hazard checklist whose rows the unit tests enforce
(``test/unit_tests/firkode/README.md`` maps each row to its test). The source
files follow Table 8.9 of the design: ``firkode.c`` (the evolve loop and the
initial step), ``firkode_step.c`` (one step attempt), ``firkode_nls.c`` (the
Newton driver), ``firkode_adapt.c`` (step control and the reuse policies),
``firkode_dense.c`` (predictor and dense output), ``firkode_stack.c`` (the
stage stack), ``firkode_ls.c`` (the linear solver and mass matrix
interfaces), ``firkode_root.c``, ``firkode_io.c``, ``firkode_cli.c`` and the
generated ``firkode_tables.c``.

The evolve loop
---------------

.. digraph:: firkode_evolve

   node [shape=box]
   splines=ortho

   init    [label="First call after Init, ReInit or Reset:\lcheck inputs, load the weights,\linitialize the solvers, set up M,\levaluate f(t0, y0), choose h0,\lcheck g at t0\l"]
   stop    [label="Later calls: roots in the last step?\lat tstop? past tout (interpolate)?\lONE_STEP output pending?\l"]
   ewt     [label="Weights at y_n (once per step)\lstep and accuracy limits\lf(t_n, y_n) if not current\l"]
   step    [label="firkStep: one step from t_n\l(see below)\l"]
   fail    [label="Return the failure code\l(firkHandleFailure)\l", style=filled, fillcolor=tomato1]
   roots   [label="Roots in (t_n, t_n + h)?\lreturn FIRK_ROOT_RETURN\l"]
   tstop   [label="At tstop? return FIRK_TSTOP_RETURN\lnext step would pass tstop? clip h\l"]
   tout    [label="NORMAL: past tout? interpolate\lwith the collocation polynomial\lONE_STEP: return y_n\l", style=filled, fillcolor=lightskyblue]

   init  -> stop [label="not the first call"]
   init  -> ewt  [label="first call"]
   stop  -> ewt  [label="keep integrating"]
   ewt   -> step
   step  -> fail [label="failure"]
   step  -> roots [label="success"]
   roots -> tstop [label="no root"]
   tstop -> tout  [label="not at tstop"]
   tout  -> ewt   [label="tout not reached"]

One step attempt
----------------

.. digraph:: firkode_step_attempt

   node [shape=box]
   splines=ortho

   policy  [label="Reuse policy (RADAU5 or PERIODIC):\lsetup? new Jacobian?\l"]
   predict [label="Predict the stage increments Z\lfrom the last collocation polynomial\l"]
   setup   [label="lsetup: M - gamma J, block solver setup\l"]
   newton  [label="Newton iteration on the stacked system\l(I_s (x) M - h A (x) J) delta = -G(Z)\lstacked FGMRES, blocks preconditioned\lby the user's solver; RADAU5's tests\l"]
   nfail   [label="Recoverable failure:\lh <- h x factor of the failure,\lretry (PREDICT_AGAIN) or give up\l", style=filled, fillcolor=tomato1]
   est     [label="Error estimate\l(M - gamma J)^-1 [gamma f(t_n, y_n) + M sum e_i Z_i]\lrefiltered on the first step\lor after a rejection\l"]
   adapt   [label="firkAdapt: the next step size\l(controller, Newton count safety,\ldead band, bounds)\l"]
   etest   [label="dsm <= 1?"]
   efail   [label="Error test failed:\lh <- h x eta, retry (TRY_AGAIN)\lor give up\l", style=filled, fillcolor=tomato1]
   done    [label="firkCompleteStep: controller history,\ly_n+1 = y_n + Z_s, swap Z and Zprev,\lt_n <- t_n + h\l", style=filled, fillcolor=lightskyblue]

   policy  -> predict
   predict -> setup   [label="setup needed"]
   predict -> newton  [label="no setup"]
   setup   -> newton
   newton  -> nfail   [label="failure"]
   newton  -> est     [label="converged (adaptive)"]
   newton  -> done    [label="converged (fixed step)"]
   nfail   -> policy  [label="retry"]
   est     -> adapt
   adapt   -> etest
   etest   -> efail   [label="no"]
   etest   -> done    [label="yes"]
   efail   -> policy  [label="retry"]

Where to look
-------------

- The Newton test, its rate estimate and step factors: ``firkode_nls.c`` and
  design Section 4.4; the automatic Newton tolerance is RADAU5's rule.
- The tolerances of the inner solves and the units conversion with a mass
  matrix: ``firkode_ls.c`` (``firkLsSolveBlock``, ``firkLsSolveStacked``) and
  design Section 5.10.
- The comparison with ``radau5.f``: ``test/unit_tests/firkode/radau5`` and
  ``firk_test_radau5_steps``.
- The formal verification of the tables and the stability function:
  ``verification/lean``.
