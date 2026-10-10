# FIRKODE v0.1.0 Design Document (October 6, 2026)

FIRKODE is a proposed SUNDIALS package for stiff initial value problems

$$
M \dot y = f(t, y), \quad y(t_0) = y_0, \quad y \in \mathbb{R}^N,
$$

integrated by the fully implicit Radau IIA Runge–Kutta methods.
This document is written for two groups:
the SUNDIALS developers and collaborators who decide whether the package belongs in the suite
(Chapters 1–3 and 5–7), and the implementer who builds it milestone by milestone (Chapters 4 and
8–11). The user guide is a deliverable of the milestones and is not duplicated here.

Status: a draft proposal based on SUNDIALS 7.9.0. The only pre-existing material is the supporting
tooling named in Appendix A (four stability figures and the script that draws them) and a Lean
pilot in `verification/lean`, cited only as evidence that the certification approach of
Milestone M7 is feasible.

Every claim that is not self-evident carries one of five evidence tags, so that the reader
can see what rests on a theorem, what on a computation, and what on nothing yet:

| Evidence tag         | Meaning                                                                                                |
| -------------------- | ------------------------------------------------------------------------------------------------------ |
| `[proved]`           | a theorem of the literature or a kernel-checked Lean statement                                         |
| `[derived]`          | computed from the method tables by a procedure described here and reproducible by the table generator  |
| `[cited]`            | a published result, with its reference                                                                 |
| `[fact: path]`       | read from the SUNDIALS sources at the given path (and line, where stable)                              |
| `[to measure]`       | no number exists; a measurement protocol and a milestone gate are given instead                        |

## 1. Summary

FIRKODE integrates $M \dot y = f(t, y)$ with constant nonsingular $M$ by the $s$-stage Radau IIA
collocation methods, which are A-stable, L-stable and stiffly accurate, with classical order $2s - 1$
and stage order $s$ `[proved]`. Version 0.1.0 ships the 3-stage method of order 5, the method of
Hairer and Wanner's RADAU5 [1996HW, 1999HW], with RADAU5's simplified Newton iteration, convergence
test, filtered error estimate and step-control heuristics. Later versions extend the same code path to
$1 \le s \le 9$ (orders 1 to 17), add RADAU's in situ order selection, and open the stage solver to
expert users.

### 1.1 What is new relative to RADAU5

RADAU5 transforms the $sN \times sN$ Newton
matrix to block-diagonal form and factors one real and $\left\lfloor \frac{s}{2} \right\rfloor$ complex $N \times N$
matrices. FIRKODE instead solves the stacked Newton system with flexible GMRES [1993S],
right-preconditioned by $I_s \otimes \left(M - h\gamma_0 J \right)$, where $\gamma_0$ is the reciprocal of the
real eigenvalue of $A^{-1}$. The preconditioner needs one real factorization per setup and is applied
by $s$ independent solves with the user's ordinary $N \times N$ `SUNLinearSolver`, or, when that
solver is iterative, by $s$ applications of its preconditioner. Every linear
solver and preconditioner written for CVODE's or ARKODE's $M - \gamma J$ therefore works unchanged,
including matrix-free Krylov solvers, algebraic multigrid and batched GPU solvers, on any `N_Vector`.
The price is more block solves per Newton iteration than a direct complex factorization needs; Chapter
7 quantifies the break-even.

### 1.2 What FIRKODE adds to SUNDIALS

SUNDIALS has no fully implicit Runge–Kutta method. CVODE's
BDF methods are A-stable only up to order 2 [1963D], and ARKODE's diagonally implicit methods have stage
order at most 2 [2019KC], so neither offers a high-order L-stable method free of stage-order reduction on
stiff problems. Radau IIA is the classical answer to that gap, and SUNDIALS' shared vectors, linear
solvers, controllers, logging and profiling are what RADAU5 never had.

### 1.3 Design decisions

| ID | Decision                         | Verdict                                                                                                 | Principle | Revisit when                                           |
| -- | -------------------------------- | ------------------------------------------------------------------------------------------------------- | --------- | ------------------------------------------------------ |
| D1 | Position in SUNDIALS             | A standalone package, not an ARKODE stepper                                                             | P1, P3    | ARKODE gains a stage-coupled linear-solver interface   |
| D2 | Stage linear algebra             | Real FGMRES with the single-shift block preconditioner; no complex arithmetic                           | P4        | measured $k_{lin}$ for $s \ge 5$ exceeds about $2s$    |
| D3 | Nonlinear iteration              | An internal simplified-Newton loop with RADAU5's test, not a `SUNNonlinearSolver`                       | P4        | a semantic conflict with the hook interface disappears |
| D4 | Stage storage                    | An `NVECTOR_MANYVECTOR` of $s$ clones of the user's vector; `MPIMANYVECTOR` under MPI                   | P1        | GPU launch overhead is measured to dominate            |
| D5 | Step control                     | `SUNAdaptController` with $p = s$ plus RADAU5's safety, limits and dead band applied by FIRKODE         | P1        | never; only the defaults                               |
| D6 | Method tables                    | Generated decimal literals, checked in CI, optionally certified in Lean                                 | P1        | adding $s > 9$                                         |
| D7 | Order                            | v0.1.0: $s = 3$; v0.2: any $1 \le s \le 9$ in EXPERT mode; v0.3: in situ odd order in BEGINNER mode     | P2        | calibration data from v0.2                             |
| D8 | Modes                            | BEGINNER exposes the essentials and rejects every other setter; EXPERT unlocks everything               | P2        | user feedback                                          |
| D9 | Jacobian and factorization reuse | RADAU5's policy (BEGINNER default) or CVODE's periodic policy; the step driver decides, FIRKLS executes | P4        | M4 shows the better policy depends on the block solver |

*Table 1.1. Decision record; D-numbers refer to Chapter 5. The principles P1–P4 are stated in
Chapter 2.*

### 1.4 Scope

#### In v0.1.0

Milestones M0–M3 of Chapter 10, that is a RADAU5-equivalent integrator with
SUNDIALS' interfaces, rootfinding, a constant mass matrix, logging, profiling, command-line options,
and a validation suite that later versions run unchanged.

#### Not in v0.1.0

$s \ne 3$, order selection,
EXPERT hooks, SUNStepper glue, time-dependent or singular $M$, sensitivities, relaxation, resizing,
and language bindings. Chapter 7 states the costs and risks, and what is asked of the SUNDIALS
maintainers.

## 2. Guiding principles

Four principles decide every trade-off in this document. Each decision in Chapter 5 and each cut in
Chapter 10 names the principle it follows.

- **P1: Plan for the long term.** FIRKODE is unfunded work. Maintainability and a bounded maintenance
  surface take precedence over features. Shared SUNDIALS components are reused wherever they fit, and
  nothing is written that a shared component already provides. Every milestone is useful on its own
  and small enough for one person to finish and one reviewer to read.
- **P2: Limit the scope of v0.1.0.** No feature is added because it looks attractive. Version 0.1.0 is
  a RADAU5-equivalent integrator for $s = 3$. The method family, order selection, expert hooks and
  SUNStepper glue are versioned follow-ups with their own acceptance criteria, and BEGINNER mode
  exposes only what a user needs who has no deep understanding of fully implicit methods.
- **P3: Do not optimize for use as an inner step.** FIRK methods are not tuned for service as the
  implicit stage of multirate, splitting or other advanced methods. The SUNStepper interface is the
  last optional milestone, the argument for a standalone package does not rest on composability, and
  no example uses FIRKODE inside another integrator.
- **P4: Nonlinear-solver performance is critical.** The cost of a fully implicit method is the cost of
  its Newton iteration. The stacked iteration, the block-preconditioned stage solve, the Jacobian and
  factorization reuse policies, and the measurement gates on Krylov iteration counts are the centre of
  the architecture and of the cost model. The expert hooks of v0.4 exist to let users change the stage
  solver.

## 3. Motivation

### 3.1 The gap in SUNDIALS

SUNDIALS has no integrator that is at the same time A-stable at order above 2, L-stable,
algebraically stable, of stage order above 2, and one-step. Table 3.1 lists the suite's integrators
and the properties that matter for stiff problems.

| Integrator              | Methods and orders                                                                 | A-stable orders                                                        | L-stable orders           | Stage order                    | Stiffly accurate   | Mass matrix                 | Dense output                                                          | Variable order   |
| ----------------------- | ---------------------------------------------------------------------------------- | ---------------------------------------------------------------------- | ------------------------- | ------------------------------ | ------------------ | --------------------------- | --------------------------------------------------------------------- | ---------------- |
| CVODE, Adams            | Adams–Moulton 1–12 `[fact: src/cvode/cvode_impl.h]`                                | 1, 2 `[cited]`                                                         | 1                         | not applicable                 | no                 | no                          | Nordsieck polynomial, degree $q$                                      | yes              |
| CVODE, BDF              | BDF 1–5 `[fact: src/cvode/cvode_impl.h]`                                           | 1, 2; $A(\alpha)$ with $\alpha = 86^\circ, 73^\circ, 52^\circ$ for 3–5 [1996HW, §V.2] | 1, 2                      | not applicable                 | no (multistep)     | no                          | Nordsieck polynomial, degree $q$                                      | yes              |
| IDA                     | BDF 1–5 on $F(t, y, \dot y) = 0$ `[fact: src/ida/ida_impl.h]`                      | as CVODE BDF                                                           | 1, 2                      | not applicable                 | no                 | implicit form, index-1 DAEs | divided differences, degree $q$                                       | yes              |
| ARKStep, DIRK/ESDIRK    | tables of orders 2–5 `[fact: include/arkode/arkode_butcher_dirk.h]`                | tables of every order 2–5                                              | tables of every order 2–5 | 1 (SDIRK), 2 (ESDIRK) [2019KC] | many tables        | yes (ARKStep only)          | Hermite or Lagrange, degree $\le 5$ `[fact: include/arkode/arkode.h]` | no               |
| ARKStep, IMEX           | ARK pairs of orders 2–5                                                            | implicit part as above                                                 | as above                  | as above                       | as above           | yes                         | as above                                                              | no               |
| ERKStep                 | explicit tables of orders 2–9 `[fact: include/arkode/arkode_butcher_erk.h]`        | none                                                                   | none                      | 1                              | no                 | no                          | as above                                                              | no               |
| LSRKStep                | RKC and RKL (order 2), SSP (orders 2–4) `[fact: include/arkode/arkode_lsrkstep.h]` | none (explicit, extended real stability)                               | none                      | 1                              | no                 | no                          | as above                                                              | no               |
| SPRKStep                | symplectic partitioned, orders 1–10 `[fact: include/arkode/arkode_sprk.h]`         | none                                                                   | none                      | not applicable                 | no                 | no                          | as above                                                              | no               |
| MRIStep                 | MRI-GARK 1–4, IMEX-MRI-SR 2–4, MERK 2–5 `[fact: include/arkode/arkode_mristep.h]`  | slow implicit part as DIRK                                             | as DIRK                   | as DIRK                        | as DIRK            | no                          | as above                                                              | no               |
| FIRKODE (proposed)      | Radau IIA, $s = 3$ (order 5) in v0.1.0; $1 \le s \le 9$ (orders 1–17) from v0.2    | all `[proved]`                                                         | all `[proved]`            | $s$ `[proved]`                 | all `[proved]`     | yes, constant nonsingular   | collocation polynomial, degree $s$                                    | v0.3             |

*Table 3.1. Integrators of SUNDIALS 7.9.0 and the properties that matter for stiff problems. The*
$A(\alpha)$ *angles are rounded.*

Three consequences for stiff problems follow from the table.

- **Oscillatory stiff spectra.** When $J$ has eigenvalues of large modulus near the imaginary axis,
  for example in wave-dominated or weakly damped mechanical problems, BDF orders above 2 lose
  stability and CVODE must reduce its order. Radau IIA is A-stable at every order `[proved]`, so its
  order is limited by accuracy alone.
- **Stage-order reduction.** On stiff problems the observed order of a Runge–Kutta method is governed
  by its stage order rather than its classical order [1974PR, 1985FSU]. A DIRK method with stage order
  1 or 2 then converges at order about 2 to 3 whatever its classical order. The remedy under study in
  the DIRK literature, high *weak* stage order [2020KSSZ, 2023BKSS], is defined for linear problems
  with time-dependent forcing and recovers the classical order there. Radau IIA has stage order $s$
  for nonlinear problems too, so the 3-stage method keeps order at least 3, and typically more, on the
  same problems; on parabolic PDEs the order depends on the boundary data (Table 4.3). Section 4.2
  states the results precisely, and the validation plan measures them on the Prothero–Robinson family
  and on a heat equation with time-dependent boundary data.
- **Nonlinear stability.** The stability theory for stiff *nonlinear* problems, those with a one-sided
  Lipschitz condition $\langle f(t, y) - f(t, z), y - z \rangle \le \nu \lVert y - z \rVert^2$, rests
  on B-stability [1975B] and its algebraic criterion $BA + A^T B - b b^T \succeq 0$,
  $B = \mathrm{diag}(b)$ [1979BB, 1984DV]. Radau IIA is algebraically stable [1996HW, §IV.12], so for
  $\nu \le 0$ its numerical solutions are contractive for every step size, and it is B-convergent
  [1981FSU] of order $s$ on nonlinear dissipative problems [1985FSU, 1987BH] and of order $s + 1$ on
  semilinear problems [1986BHV] `[cited]`. The multistep analogue, G-stability, is equivalent to
  A-stability [1996HW, Ch. V], so BDF above order 2 has no such property, and a DIRK table has it only
  by construction [1984DV].

### 3.2 Why fully implicit Runge–Kutta? Why Radau IIA?

Among implicit Runge–Kutta families, Radau IIA alone combines L-stability, stiff accuracy, stage
order $s$ and classical order $2s - 1$ [1996HW, §IV.5]. Table 3.2 lists the alternatives and what
rules each out as the primary method.

| Family                 | Order                | Stage order | A-stable               | L-stable                          | Algebraically stable               | Stiffly accurate       | Stage coupling                          | Reason not chosen                                                            |
| ---------------------- | -------------------- | ----------- | ---------------------- | --------------------------------- | ---------------------------------- | ---------------------- | --------------------------------------- | ---------------------------------------------------------------------------- |
| Radau IIA              | $2s - 1$             | $s$         | yes                    | yes                               | yes                                | yes                    | full, complex eigenvalues of $A$        | chosen                                                                       |
| Radau IA               | $2s - 1$             | $s$         | yes                    | yes                               | yes                                | no                     | full                                    | not stiffly accurate; the update needs a quadrature of the stage derivatives |
| Gauss                  | $2s$                 | $s$         | yes                    | no, $\lvert R(\infty) \rvert = 1$ | yes                                | no                     | full                                    | no damping of stiff components; not stiffly accurate                         |
| Lobatto IIIC           | $2s - 2$             | $s - 1$     | yes                    | yes                               | yes                                | yes                    | full                                    | one order and one stage order less than Radau IIA for the same $s$           |
| SDIRK, ESDIRK          | $\le 5$ in practice  | 1, 2        | by construction        | by construction                   | some tables                        | by construction        | diagonal, one real eigenvalue           | stage-order reduction [2019KC]; already in ARKODE                            |
| Singly implicit (SIRK) | $\le s + 1$ [1977NW] | $s$         | for suitable $\lambda$ | for suitable $\lambda$            | not considered                     | for suitable $\lambda$ | one real eigenvalue of multiplicity $s$ | order barrier $s + 1$; optional later milestone (M10)                        |
| Rosenbrock, W-methods  | $\le 4$ in practice  | low         | by construction        | by construction                   | not applicable (linearly implicit) | no                     | linearly implicit                       | exact Jacobian every step, or order reduction with inexact $J$               |

*Table 3.2. Implicit Runge–Kutta families considered. Properties are classical results [1996HW, Ch.
IV; 1984DV]; the Radau IIA row is machine-checked for* $s \le 9$ *`[proved]`, except algebraic
stability, which Section 11.5 lists as open.*

Two classes lie outside the table. Exponential integrators [2010HO] treat semilinear stiff PDEs
through the matrix exponential; SUNDIALS has none, the MERK coupling tables of MRIStep being the only
exponential element in the suite `[fact: include/arkode/arkode_mristep.h]`. Spectral deferred
correction [2000DGR] iterates on a collocation system and, on the Radau nodes, converges to the Radau
IIA solution; its stage-parallel sweeps reappear in Section 5.2 as candidate preconditioners.

A one-step method brings further practical gains over the suite's multistep integrators: no startup
phase, step changes at no cost, and a collocation polynomial that serves dense output, rootfinding and
the Newton predictor without extra work (Section 4.9). The honest counterweight appears in Figure 3 of
Appendix A: the damping of stiff components per step, $\lvert R(z) \rvert \sim \frac{s}{\lvert z \rvert}$ as
$\lvert z \rvert \to \infty$, weakens as $s$ grows, which is one reason the default order range of
RADAU, and of FIRKODE from v0.3, stops at $s = 7$.

RADAU5 has a thirty-year production record and is the reference against which stiff solvers are still
compared [1996HW, §IV.10; 2025ESR]. All adaptive Radau codes in Chapter 6 keep its core: simplified
Newton with one Jacobian per step, a filtered embedded error estimate and Gustafsson's predictive
controller [1994G]. FIRKODE keeps that core and replaces the complex transformation.

### 3.3 Why now?

Two developments make a fully implicit package natural today, where it was not when CVODE and ARKODE
were designed.

First, the shared infrastructure that a block-preconditioned stage solver needs now exists in SUNDIALS.
Table 3.3 lists what FIRKODE reuses and the release that introduced it.

| Component                                                                                                    | Release `[fact: CHANGELOG.md]` | FIRKODE uses it for                                                 | What would otherwise be written       |
| ------------------------------------------------------------------------------------------------------------ | ------------------------------ | ------------------------------------------------------------------- | ------------------------------------- |
| `NVECTOR_MANYVECTOR`, `MPIMANYVECTOR`                                                                        | 5.0.0                          | the stacked vector of $s$ stage increments on any backend           | a stage multivector for every backend |
| `SUNLinSol_SPFGMR`                                                                                           | 3.0.0                          | the Krylov solver for the $sN$ Newton system                        | a flexible GMRES on stacked vectors   |
| `SUNContext`, `SUNProfiler`                                                                                  | 6.0.0                          | error handlers, profiling regions                                   | ad hoc error and timing code          |
| `SUNLogger`                                                                                                  | 6.2.0                          | per-step, per-iteration records in ARKODE's format                  | a private log format                  |
| `SUNAdaptController`                                                                                         | 6.7.0                          | step control, with every built-in controller selectable             | the controllers of Section 4.7        |
| command-line option helpers (`SetOptions`)                                                                   | 7.5.0                          | `firkode.key value` parsing, as ARKODE                              | an option parser                      |
| `SUNDomEigEstimator`                                                                                         | 7.5.0                          | spectral information for the v0.4 hooks                             | a power iteration                     |
| `SUNStepper`                                                                                                 | 7.2.0                          | the optional glue of M8 only (P3)                                   | nothing in v0.1.0                     |
| batched block solvers: `SUNLinSol_KokkosDense`, `SUNLinSol_GinkgoBatch` and the other block-diagonal solvers | 6.4.0, 7.5.0                   | the block $M - \gamma J$ for many independent systems (Section 3.3) | a batched stage solver                |
| `SUNAdjointStepper`                                                                                          | 7.3.0                          | the optional discrete adjoint of M13 only                           | nothing in v0.1.0                     |

*Table 3.3. Shared SUNDIALS infrastructure reused by FIRKODE.*

Second, the solver technology that makes fully implicit methods competitive at scale is recent and is
not available in any general-purpose C library: block preconditioners with mesh- and step-independent
bounds [2006SMN, 2007MNS, 2024GO], stage-parallel Radau IIA at the strong-scaling limit [2024M+],
monolithic multigrid for the stage system [2005VV, 2024K], real preconditioners for the complex
conjugate pairs [2022SKPD], and automated Radau IIA for finite elements in Irksome [2021FKM, 2025KM].
A fully adaptive Radau code in Julia reports gains over RADAU5 at tight tolerances [2025ESR]. FIRKODE
brings this class of method to the vectors, linear solvers and preconditioners that SUNDIALS users
already have.

Three further developments shape the design. Many independent stiff systems, one per cell of a
reacting-flow simulation, are now integrated on GPUs with batched block solvers, in SUNDIALS itself
[2025B+] and in dedicated kinetics codes [2018SAN]; FIRKODE applies its preconditioner block by block,
so the suite's batched `SUNLinearSolver` implementations serve it unchanged, and the added cost over
CVODE is the $s$ RHS evaluations and $s k_{lin}$ block solves per Newton iteration of Table 7.1.
Monolithic multigrid for fully implicit Runge–Kutta discretizations of incompressible flow [2023AMF],
an index-2 DAE, shows the demand for the DAE extension that Milestone M9 leaves open. And the
spectral-deferred-correction literature has produced analysed stage-parallel diagonal preconditioners
on the Radau nodes [2018S, 2025CLGR] and an exact low-rank correction for stage-parallel solves
[2026DM], which Section 5.2 lists as EXPERT stage solvers.

### 3.4 Why this design?

Each decision below is argued in Chapter 5; this section states the thread that connects them.

- A standalone package, because the stage-coupled linear algebra and RADAU5's coupling of Newton rate,
  step control and order control do not fit ARKODE's per-stage interfaces (D1).
- One real shift and flexible GMRES, because SUNDIALS has no stable complex types and because the user's
  existing $N \times N$ solver and preconditioner are the asset to protect (D2).
- An internal Newton loop with RADAU5's three-way test, because the test's step-size recommendation
  has no channel through the `SUNNonlinearSolver` interface (D3).
- A ManyVector stack, because it costs no new vector code and runs on every backend (D4).
- RADAU5's heuristics generalized to any $s$ through `SUNAdaptController`, so that the suite's
  controllers are available and RADAU5's behavior is the BEGINNER default (D5).
- Generated, certified tables, because platform-independent literals are reviewable and provable
  (D6).

### 3.5 Scope and non-goals

#### In scope for v0.1.0

$M \dot y = f(t, y)$ with constant nonsingular $M$, forward and backward
integration, the 3-stage method, adaptive and fixed steps, rootfinding, dense output, the complete
CVODE and ARKODE linear-solver interface for the block $M - \gamma J$, logging, profiling, statistics,
command-line options and a validation suite.

#### Not in scope

- **Singular** $M$ **and DAEs.** Radau IIA is a classical DAE method up to index 3 [1989HLR], but
  consistent initialization and index-aware error scaling are a project of their own. IDA covers
  index-1 DAEs today. Milestone M9 sketches time-dependent $M$; singular $M$ remains open.
- **Resizing.** PDE codes with adaptive meshes change the size of the state between steps;
  `ARKodeResize` and `CVodeResizeHistory` are the suite's precedents
  `[fact: include/arkode/arkode.h, include/cvode/cvode.h]`. A one-step method carries no history
  across the resize beyond the last collocation polynomial, so the cost is a reallocation of the
  stage stack and one step with the trivial predictor. Deferred to the optional Milestone M11 (P2).
- **Relaxation.** Relaxation Runge–Kutta methods [2019K, 2020R+] rescale the step so that a convex
  functional is conserved exactly; ARKODE supports them
  `[fact: include/arkode/arkode.h, ARKodeSetRelaxFn]`. For Radau IIA the relaxation direction
  $h \sum_j b_j F_j$ equals $\sum_j d_j Z_j$ (Section 4.3), so no RHS evaluation is added, and
  algebraic stability already makes the unrelaxed method contractive (Section 4.2). Deferred to the
  optional Milestone M12 (P2).
- **Sensitivities and adjoints.** `SUNAdjointStepper` provides discrete adjoints for explicit
  Runge–Kutta methods `[fact: include/sundials/sundials_adjointstepper.h]`. The discrete adjoint of
  an implicit Runge–Kutta method is again a Runge–Kutta-type method with the transposed stage coupling
  [2006S]; it needs checkpoints of the stage increments and a block solver for the transposed system.
  Forward sensitivities are not planned. Deferred to the optional Milestone M13, after the API freeze
  (P2).
- **Constraints.** CVODE, ARKODE and IDA enforce inequality constraints by step rejection
  `[fact: include/cvode/cvode.h, CVodeSetConstraints]`. Like every Runge–Kutta method of order above
  one, Radau IIA does not preserve positivity unconditionally [2021NRK]; the step-rejection approach
  of [2021NRK] is the route if users ask. Not planned (P1).
- **Language bindings.** The C API is kept bindable (Section 9.4); SWIG-Fortran and `sundials4py`
  bindings are a follow-up once the API has settled.
- **Device-specific kernels.** Any `N_Vector` works through its own operations. No GPU-specific code.

### 3.6 Migrating to FIRKODE

The target user integrates a stiff system with CVODE or ARKStep today and either needs higher order
with L-stability, has an oscillatory stiff spectrum, works at tight tolerances, or needs a mass
matrix with a stiffly accurate method. The interface is similar to CVODE and ARKODE: the same
callbacks, the same linear-solver attachment, the same tolerance and stop-time semantics.

| CVODE call                                                                                                                                                                       | FIRKODE call                                           | Note                                                                                                                                         |
| -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------- |
| `CVodeCreate(CV_BDF, sunctx)`                                                                                                                                                    | `FIRKodeCreate(sunctx)`                                | one method family; `SetMaxOrd` has no analogue in v0.1.0                                                                                     |
| `CVodeInit(mem, f, t0, y0)`                                                                                                                                                      | `FIRKodeInit(mem, f, t0, y0)`                          | identical RHS signature                                                                                                                      |
| `CVodeSStolerances`, `SVtolerances`, `WFtolerances`                                                                                                                              | same names                                             | WRMS weights at $y_n$; no RADAU5 tolerance transformation                                                                                    |
| `CVodeSetLinearSolver(mem, LS, A)`                                                                                                                                               | `FIRKodeSetLinearSolver(mem, LS, A)`                   | the solver targets the block $M - \gamma J$; a Krylov `SUNLinearSolver` contributes only its preconditioner to the stage solve (Section 8.6) |
| `CVodeSetJacFn`, `SetJacTimes`, `SetPreconditioner`, `SetLinSysFn`, `SetEpsLin`, `SetLSetupFrequency`, `SetJacEvalFrequency`, `SetDeltaGammaMaxLSetup`, `SetDeltaGammaMaxBadJac` | same names; the last four are EXPERT options           | identical callback signatures; all evaluated at $(t_n, y_n)$ with $\gamma = h\gamma_0$                                                       |
| `CVodeSetMaxNonlinIters(3)`, `SetNonlinConvCoef(0.1)`                                                                                                                            | same names, defaults 7 and RADAU5's rule (Section 4.4) | EXPERT mode                                                                                                                                  |
| `CVodeSetNonlinearSolver`                                                                                                                                                        | none                                                   | the Newton iteration is internal (D3)                                                                                                        |
| `CVode(mem, tout, y, &t, CV_NORMAL)`                                                                                                                                             | `FIRKodeEvolve(mem, tout, y, &t, FIRK_NORMAL)`         | same `tstop`, root and one-step semantics                                                                                                    |
| `CVodeGetDky(mem, t, k, dky)`                                                                                                                                                    | `FIRKodeGetDky(mem, t, k, dky)`                        | $0 \le k \le s$, collocation polynomial of the last step                                                                                     |
| `CVodeRootInit`, `SetRootDirection`, `GetRootInfo`                                                                                                                               | same names                                             | ARKODE's Illinois rootfinder on the collocation polynomial                                                                                   |
| `CVodeSetStabLimDet`, `SetConstraints`, `SetUseIntegratorFusedKernels`, `CVodeResizeHistory`                                                                                     | none                                                   | no analogue                                                                                                                                  |
| `CVodeSetOptions(mem, "cvode", NULL, argc, argv)`                                                                                                                                | `FIRKodeSetOptions(mem, "firkode", NULL, argc, argv)`  | one key per scalar setter                                                                                                                    |
| `CV_*` return codes                                                                                                                                                              | `FIRK_*`                                               | CVODE numbering where the meaning is shared (Table 8.2)                                                                                      |

*Table 3.4. Migrating from CVODE.*

| ARKStep call                                                                                                                                                         | FIRKODE call                                                                 | Note                                                                                                                                                 |
| -------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------- |
| `ARKStepCreate(NULL, fi, t0, y0, sunctx)`                                                                                                                            | `FIRKodeCreate`, `FIRKodeInit`                                               | implicit RHS only                                                                                                                                    |
| `ARKodeSetMassLinearSolver(mem, LS, M, time_dep)`                                                                                                                    | `FIRKodeSetMassLinearSolver(mem, LS, M, SUNFALSE)`                           | `time_dep = SUNTRUE` is rejected in v0.1.0 (M9)                                                                                                      |
| `ARKodeSetMassFn`, `SetMassTimes`, `SetMassPreconditioner`, `SetMassEpsLin`                                                                                          | same names                                                                   | identical signatures; the `SetMassEpsLin` factor multiplies the norm of the right-hand side (Sections 4.8, 5.10), not $\epsilon_{\text{nls}}$        |
| `ARKodeSetAdaptController(mem, C)`, `SetAdaptControllerByName`                                                                                                       | same names                                                                   | EXPERT mode; the controller must be of type `SUN_ADAPTCONTROLLER_H` and receives $p = s$, not an embedding order; `ByName` does not accept `ImExGus` |
| `ARKodeSetSafetyFactor`, `SetMaxGrowth`, `SetMinReduction`, `SetFixedStepBounds`, `SetMaxFirstGrowth`, `SetMaxEFailGrowth`, `SetSmallNumEFails`, `SetMaxCFailGrowth` | same names                                                                   | EXPERT mode; same semantics                                                                                                                          |
| `ARKStepSetTables`, `ARKodeSetOrder`                                                                                                                                 | `FIRKodeSetTable` (v0.4), `FIRKodeSetOrder` (EXPERT; $s = 3$ only in v0.1.0) | `SetOrder(q)` picks the smallest $s$ with $2s - 1 \ge q$                                                                                             |
| `ARKodeSetInterpolantType`, `SetInterpolantDegree`                                                                                                                   | none                                                                         | the collocation polynomial is the interpolant                                                                                                        |
| `ARKodeSetNonlinCRDown`, `SetNonlinRDiv`                                                                                                                             | `FIRKodeSetNonlinDivergenceRate`                                             | RADAU5's test replaces the CVODE/ARKODE rate test                                                                                                    |
| `ARKodeResize`, `SetRelaxFn`, `SetPostprocessStageFn`, `SetStagePredictFn`, `SetDeduceImplicitRhs`                                                                   | none                                                                         | no analogue                                                                                                                                          |
| `ARKodeCreateSUNStepper`                                                                                                                                             | `FIRKodeCreateSUNStepper` (optional M8)                                      | not in v0.1.0 (P3)                                                                                                                                   |

*Table 3.5. Migrating from ARKStep.*

Inside the suite, Table 6.3 says which package to reach for; FIRKODE is the answer for stiff problems
with oscillatory spectra, tight tolerances, or a mass matrix that needs a stiffly accurate high-order
method.

## 4. Mathematical formulation

This chapter defines the algorithms.
Decisions and their alternatives are in Chapter 5;
the data structures and control flow for the algorithms are in Chapter 8.

### 4.1 Problem class

$$
M \dot y = f(t, y),
\quad y(t_0) = y_0,
\quad y \in \mathbb{R}^N.
$$

- $M$ is constant and nonsingular. It is the identity unless a mass matrix is attached.
- Integration runs forward or backward in $t$.
- $J = \frac{\partial f}{\partial y}$ is the Jacobian, $u$ the unit roundoff, and $\lVert \cdot \rVert$ the
  WRMS norm with the weights $w_i = \frac{1}{\mathrm{rtol} \lvert y_i \rvert + \mathrm{atol}_i}$,
  evaluated at the start of each step.

### 4.2 Radau IIA methods

The $s$-stage Radau IIA method is the collocation method on the nodes $0 < c_1 < \dots < c_s = 1$,
the zeros of $P_s(2c - 1) - P_{s-1}(2c - 1)$ with $P_k$ the Legendre polynomial of degree $k$
[1964Bb, 1969A]. Its coefficients satisfy $C(s)$,

$$
\sum_{j=1}^{s} a_{ij} c_j^{k-1} = \frac{c_i^k}{k},
\quad i, k = 1, \dots, s,
$$

and $b_j = a_{sj}$, so the method is stiffly accurate.

- Order $2s - 1$, from $B(2s-1)$, $C(s)$ and $D(s-1)$ [1964Ba, 1993HNW]; stage order $s$ `[proved]`.
- The stability function is the $(s-1, s)$ Padé approximant of $e^z$ `[proved]`. The method is
  A-stable and L-stable [1978WHN, 1996HW] `[proved]`; Figures 1–3 in Appendix A plot $R(z)$ for
  $1 \le s \le 9$.
- On stiff problems the observed order lies between the stage order and the classical order. On the
  Prothero–Robinson problem it is at least $s$, and $s + 1$ in the stiff limit, where a DIRK method
  is limited to about its stage order plus one [1974PR, 1985FSU; 1996HW, §IV.15].
- $s = 3$ is RADAU5 and the only method of v0.1.0. Later versions provide $1 \le s \le 9$, orders 1
  to 17; $s = 1$ is the implicit Euler method.
- The method is algebraically stable: $BA + A^T B - b b^T \succeq 0$ with $B = \mathrm{diag}(b)$ and
  $b_i > 0$ [1979BB; 1996HW, §IV.12] `[cited]`. Hence it is B-stable [1975B]: for $M = I$ and a
  one-sided Lipschitz constant $\nu \le 0$, two numerical solutions satisfy
  $\lVert y_{n+1} - \tilde y_{n+1} \rVert \le \lVert y_n - \tilde y_n \rVert$ in the Euclidean norm
  for every $h$ [1984DV]. It is B-convergent of order $s$ on nonlinear dissipative problems [1985FSU,
  1987BH] and of order $s + 1$ on semilinear problems whose stiffness sits in the linear part
  [1986BHV] `[cited]`. The criterion is a finite check on the table; Section 11.5 lists it as open.
- The $s$-stage method coincides with the discontinuous Galerkin method in time of degree $s - 1$
  with the right Radau quadrature, a correspondence attributed to [1974LR] and used by [2022SKPD];
  preconditioners analysed for DG in time therefore apply to the stage system `[cited]`.
- On parabolic PDEs the order is set by the boundary data rather than by the classical order
  [1986SVH, 1992OR]. For a strongly $A(\theta)$-stable method of stage order $q$ and order
  $p \ge q + 1$ on the linear heat equation, the error is $O(h^{\min(p,\, q + 1 + \beta)})$ under the
  regularity $y^{(q+1)}(t) \in D(\mathcal{A}^\beta)$ of the spatial operator $\mathcal{A}$, with
  $\beta$ determined by the boundary conditions [1993LO, Thm. 3.3]; semilinear and quasilinear
  problems behave alike [1993LO, Thm. 4.1; 1995LO], the quasilinear energy argument needing algebraic
  stability. Table 4.3 evaluates the bound for Radau IIA ($q = s$, $p = 2s - 1$): in the worst case
  the order is the stage order plus one, so $s = 3$ gives about $4.25$ on a heat equation with
  time-dependent Dirichlet data and the full 5 only with homogeneous Dirichlet data `[cited]`. The
  textbook account is [2003HV]; the validation plan measures the slopes (Table 11.4).

| Boundary data of the heat equation | $\beta$ [1993LO]         | Order of the $s$-stage Radau IIA method | $s = 3$           |
| ---------------------------------- | ------------------------ | --------------------------------------- | ----------------- |
| time-dependent Dirichlet           | $\frac{1}{4} - \epsilon$ | $\min(2s - 1,\ s + 1.25 - \epsilon)$    | $4.25 - \epsilon$ |
| time-dependent Neumann             | $\frac{3}{4} - \epsilon$ | $\min(2s - 1,\ s + 1.75 - \epsilon)$    | $4.75 - \epsilon$ |
| homogeneous Dirichlet              | $\frac{5}{4} - \epsilon$ | $\min(2s - 1,\ s + 2.25 - \epsilon)$    | 5                 |
| homogeneous Neumann                | $\frac{7}{4} - \epsilon$ | $\min(2s - 1,\ s + 2.75 - \epsilon)$    | 5                 |

*Table 4.3. Temporal order of Radau IIA on the heat equation in the discrete $L^2$ norm in space,
from the exponents* $\beta$ *of [1993LO, §3] `[cited]`; the maximum norm loses* $\frac{1}{4}$. *With
smooth, compatible data the classical order* $2s - 1$ *is recovered.*

### 4.3 Stage equations

With step size $h$ and stage increments $Z_i = Y_i - y_n$, one step solves the $sN$ equations

$$
M Z_i = h \sum_{j=1}^{s} a_{ij} f(t_n + c_j h, y_n + Z_j),
\quad i = 1, \dots, s,
$$

and then sets

$$
y_{n+1} = y_n + \sum_{j=1}^{s} d_j Z_j,
\quad d^T = b^T A^{-1}.
$$

#### Where does $d$ come from?

Write $F_j = f(t_n + c_j h, Y_j)$. The stage equations read
$(I_s \otimes M) Z = h (A \otimes I_N) F$, and $A$ is nonsingular for a collocation method
$\left(\det A = \frac{\prod_i c_i}{s!} \right)$. Hence $h F = (A^{-1} \otimes M) Z$, and the quadrature update
$M (y_{n+1} - y_n) = h \sum_j b_j F_j$ becomes $M \sum_j d_j Z_j$. The update needs neither the
$F_j$ nor $M^{-1}$.

#### Why is $d = \mathbf{e}_s$ for Radau IIA?

Let $\mathbf e_s = (0, \dots, 0, 1)^T$ be the last unit
vector of $\mathbb{R}^s$. Stiff accuracy means $b^T$ is the last row of $A$, that is
$b^T = \mathbf e_s^T A$, so $d^T = \mathbf e_s^T$ and $y_{n+1} = y_n + Z_s$ `[proved]`.
Equivalently, $c_s = 1$ puts the last collocation point at the end of the step, so the last stage
value is the collocation polynomial there. The step ends with one vector addition.

#### Why use increments?

The Newton iteration determines $Z$ to the tolerance $\epsilon_{\text{nls}}$. Any error
in $Y_j$ is amplified by $\lVert J \rVert$ in $F_j$, and $h \lVert J \rVert \gg 1$ on stiff
components, so the form $y_n + h \sum_j b_j F_j$ inherits that amplification while
$y_n + \sum_j d_j Z_j$ does not [1996HW, §IV.8]. The same holds for roundoff in $F_j$.

#### General tables

The integrator stores $d$ with every table and applies the sum above in general
form, so a table that is not stiffly accurate, a Gauss or a singly implicit one, runs on the same
core (Milestone M10). For Radau IIA the sum collapses to a copy of $Z_s$.

#### Notation notes

The bold $\mathbf{e}_s$ is the unit vector. The italic $e$ of Section 4.6 is the
estimate weight vector $e^T = (\hat b - b)^T A^{-1}$; its last entry is $e_s = (-1)^s \cdot \frac{\gamma_0}{s}$
`[proved]`.

### 4.4 Simplified Newton iteration

Stack the increments in $Z \in \mathbb{R}^{sN}$ and let $G(Z)$ be the residual of the stage
equations, $G_i(Z) = M Z_i - h \sum_j a_{ij} F_j$. Starting from the predictor $Z^{(0)}$ of Section
4.9, iteration $k = 1, 2, \dots, m$ computes

$$
(I_s \otimes M - hA \otimes J) \delta^{(k)} = -G(Z^{(k-1)}),
\quad Z^{(k)} = Z^{(k-1)} + \delta^{(k)},
$$

with $J$ evaluated at $(t_n, y_n)$ and reused over several steps. The convergence test is RADAU5's
[1996HW, §IV.8], in the WRMS norm over all $sN$ components, with every stage weighted by the weights of
$y_n$. Write $\hat\theta_k = \frac{\lVert \delta^{(k)} \rVert}{\lVert \delta^{(k-1)} \rVert}$ for the
ratio of successive increments.

#### Rate

$\theta_2 = \hat\theta_2$; for $3 \le k < m$, $\theta_k = \sqrt{\hat\theta_k \hat\theta_{k-1}}$,
the geometric mean of the last two ratios. The convergence factor is
$\eta_k = \frac{\theta_k}{1 - \theta_k}$.

#### First iteration

No rate exists, so $\eta_1 = \max(\eta_{\mathrm{prev}}, u)^{0.8}$, where
$\eta_{\mathrm{prev}}$ is the final factor of the previous Newton solve, 1 before the first.

#### Divergence

If $\theta_k \ge 0.99$ the solve fails and the step is halved.

#### Predicted failure

If $\eta_k \theta_k^{m-1-k} \lVert \delta^{(k)} \rVert > \epsilon_{\text{nls}}$,
convergence is not expected within the remaining iterations; the solve fails and the step is
multiplied by $0.8 q^{-\frac{1}{4 + m - 1 - k}}$ with
$q = \mathrm{clamp}\Bigl(\frac{\eta_k \theta_k^{m-1-k} \lVert \delta^{(k)} \rVert}{\epsilon_{\text{nls}}},\ 10^{-4},\ 20\Bigr)$.
The exponent $m - 1 - k$ is RADAU5's; it predicts at the second-to-last iteration.

#### Stop

The update is applied, and the iteration stops when
$\eta_k \lVert \delta^{(k)} \rVert \le \epsilon_{\text{nls}}$. By default $\epsilon_{\text{nls}}$ is
RADAU5's `FNEWT` [1996HW, §IV.8], evaluated with the relative tolerance as given,

$$
\epsilon_{\text{nls}} = \max\Bigl(\frac{10u}{\mathrm{rtol}}, \min(0.03, \sqrt{\mathrm{rtol}})\Bigr),
$$

recomputed whenever the tolerances change; before the tolerances are set, for $\mathrm{rtol} = 0$
and with a user weight function, which give no relative accuracy to derive it from, it is the cap
0.03. The EXPERT setter fixes a positive value; a non-positive value restores the rule. With the
fixed $\epsilon_{\text{nls}} = 0.1$ of CVODE the final-time errors of the M2 comparison were 6–330
times RADAU5's at equal step counts
`[fact: test/unit_tests/firkode/CXX_serial/firk_test_radau5_steps.cpp]`.

#### Exhaustion

If $k = m$ without convergence the solve fails and the step is halved. Default $m = 7$.

#### Order of the tests

The divergence and prediction tests precede the update of $Z$, and the stop test follows it, as in
RADAU5. The linear systems are solved inexactly, within the theory of [1982DES, 2000J]; Section 5.10
gives the inner tolerances, which are fixed forcing terms in the sense of [1996EW], an adaptive
forcing term being an EXPERT option of Milestone M6 (Open question 10). With difference-quotient
products $Jv$ the attainable relative residual of the stacked solve is of the order of $\sqrt u$
[1989BH, 2004KK], so EXPERT values of $\epsilon_{stk} \epsilon_{\text{nls}}$ below that bound stall
rather than converge `[cited]`. A linear problem converges in at most two iterations when $J$ is
current, which the Dahlquist test of Chapter 11 checks.

### 4.5 Single-shift block preconditioning

The Newton matrix $K = I_s \otimes M - hA \otimes J$ couples all stages. RADAU5 transforms $A^{-1}$ to
real block-diagonal form and factors one real and $\left\lfloor \frac{s}{2} \right\rfloor$ complex $N \times N$ matrices
[1976B, 1977B]. FIRKODE instead solves $K \delta = r$ with right-preconditioned flexible GMRES
[1993S] and the preconditioner

$$
P = I_s \otimes (M - \gamma J),
\quad \gamma = h \gamma_0.
$$

- Odd $s$: $\gamma_0 = \frac{1}{U_1}$, with $U_1$ the unique real eigenvalue of $A^{-1}$ `[proved]`. One
  mode is preconditioned exactly. For $s = 3$ this is RADAU5's real factor.
- Even $s$: $A^{-1}$ has no real eigenvalue `[proved]`. $\gamma_0 = \frac{1}{\mathrm{Re} \lambda_\ast}$, with $\lambda_\ast$ the eigenvalue of $A^{-1}$ of smallest modulus. This is a
  heuristic without an exactly preconditioned mode; even $s$ is an EXPERT-mode option (D7).
- Applying $P^{-1}$ costs $s$ independent solves with one matrix. For $s = 1$, $P = K$ and the
  Krylov iteration is skipped.

Assume $M^{-1} J = V \Lambda V^{-1}$ with $\mathrm{Re} \lambda \le 0$, and let $\mu_j$ be the
eigenvalues of $A$. Then $K P^{-1}$ is similar to a diagonal matrix with the entries

$$
w_j(z) = \frac{1 - z \mu_j}{1 - \gamma_0 z},
\quad z = h \lambda.
$$

Each $w_j$ is a Möbius map with its pole in the right half-plane, so it maps $\mathrm{Re} z \le 0$
onto a closed disk $D_j$ through $1$ and $\frac{\mu_j}{\gamma_0}$. The spectrum lies in $\cup_j D_j$ for
every $h$ and $N$. With $p(w) = \left(1 - \frac{w}{c}\right)^k$ in the GMRES bound [1986SS], exact block solves give

$$
\frac{\lVert r_k \rVert}{\lVert r_0 \rVert} \le \kappa \rho_\ast(\gamma_0)^k,
\quad
\rho_\ast(\gamma) = \min_{c > 0} \max_{w \in \cup_j D_j} \left\lvert 1 - \frac{w}{c} \right\rvert,
$$

where $\kappa$ collects the conditioning of the eigenvector matrices. The bound is independent of the
mesh and the step, the property that makes block preconditioners order-optimal for parabolic problems
[2006SMN, 2007MNS, 2024GO]. For a non-normal $J$ the dependence enters through $\kappa$.

| $s$                          | 2       | 3       | 4       | 5       | 6       | 7       | 8       | 9       |
| ---------------------------- | ------- | ------- | ------- | ------- | ------- | ------- | ------- | ------- |
| $\rho_\ast(\gamma_0)$        | 0.577   | 0.766   | 0.847   | 0.897   | 0.925   | 0.946   | 0.959   | 0.969   |
| $k$ for a 20-fold reduction  | 6       | 12      | 19      | 28      | 39      | 54      | 71      | 96      |

*Table 4.1. Worst-case contraction per FGMRES iteration and the iteration count* $k$ *that
guarantees a 20-fold residual reduction when* $\kappa = 1$ *`[derived]`. Appendix B describes the
computation; the table generator of Milestone M0 reproduces it with `--table1`.*

The bound is pessimistic: GMRES exploits the clustering of the spectrum on $s$ curves, and the
block solves are inexact in practice. Actual iteration counts for every $s$ are `[to measure]`
(Milestone M4, Section 11.4); until then the Krylov defaults of Section 5.10 are provisional. With
$\mathrm{maxl} = 3s$ and one restart, FGMRES allows $6s$ iterations, which falls below the
worst-case $k$ of Table 4.1 for $s \ge 6$ ($36 < 39$).

### 4.6 Temporal error estimate

Let $\hat b$ be the embedded quadrature of order $s$ that adds the left end point with weight
$\gamma_0$, and let $e^T = (\hat b - b)^T A^{-1}$, so that $e_j = -\frac{\gamma_0 l_j(0)}{c_j}$ with $l_j$
the Lagrange basis on the nodes `[proved]`. The estimate is

$$
\mathrm{err} = (M - \gamma J)^{-1} \Bigl[ \gamma f\left(t_n, y_n\right) + \sum_{i=1}^{s} e_i M Z_i \Bigr],
\quad \gamma = h \gamma_0,
$$

of order $p = s$, that is with local error $O(h^{s+1})$ [1996HW, 1999HW] `[proved]` for the
polynomial-exactness identity; the $O(h^{s+1})$ bound for smooth $f$ is the standard argument. For
$s = 3$ the weights are RADAU5's `[proved]`. The construction of embedded estimators for implicit
Runge–Kutta methods in general is analysed in [1997SS].

- The filter $(M - \gamma J)^{-1}$ reuses the preconditioner's matrix and keeps the estimate bounded
  for stiff components as $h \to \infty$ (Figure 4).
- $\mathrm{dsm} = \max\left(\lVert \mathrm{err} \rVert, 10^{-10}\right)$; the floor is RADAU5's and also keeps
  the controller away from $\mathrm{dsm} = 0$. The step is accepted if $\mathrm{dsm} \le 1$.
- On the first step (after `Init`, `ReInit` or `Reset`) and after a rejected attempt of the same
  step, whether by the error test or by the Newton iteration (RADAU5's `REJECT` flag covers both),
  a failing estimate is filtered again with $f(t_n, y_n + \mathrm{err})$ in place of $f(t_n, y_n)$,
  as in RADAU5. The refilter is on by default; turning it off is an EXPERT option.

### 4.7 Step-size control

The controllers of [2003S, 2025KC] have the form

$$
h_{n+1} = \kappa h_n
\Bigl( \frac{\epsilon}{\lVert \delta_n \rVert} \Bigr)^{\alpha}
\Bigl( \frac{\lVert \delta_{n-1} \rVert}{\epsilon} \Bigr)^{\beta}
\Bigl( \frac{\epsilon}{\lVert \delta_{n-2} \rVert} \Bigr)^{\gamma}
\Bigl( \frac{h_n}{h_{n-1}} \Bigr)^{a}
\Bigl( \frac{h_{n-1}}{h_{n-2}} \Bigr)^{b},
$$

with $\kappa$ the safety factor, $\delta_n$ the local error estimate of step $n$, and $k = p + 1$ in
the coefficients for error per step [2025KC, Table 2]:

| Controller | $\alpha$        | $\beta$          | $\gamma$         | $a$            | $b$            | In SUNDIALS                                                 |
| ---------- | --------------- | ---------------- | ---------------- | -------------- | -------------- | ----------------------------------------------------------- |
| I          | $\frac{1}{k}$   | 0                | 0                | 0              | 0              | `SUNAdaptController_I`                                      |
| H211       | $\frac{1}{4k}$  | $-\frac{1}{4k}$  | 0                | $-\frac{1}{4}$ | 0              | `SUNAdaptController_H211`                                   |
| PC         | $\frac{2}{k}$   | $\frac{1}{k}$    | 0                | 1              | 0              | `SUNAdaptController_ImpGus` with $k_1 = k_2 = 1$ (RADAU5's) |
| PID        | $\frac{1}{18k}$ | $-\frac{1}{9k}$  | $\frac{1}{18k}$  | 0              | 0              | `SUNAdaptController_PID`                                    |
| H312       | $\frac{1}{8k}$  | $-\frac{1}{4k}$  | $\frac{1}{8k}$   | $-\frac{3}{8}$ | $-\frac{1}{8}$ | `SUNAdaptController_H312`                                   |
| PPID       | $\frac{6}{20k}$ | $-\frac{1}{20k}$ | $-\frac{5}{20k}$ | 1              | 0              | `SUNAdaptController_Soderlind` with `SetParams_Soderlind`   |
| H321       | $\frac{1}{3k}$  | $-\frac{1}{18k}$ | $-\frac{5}{18k}$ | $\frac{5}{6}$  | $\frac{1}{6}$  | `SUNAdaptController_Soderlind` with `SetParams_Soderlind`   |

*Table 4.2. Step-size controllers and their SUNDIALS realizations. The generic
`SUNAdaptController_Soderlind` evaluates*
$h_{n+1} = h_n \varepsilon_n^{-\frac{k_1}{k}} \varepsilon_{n-1}^{-\frac{k_2}{k}} \varepsilon_{n-2}^{-\frac{k_3}{k}} \left(\frac{h_n}{h_{n-1}}\right)^{k_4} \left(\frac{h_{n-1}}{h_{n-2}}\right)^{k_5}$
*with* $\varepsilon = \mathrm{bias} \cdot \mathrm{dsm}$ *and* $k = p + 1$
*`[fact: src/sunadaptcontroller/soderlind/sunadaptcontroller_soderlind.c]`; the mapping is*
$(k_1, \dots, k_5) = (\alpha k, -\beta k, \gamma k, a, b)$. *There are no built-in PPID or H321
constructors.*

FIRKODE calls `SUNAdaptController_EstimateStep(C, h, p, dsm, &hnew)` with $p = s$, the order of the
estimate, at most once per step attempt (not for a rejected attempt in BEGINNER mode) and before the
history update, and applies RADAU5's heuristics to the proposal afterwards, because the controller
interface has no hook for them:

- the safety factors: RADAU5 multiplies the predictive proposal by `SAFE` $= 0.9$ and the
  I-controller proposal by `FAC` $= \kappa = 0.9 \min\Bigl(1, \frac{1 + 2m}{k_{\text{nls}} + 2m}\Bigr)$,
  with $k_{\text{nls}}$ the Newton iterations of the step, so that a step that needed many iterations
  grows less; FIRKODE does the same in BEGINNER mode, and in EXPERT mode, where no I-controller
  proposal exists, applies $\kappa$ to the controller's proposal (`SetNewtonCountSafety` turns the
  Newton factor off);
- in BEGINNER mode the smaller of the controller's proposal and the I-controller's
  $\kappa h_n \mathrm{dsm}^{-\frac{1}{s+1}}$, which is RADAU5's `QUOT = MAX(QUOT, FACGUS)`; on a
  rejected step BEGINNER mode uses the I-controller proposal alone, because RADAU5 computes `FACGUS`
  on accepted steps only, and EXPERT mode the smaller of $\kappa$ times the controller's proposal and
  the I-controller's, so that a rejected step always shrinks, whatever the controller's history;
- the ratio limits $\frac{h_{n+1}}{h_n} \in [0.2, 8]$; $10^4$ on the first step (after `Init`,
  `ReInit` or `Reset`); at most 1 on the step accepted after a rejection, RADAU5's
  `HNEW = MIN(HNEW, H)` after `REJECT`; $0.1$ when the first step is rejected; $0.3$ after two
  consecutive error-test failures (an ARKODE safeguard, not in RADAU5); $0.5$ after a Newton
  divergence; `hmin` and `hmax` as set by the user, applied last, with RADAU5's minimum step
  $10 u |t_n|$ as a floor of `hmin` for adaptive steps (R18): a failure at that size ends the
  integration;
- the dead band: an accepted step with $\frac{h_{n+1}}{h_n} \in [1, 1.2]$ and $\theta \le 10^{-3}$ keeps
  $h$, $\gamma$ and the factorization, and keeps the step sequence smooth, the property [2006SW] call
  computational stability.

The default controller is `SUNAdaptController_ImpGus` with $k_1 = k_2 = 1$, Gustafsson's predictive
controller as RADAU5 uses it [1994G], which falls back to the I-controller until one accepted step of
history exists (the controller's history size is 1
`[fact: src/sunadaptcontroller/soderlind/sunadaptcontroller_soderlind.c]`), as RADAU5 uses `FACGUS`
from the second accepted step on. The controller's history is updated only on accepted steps and only
after the proposal has been computed, as the `SUNAdaptController` contract and ARKODE's order of calls
require `[fact: src/arkode/arkode_adapt.c, src/arkode/arkode.c]`, and it is reset at an order change
because estimates of different order are not comparable. [2025KC] recommends H321 or PPID; EXPERT
mode may attach any controller.

### 4.8 Initial step

Without a user-supplied `h0`, FIRKODE uses CVODE's heuristic `[fact: src/cvode/cvode.c, cvHin]`: a
lower bound of $100 u \max(\lvert t_0 \rvert, \lvert t_{out} \rvert)$, an upper bound from
$0.1 \lvert t_{out} - t_0 \rvert$ and from $\lVert \dot y \rVert$ against $0.1 \lvert y \rvert + \frac{1}{w}$, a
geometric-mean trial, a finite-difference estimate of $\ddot y$, at most four refinements, and half the
resulting value. With a mass matrix, $\dot y = M^{-1} f$ needs one mass solve at $(t_0, y_0)$ and one
per trial estimate of $\ddot y$, at most four more, as `arkHin` does
`[fact: src/arkode/arkode.c, arkHin and arkYddNorm]`. These solves only feed heuristics, so their
tolerance is relative, $\epsilon_{L,M} \lVert b \rVert$, and a solve that reduced the residual
without reaching it (`SUNLS_RES_REDUCED`) is accepted: an absolute tolerance would compare a residual
in the units of $f$ with weights in the units of $y$, which nonsmooth initial data can make
unattainable. An iterative mass solver is applied once to the difference quotient of $f$ for
$\ddot y$, because the errors of two separate inexact solves for $\dot y$ need not cancel and their
difference is divided by the trial step; a direct solver solves for $\dot y$ at the trial point as
`arkYddNorm` does. Besides the initial step and rootfinding, $M^{-1}$ is applied only in the units
conversion of the inner solves (Section 5.10). RADAU5 instead starts from a fixed default of
$10^{-6}$; this is a deliberate deviation (Table 6.2).

### 4.9 Dense output and predictor

Within a step, the collocation polynomial

$$
u(t_n + \vartheta h) = y_n + \sum_{j=1}^{s} L_j(\vartheta) Z_j
= y_{n+1} + \sum_{j=1}^{s} \bigl(L_j(\vartheta) - \delta_{js}\bigr) Z_j
$$

is available at no extra cost. $L_j$ is the Lagrange basis on $0, c_1, \dots, c_s$ with $L_j(0) = 0$,
and the second form uses $y_{n+1} = y_n + Z_s$, so that the history of a step consists of $y_{n+1}$
and the $s$ increments and no copy of $y_n$ is kept. Derivatives are
$u^{(k)}(t_n + \vartheta h) = h^{-k} \sum_j L_j^{(k)}(\vartheta) Z_j$ for $0 \le k \le s$. The
interpolation error at interior points is $O(h^{s+1})$ [1993HNW, §II.7], and each derivative loses one
order. The polynomial serves

- `FIRKodeGetDky`, derivatives up to order $s$ on the last step `[proved]` for the interpolation
  identities;
- rootfinding (Section 4.10);
- the Newton predictor, by extrapolation of the previous step's polynomial to the nodes of the new
  step, $Z_i^{(0)} = u_{\mathrm{prev}}(t_n + c_i h) - y_n$; the trivial predictor $Z^{(0)} = 0$ is an
  option;
- from v0.3, the predictor across an order change, by evaluation at the new nodes. RADAU restarts
  from zero.

### 4.10 Rootfinding

FIRKODE locates sign changes of user functions $g_i(t, y)$ on the collocation polynomial with ARKODE's
algorithm `[fact: src/arkode/arkode_root.c]`: the functions are evaluated at the end of each step, a
sign change is bracketed, and the root is located by the Illinois method, a modified secant iteration,
to a tolerance of $100 u (\lvert t \rvert + \lvert h \rvert)$. Root directions, inactive-root
warnings and a root at the initial time follow CVODE's and ARKODE's semantics. Section 8.10 lists the
entry points.

### 4.11 In situ order selection (v0.3)

From v0.3 the order can be chosen at run time from what the current order costs, after RADAU
[1999HW, §6]. The candidates are odd, $s \in \{1, 3, 5, 7, 9\}$, within a range
$s_{\min} \le s \le s_{\max}$ (default $3 \le s \le 7$, as RADAU). Even $s$ have no exactly
preconditioned mode and are excluded.

Signals, recorded per accepted step:

- the Newton rate $\theta$ and its filtered value $\tilde\theta \leftarrow \min\Bigl(10, \max\bigl(\theta, \frac{\tilde\theta}{2}\bigr)\Bigr)$;
- the Newton count $k_{\text{nls}}$ and the FGMRES count per Newton iteration $k_{lin}$;
- the step ratio $\eta = \frac{h_{n+1}}{h_n}$ proposed by the controller;
- two flags: a Newton failure that reduced the step, and an unexpected rejection (singular
  $M - \gamma J$ or divergence).

Rule, evaluated whenever $M - \gamma J$ is refactored and at the latest after 10 kept steps:

- $s \to s + 2$ if $k_{\text{nls}} > 1$, $\tilde\theta \le 0.002$, $0.8 < \eta < 1.2$, and (FIRKODE only)
  $k_{lin} \le 2s$: the block preconditioner degrades with $s$ (Table 4.1), so slow stage solves veto
  a raise. No raise within 10 steps of the last change.
- $s \to s - 2$ if $\tilde\theta \ge 0.8$, or either flag is set.

The work per unit step,

$$
W(s) = \frac{k_{\text{nls}}\bigl[ s C_f + k_{lin} s (C_{sol} + C_{Jv}) \bigr]}{h},
$$

with $C_f$, $C_{sol}$, $C_{Jv}$ the running average wall-clock costs of one RHS evaluation, block solve
and Jacobian product, is kept per order as an exponentially weighted average. In v0.3 it is logged and
breaks ties; the reuse-policy hook of v0.4 sees it. The thresholds are RADAU's; the $k_{lin}$ veto is
uncalibrated until the measurements of Milestone M4 exist `[to measure]`.

At a change: the table and $\gamma_0$ switch, the stage stack and the Krylov basis are reallocated,
the predictor evaluates the previous collocation polynomial at the new nodes, the controller is reset
(`SUNAdaptController_Reset`) and $h$ is kept.

## 5. Design decisions and alternatives considered

Each section states the context, the options, the verdict, the decisive reasons, the cost accepted,
the principle followed and the condition under which the decision is revisited. Each decision
carries an identifier by which other chapters cite it: D1–D9 have one section each (5.1–5.9);
D10–D14 are collected in Table 5.8.

### 5.1 D1: a standalone package, not an ARKODE stepper

#### Context

Every ARKODE method added since 2019 (MRIStep, SPRKStep, LSRKStep, SplittingStep,
ForcingStep) is an in-tree stepper behind ARKODE's private `step_*` function-pointer table
`[fact: src/arkode/arkode_impl.h:412–477]`, sharing the evolve loop, rootfinding, interpolation,
adaptivity, statistics, options and bindings. A SUNDIALS reviewer will ask why FIRKODE is not a
`FIRKStep`.

| ARKODE component              | ARKODE today `[fact]`                                                                                                                         | What Radau IIA needs                                                                                            | Core change a FIRKStep would require                                |
| ----------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------- |
| Butcher table check           | `arkStep_CheckButcherTables` rejects implicit tables with entries above the diagonal (`src/arkode/arkode_arkstep.c:2682`, message at `:2804`) | a dense $A$                                                                                                     | a new stepper, or a new table class                                 |
| Linear-solver interface       | ARKLS forms one shifted matrix per stage with $\gamma = h a_{ii}$ (`arkode_arkstep.c:3122`); the factorization is reused across stages        | one stage-coupled $sN \times sN$ system per Newton iteration, solved iteratively with $N \times N$ block solves | a second linear-solver interface inside ARKODE, or an ARKLS rewrite |
| Nonlinear solve               | per stage, through `SUNNonlinearSolver`, with the CVODE-style rate test                                                                       | one stacked iteration with RADAU5's three-way test and step-size recommendation                                 | bypassing `SUNNonlinearSolver` inside a stepper                     |
| Adaptivity                    | the controller receives the embedding order; safety and growth clamps are applied by ARKODE                                                   | $p = s$; RADAU5's Newton-count safety, I-controller minimum and $\theta$-dependent dead band                    | changes to `arkode_adapt.c` or duplication inside the stepper       |
| Interpolation                 | Hermite and Lagrange modules of degree $\le 5$                                                                                                | the collocation polynomial of degree $s$                                                                        | a third interpolation module                                        |
| Mass matrix                   | ARKStep only (`step_supports_massmatrix`, `src/arkode/arkode.c:1598`)                                                                         | yes                                                                                                             | enabling and testing the mass path in a new stepper                 |
| Variable order                | none in ARKODE                                                                                                                                | v0.3                                                                                                            | reallocation of stage storage inside a stepper                      |

*Table 5.1. What a FIRKStep would need from ARKODE.*

| Duplicated by a standalone package                      | Reference implementation `[fact]`                       | Lines      | Mitigation                                                     |
| ------------------------------------------------------- | ------------------------------------------------------- | ---------- | -------------------------------------------------------------- |
| evolve loop, stop tests, failure handling, initial step | `src/cvode/cvode.c`                                     | 5003       | the loop is small once the step attempt is separated           |
| rootfinding                                             | `src/arkode/arkode_root.c`                              | 961        | adapted verbatim; a shared module is proposed for after 1.0    |
| setters, getters, statistics                            | `src/cvode/cvode_io.c`                                  | 1841       | generated documentation tables keep them in step               |
| linear-solver and mass interface                        | `src/cvode/cvode_ls.c`, `src/arkode/arkode_ls.c`        | 1976, 4159 | the stacked solve is new anyway                                |
| command-line options                                    | `src/arkode/arkode_cli.c`                               | 328        | shared `sundials_cli.h` helpers                                |

*Table 5.2. What a standalone package duplicates, with the sizes of the reference implementations.*

#### Verdict

A standalone package `libsundials_firkode` with its own headers, examples, tests and
guide, built from `sundials_core` and the shared object libraries for ManyVector, SPFGMR and the
controllers. No ARKODE or CVODE source is included; callback types are redeclared with identical
signatures so that users pass the same functions.

#### Reasons

1. The stage solver is stage-coupled where ARKLS is per-stage. A FIRKStep would need a second
   linear-solver interface inside ARKODE, and the mass-matrix path that only ARKStep supports.
2. RADAU5's control logic couples the Newton rate, the step controller and the order rule. ARKODE
   splits these into modules with fixed contracts (Table 5.1); honoring the contracts would weaken
   the algorithm, bypassing them would weaken ARKODE.
3. The package touches no ARKODE file and decouples review and release cadence (P1). The real cost is
   the duplication of mature code in Table 5.2, which the milestone plan bounds by adapting rather
   than redesigning.

Composability through `SUNStepper` is not an argument for either side (P3): FIRKODE is not designed as
an inner stepper, and the glue is the last optional milestone (M8).

#### Cost

About 10 000 lines that track fixes in CVODE and ARKODE; a separate documentation tree; separate
binding work.

Revisit when ARKODE gains a stage-coupled linear-solver interface or a collocation interpolant, and in
any case before 1.0 (Chapter 12).

### 5.2 D2: one real shift and flexible GMRES, not the complex transformation

#### Context

RADAU5's transformation yields exact Newton solves with one real and $\lfloor \frac{s}{2} \rfloor$
complex factorizations per setup. SUNDIALS has real vectors, matrices and solvers only.

| Option                                                                                                                                           | Newton solve exact        | Factorizations per setup                            | Solves per Newton iteration                                  | Works with                                                            | New SUNDIALS types                                    | Status                                         |
| ------------------------------------------------------------------------------------------------------------------------------------------------ | ------------------------- | --------------------------------------------------- | ------------------------------------------------------------ | --------------------------------------------------------------------- | ----------------------------------------------------- | ---------------------------------------------- |
| (a) RADAU5 transformation, complex LU [1976B, 1977B]                                                                                             | yes                       | 1 real + $\lfloor \frac{s}{2} \rfloor$ complex      | $1 + \lfloor \frac{s}{2} \rfloor$ (complex counted once)     | dense and banded direct solvers                                       | complex `sunrealtype`, `SUNMatrix`, `SUNLinearSolver` | rejected                                       |
| (b) transformation with real $2N \times 2N$ blocks                                                                                               | yes                       | 1 real + $\lfloor \frac{s}{2} \rfloor$ of size $2N$ | as (a)                                                       | direct solvers; defeats $N \times N$ user preconditioners             | a $2N$ block interface                                | rejected                                       |
| (c) FGMRES with $I_s \otimes (M - h\gamma_0 J)$                                                                                                  | no, to tolerance          | 1 real                                              | $s k_{lin}$ block solves and $J$-products                    | every `SUNLinearSolver`, preconditioner-only Krylov, AMG, batched GPU | none                                                  | chosen                                         |
| (d) $s$ real factorizations with a real-spectrum or diagonal $B \approx A$ [1991HS, 1997HS, 1999JB]; the SDC sweep matrices of [2018S, 2025CLGR] | no, inner iteration       | $s$ real                                            | $s$ independent                                              | direct solvers, stage-parallel                                        | none                                                  | possible EXPERT stage preconditioner (M6, M10) |
| (e) one real preconditioner per conjugate pair [2022SKPD, 2022SKP]                                                                               | no, to tolerance          | $\lceil \frac{s}{2} \rceil$ real                    | $\lceil \frac{s}{2} \rceil$ real solves per Krylov iteration | preconditioner-only solvers                                           | none                                                  | possible EXPERT stage preconditioner (M6)      |
| (f) user-supplied monolithic solver, e.g. multigrid [2005VV, 2024K]                                                                              | user's                    | user's                                              | user's                                                       | PDE codes with a stage-system hierarchy                               | none                                                  | EXPERT hook `SetStageLinearSolver` (M6)        |
| (g) stage-parallel solve with an exact low-rank Sylvester correction [2026DM]                                                                    | yes, after the correction | $s$ real                                            | $s$ independent, plus one $s \times s$ matrix equation       | direct solvers, stage-parallel                                        | none                                                  | possible EXPERT stage solver (M6)              |

*Table 5.3. Options for the Newton systems.*

#### Verdict

Option (c), with (d), (e), (f) and (g) reachable through the v0.4 hooks.

#### Reasons

1. SUNDIALS has no stable complex types, and a real $2N \times 2N$ embedding cannot be served by a user's
   $N \times N$ solver or preconditioner. Option (c) needs nothing new.
2. The user's CVODE and ARKODE linear-solver interface is reused unchanged: one block $M - \gamma J$,
   the same callbacks, the same reuse controls. Algebraic multigrid, batched GPU solvers and
   matrix-free Krylov solvers with user preconditioners all work (P4: the stage solve is where the
   time goes, and the user's best solver must be usable there).
3. One real factorization per setup for any $s$, and a residual bound independent of mesh and step
   (Section 4.5).

#### Cost

$s k_{lin}$ block solves and $J$-products per Newton iteration instead of
$1 + \lfloor \frac{s}{2} \rfloor$ exact solves; for small dense systems with infrequent setups the complex
factorization is cheaper in flops (Table 7.2). FGMRES stores $2\mathrm{maxl} + 4$ stacked vectors,
which dominate the memory footprint (Table 8.7).

#### Details

Flexible GMRES is required because iterative block solves are inexact and vary between
iterations [1993S]. With a matrix-free user solver, only the user's preconditioner is applied block
by block; the user's Krylov method serves the $s = 1$ Newton solve and the error estimate. The
products $J v$ use the saved matrix, the user's $Jv$ routine, or difference quotients [1989BH]; the
matrix-free path is a Jacobian-free Newton–Krylov method in the sense of [2004KK], with the user's
preconditioner for $M - \gamma J$ in the role of the physics-based preconditioner. A spectral deferred
correction sweep with a diagonal matrix $Q_\Delta$ is a Richardson iteration on the stage system
preconditioned by $I_s \otimes M - h Q_\Delta \otimes J$ [2018S]; the single shift is the case
$Q_\Delta = \gamma_0 I_s$, and the MIN-SR matrices of [2025CLGR], optimized on the Radau nodes for the
non-stiff and for the stiff limit, are candidate data for `SetStagePreconditioner` (M6), to be
measured against Table 4.1 `[to measure]`.

Revisit when the measured $k_{lin}$ for $s \ge 5$ exceeds about $2s$ on the PDE tests of Section 11.4,
or when small dense problems dominate the user base.

### 5.3 D3: an internal Newton loop, not a `SUNNonlinearSolver`

#### Context

CVODE, ARKStep, MRIStep and IDA drive their Newton iterations through
`SUNNonlinearSolver`, and `SUNNonlinSolSetConvTestFn` and `SUNNonlinSolSetNormFn` exist so that an
integrator can supply its own test and norm `[fact: include/sundials/sundials_nonlinearsolver.h]`.
Table 5.4 maps RADAU5's test onto that interface.

| RADAU5 element                                       | Through `SUNNonlinSol_Newton`                                                                                                                                     | Problem                                                                                                                              |
| ---------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------ |
| $\theta_k$, geometric mean, $\eta_k$                 | state in the test's data (the package memory)                                                                                                                     | none                                                                                                                                 |
| stop test $\eta_k \lVert \delta \rVert \le \epsilon$ | return `SUN_SUCCESS`                                                                                                                                              | none                                                                                                                                 |
| all-stage WRMS norm                                  | `SetNormFn` on the stacked vector                                                                                                                                 | none                                                                                                                                 |
| divergence and predicted failure                     | return `SUN_NLS_CONV_RECVR`, write the step factor into the package memory                                                                                        | the test is consulted after the update $Z \leftarrow Z + \delta$, RADAU5 tests before it; harmless but a second ordering to document |
| step-size recommendation                             | a side channel, since the return code carries none                                                                                                                | a hidden contract between two modules                                                                                                |
| no retry at the same $h$                             | `LSetup` must always report `jcur = TRUE`, or the module retries once at the same $h$ with a fresh Jacobian (`src/sunnonlinsol/newton/sunnonlinsol_newton.c:361`) | the module's retry policy must be neutralized to get RADAU5's policy of reducing $h$ first                                           |
| iteration count $m = 7$ and the index $k$            | `SetMaxIters(7)`; `GetCurIter` is 0-based and read before the increment                                                                                           | off-by-one hazards in the prediction exponent                                                                                        |

*Table 5.4. RADAU5's convergence test mapped onto `SUNNonlinSol_Newton`.*

#### Verdict

An internal loop in `firkode_nls.c` that owns the three exits and the step-size
recommendation. The hook signatures of v0.4 (`SetNlsConvTestFn`, `SetNlsNormFn`) mirror those of
`SUNNonlinearSolver`, so that a later migration to the shared module remains possible without an API
change. No `SetNonlinearSolver` is offered.

#### Reasons

1. RADAU5's test has three exits, two of which carry a step-size factor to the step driver and one of
   which is evaluated before the update. The `SUNNonlinearSolver` interface returns only
   `SUN_SUCCESS`, `SUN_NLS_CONTINUE` or `SUN_NLS_CONV_RECVR`, so the factor would travel through a
   hidden side channel.
2. The shared Newton module retries once at the same $h$ with a fresh Jacobian when the Jacobian is
   stale. RADAU5 reduces $h$ first and recomputes $J$ only if it was stale. Making the module follow
   RADAU5 requires neutralizing its retry, which leaves a pass-through object (P4: the Newton loop is
   the performance-critical component, and its policy must be explicit).
3. The stacked linear solve needs the iteration index anyway, because `SUNLS_RES_REDUCED` is accepted
   only at the first iteration (Section 8.6).

#### Cost

FIRKODE cannot host a user-supplied nonlinear solver, and its Newton statistics and log labels
are produced by the package rather than by the shared module. The counters and labels follow the
module's names so that `suntools` parses them.

Revisit when the hook interface can express the three exits and the retry policy, for example a
`SUNNonlinearSolver` return code that carries a step factor.

### 5.4 D4: the stage stack is an `NVECTOR_MANYVECTOR`

#### Context

The $s$ increments must live in one vector for FGMRES, and each stage must remain a vector
of the user's type for the RHS, Jacobian and block-solve callbacks.

| Option                                     | New vector code          | SPFGMR reuse         | GPU kernels per operation | Reductions per dot product (MPI subvectors) | Hazards `[fact: src/nvector/manyvector/nvector_manyvector.c]`                                                                                                                                                                                                                                                                                                   |
| ------------------------------------------ | ------------------------ | -------------------- | ------------------------- | ------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| ManyVector of $s$ clones                   | none                     | unchanged            | $s$                       | $s$, one per subvector                      | over MPI subvectors (`NVECTOR_PARALLEL`, `MPIPLUSX`) the plain kind's `N_VDotProd` and `N_VWrmsNorm` sum the subvectors' global operations, one reduction each, while `N_VMaxNorm`, `N_VMin` and `N_VL1Norm` use the subvectors' `*Local` routines and are rank-local; `N_VGetLocalLength` has no op; `N_VLinearCombination` mallocs one pointer array per call |
| MPIManyVector over the user's communicator | none                     | unchanged            | $s$                       | 1                                           | subvector communicators must be congruent; `N_VNew_MPIManyVector` returns NULL silently when no subvector has one; every clone performs a collective `MPI_Comm_dup`                                                                                                                                                                                             |
| contiguous $N \times s$ multivector        | a new module per backend | needs the new vector | 1                         | 1                                           | large new code surface                                                                                                                                                                                                                                                                                                                                          |
| $s$ bare vectors with hand-written loops   | none                     | lost                 | $s$                       | $s$                                         | a private Krylov solver                                                                                                                                                                                                                                                                                                                                         |

*Table 5.5. Options for stage storage.*

#### Verdict

A ManyVector of $s$ clones of the user's vector [2022G+]. When SUNDIALS is built with MPI and
`N_VGetCommunicator(y0)` is not `SUN_COMM_NULL`, FIRKODE builds an `NVECTOR_MPIMANYVECTOR` over that
communicator instead, so that each dot product and norm is one allreduce instead of $s$; the plain
kind computes them correctly too, but its max norm, min and L1 norm are rank-local over distributed
subvectors (Table 5.5) `[fact: src/nvector/manyvector/nvector_manyvector.c:330–373, 754–964]`.
Because `N_VGetCommunicator` returns `SUN_COMM_NULL` for a vector without that operation
`[fact: src/sundials/sundials_nvector.c]`, a distributed `y0` must implement `nvgetcommunicator`;
this is a documented requirement that FIRKODE cannot verify. What it can verify it does: a `y0` with
a communicator must also implement `nvdotprodlocal` and `nvwsqrsumlocal`, or `Init` returns
`FIRK_ILL_INPUT`. The check is a cost guard: without them the MPI kind falls back to the
subvectors' global operations, which stays correct but reduces once per stage and once more. Stage
access goes through one internal accessor that dispatches on the kind chosen in `Init`, not on
`N_VGetVectorID`, which is that of a ManyVector for both kinds and also for a ManyVector `y0`. For
$s = 1$ no stack exists. The stage solver uses only `N_VConst`, `N_VScale`, `N_VLinearSum`,
`N_VProd`, `N_VDiv`, `N_VDotProd`, `N_VWrmsNorm`, and inside SPFGMR `N_VDotProdMulti` and
`N_VLinearCombination`, which are correct on both kinds; the forbidden operations of Table 8.8 and
the `Init` check are enforced by unit tests with mock vectors (H9). Clones happen only at
initialization and at order changes.

#### Reasons

1. No new vector code; SPFGMR runs unchanged on the $sN$ system (P1).
2. Any `N_Vector` implementation, serial, MPI or GPU, serves as the subvector.
3. The known hazards are avoidable by construction.

#### Cost

Each Newton iteration evaluates $f$ $s$ times; every stacked operation launches $s$ kernels on a
GPU; the plain ManyVector's fused operations allocate host memory per call, so the stage solver avoids
them outside SPFGMR and uses modified Gram–Schmidt by default.

Revisit when GPU launch overhead is measured to dominate, in favor of a contiguous multivector.

### 5.5 D5: `SUNAdaptController` with RADAU5's heuristics applied outside it

#### Context

RADAU5's step control is Gustafsson's predictive controller with safety, limits and a dead
band woven into one routine. SUNDIALS separates the controller from the clamps.

| Element                                | RADAU5                                                                                                              | `SUNAdaptController_ImpGus` `[fact: src/sunadaptcontroller/soderlind/sunadaptcontroller_soderlind.c]` | FIRKODE                                                                                                                |
| -------------------------------------- | ------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------- |
| exponents                              | $\frac{1}{4}$ on both errors, $s = 3$                                                                               | $k_1 = 0.98$, $k_2 = 0.95$ by default, with `ord` $= p + 1$                                           | $k_1 = k_2 = 1$ with $p = s$ in BEGINNER mode                                                                          |
| minimum with the I-controller proposal | yes                                                                                                                 | no                                                                                                    | applied by FIRKODE in BEGINNER mode                                                                                    |
| floor on the stored error              | $\max(10^{-2}, \mathrm{dsm})$ for the history                                                                       | none                                                                                                  | current `dsm` floored at $10^{-10}$; the default stores $\max(10^{-2}, \mathrm{dsm})$ on acceptance                             |
| first steps                            | I-controller until history exists                                                                                   | I-controller until history exists                                                                     | same                                                                                                                   |
| safety factor                          | $0.9 \min\Bigl(1, \frac{1 + 2m}{k_{\text{nls}} + 2m}\Bigr)$ on the I-controller, `SAFE` $= 0.9$ on `FACGUS`, inside | none                                                                                                  | $0.9$ on the controller's proposal, $\kappa$ on the I-controller minimum (BEGINNER); $\kappa$ on the proposal (EXPERT) |
| ratio limits, dead band                | inside                                                                                                              | none                                                                                                  | applied to the proposal                                                                                                |
| history update                         | on accepted steps                                                                                                   | `UpdateH` on accepted steps only                                                                      | same                                                                                                                   |
| `dsm = 0`                              | floored                                                                                                             | `pow(0, negative)` is infinite                                                                        | floored                                                                                                                |

*Table 5.6. RADAU5's controller against the SUNDIALS implementation.*

#### Verdict

Delegate the proposal to `SUNAdaptController` with $p = s$ and apply RADAU5's rules
afterwards (Section 4.7). BEGINNER default: `ImpGus` with $k_1 = k_2 = 1$. The default
controller stores $\max(10^{-2}, \mathrm{dsm})$ on accepted steps, as RADAU5 does. This
history floor does not affect acceptance or the current-step proposal. Explicitly supplied
controllers, including controllers selected by name, receive the actual accepted error.

#### Reasons

1. Shared, tested, user-selectable controllers (P1).
2. ARKODE precedent for applying the clamps outside the controller.
3. The differences from RADAU5 are enumerable and recoverable through controller parameters.

#### Cost

The RADAU5 comparison of Milestone M2 must account for the one remaining difference.

Revisit never; only the defaults.

### 5.6 D6: tables as generated literals

#### Context

The coefficients of the $s$-stage method are algebraic numbers that can be computed at run
time or stored as literals.

#### Verdict

A Python generator `scripts/firkode_radau_tables.py` (Milestone M0) computes $c$, $A$,
$A^{-1}$, $b$, $d$, $e$, $\gamma_0$ and the dense-output coefficients $P$ for $1 \le s \le 9$ with
`mpmath` at 50 digits and emits `src/firkode/firkode_tables.c` as `SUN_RCONST` decimal literals with
40 significant digits, as ARKODE's Butcher tables carry `[fact: src/arkode/arkode_butcher_dirk.def]`:
enough for every `long double` in use, including the 113-bit binary128 of aarch64 Linux and POWER,
where 21 digits would fail the extended-precision table checks of M0. Nothing is computed at run
time. The generator
has `--emit-c`, `--check` (fails CI when the emitted file is stale), `--table1` (Table 4.1) and
`--emit-lean` (Milestone M7). The stability-plot script of Appendix A imports the same generator as its
reference.

#### Reasons

1. Run-time construction in `long double` has no guard digits on Apple arm64 and MSVC, where
   `long double` is `double`; literals are identical on every platform.
2. Literals are reviewable and can be certified directly in Lean (Section 11.5), so the generator
   itself needs no trust (P1).
3. `FIRK_MAX_STAGES = 9` is fixed by the data. For $s = 9$ the smallest violated order condition,
   $B(18)$, has residual $9.4 \cdot 10^{-11}$ `[proved]`, about 5 times the double-precision
   tolerance $10^3 s^2 u = 1.8 \cdot 10^{-11}$ of `FIRKodeTable_CheckOrder` ($u$ = `SUN_UNIT_ROUNDOFF`);
   beyond $s = 9$ the order check itself becomes unreliable
   in double precision. In single precision the conditions hold to roundoff for $s \le 4$, but the
   tolerance separates satisfied from violated conditions only for $s \le 3$: for $s = 4$ the
   residual of $B(8)$, $\frac{(4!\,3!)^2}{8!\,7!} = 1.0 \cdot 10^{-4}$, lies below
   $10^3 \cdot 16 u = 1.9 \cdot 10^{-3}$ `[derived]`.

#### Cost

A generated file and a CI check.

Revisit when adding $s > 9$.

### 5.7 D7: order, by version

#### Context

RADAU selects its order at run time and CVODE is variable-order, but the block preconditioner
degrades with $s$ (Table 4.1), so exposing every stage count before its Krylov counts are measured
would violate P2.

#### Verdict

- v0.1.0: $s = 3$ only (P2). `SetNumStages` and `SetOrder` exist in the API but accept only
  $s = 3$ (`SetOrder`: $q \in \{4, 5\}$, or $q \le 0$ for the default) and return
  `FIRK_ILL_INPUT` otherwise, so that the signature is stable. The general-$s$ code path stays; the
  unit tests reach it through an internal setter, `firkSetNumStagesAny`, which becomes the body of
  `SetNumStages` in v0.2.
- v0.2 (M4): any $1 \le s \le 9$ in EXPERT mode. Even $s$ use the heuristic $\gamma_0$ of Section
  4.5 and are EXPERT-only.
- v0.3 (M5): in situ odd order within $3 \le s \le 7$ becomes the BEGINNER default; EXPERT may fix
  $s$ or set another range.

#### Reasons

1. RADAU's precedent for variable order, and CVODE's within SUNDIALS.
2. The cost signals exist anyway ($\theta$, $k_{\text{nls}}$, $k_{lin}$).
3. The degradation of the preconditioner with $s$ makes a veto necessary (P4).

#### Cost

Three versions before the full family is the default; `SetNumStages` and `SetOrder` reject every
input but $s = 3$ until v0.2.

Revisit after the calibration data of M4.

### 5.8 D8: BEGINNER and EXPERT modes

#### Context

A user without a deep understanding of fully implicit methods should see a small, safe API.
No SUNDIALS package gates setters by mode today. Three options were considered: (i) runtime gating
with rejection; (ii) modes that change defaults and warn; (iii) documentation tiers only, as
ARKODE's guide.

#### Verdict

Option (i), runtime gating with rejection:

- `FIRKodeSetMode(mem, FIRK_BEGINNER | FIRK_EXPERT)`, default BEGINNER, allowed while no step has
  been attempted since `Create` or `ReInit` (not after `Reset`, which keeps the history). Returning
  from EXPERT to BEGINNER mode restores the default of every EXPERT option, so that the guarantee
  below also holds after a change of mind.
- BEGINNER accepts the essentials: `Create`, `Init`, `ReInit`, `Reset`, `Free`, the tolerance
  functions, `SetUserData`, `SetLinearSolver` with the Jacobian, $Jv$, preconditioner and linear-system
  callbacks, the mass-matrix interface, `SetStopTime`, `SetInterpolateStopTime`, `ClearStopTime`,
  `SetMaxNumSteps`, `SetInitStep`, `SetMinStep`, `SetMaxStep`, `SetFixedStep`, rootfinding,
  `Evolve`, `GetDky`, every getter and statistic, `WriteParameters`, `PrintAllStats`, `SetOptions`
  for the accepted keys, and `SetMode`. The minimal program is five calls: `Create`, `Init`,
  `SStolerances` or `SVtolerances`, `SetLinearSolver`, `Evolve`.
- BEGINNER also accepts the inputs that do not change the method: `SetMaxHnilWarns`, the tolerance
  and scaling of the user's solvers (`SetEpsLin`, `SetLSNormFactor`, `SetLinearSolutionScaling`,
  `SetMassEpsLin`, `SetMassLSNormFactor`) and `SetJacTimesRhsFn` (Table 9.1).
- Every other setter (Newton, adaptivity, controller, reuse policy and its frequencies, stage
  solver, method and table selection, hooks) returns `FIRK_ILL_INPUT` with a message naming the
  function and the remedy `FIRKodeSetMode(mem, FIRK_EXPERT)`; Table 9.1 lists them.
- `WriteParameters` prints the mode and every effective value.

#### Reasons

1. A rejection is a return code, so bindings and command-line options are unaffected.
2. The guarantee that a BEGINNER run uses the validated defaults is absolute rather than advisory
   (P2).
3. Alternative (ii) gives the information without the guarantee; (iii) has precedent but no
   guarantee.

#### Cost

One mode check per setter, a user-guide column stating the mode of every function, and a precedent
the maintainers must accept; option (ii), the same API with warnings instead of rejections, is the
fallback if they decline.

Revisit on user feedback.

### 5.9 D9: reuse policies

#### Context

RADAU5 recomputes $J$ after almost every step and refactors whenever $h$ changes; CVODE reuses both
for many steps. The right choice depends on the cost of a setup relative to a Newton iteration.

#### Verdict

A reuse policy decides, at the start of each step attempt, whether to re-evaluate $J$
and whether to refactor $M - \gamma J$. It sees $\theta$, $k_{\text{nls}}$, $k_{lin}$, dsm, $\eta$,
$\frac{\Delta\gamma}{\gamma}$, the steps since the last $J$ and setup, and the failure flags. FIRKLS
executes the decision; it does not make it.

- `FIRK_POLICY_RADAU5`, the default in both modes: $J$ is recomputed after every accepted step
  unless $\theta \le \theta_{jac}$ ($10^{-3}$), and after a failed attempt only if the current
  point has none (RADAU5's `CALJAC`); the factorization is redone whenever $J$ or $\gamma$ changed,
  that is, when $h$ leaves the dead band, after any failure, and at an order change; a fixed step
  reuses it.
- `FIRK_POLICY_PERIODIC`, CVODE's: setup on the first step, after failures, every 20 steps and when
  $\gamma$ changed by more than 30%; $J$ on the first step, every 51 steps and after a convergence
  failure with a Jacobian that was not current `[fact: src/cvode/cvode_ls.c, cvLsSetup]`, if the
  failed attempt's $\gamma$ was within `dgmax_jbad` (0.2, `SetDeltaGammaMaxBadJac`) of its setup's:
  CVLS's attribution of the failure to the Jacobian rather than to $\gamma$ (R9). CVLS makes it at
  the shared module's retry at the same $h$, which FIRKODE does not have (D3).
- Under both policies the first step and a Newton rate above `SetJacBadConvRate`'s on the last
  accepted step (off by default; H1) ask for $J$ and a setup at the first attempt of a step.
- After a convergence failure the step is reduced by the Newton driver's factor and $J$ is
  re-evaluated only if it was stale at the failed setup (D14). A Jacobian is never re-evaluated at
  the same $(t_n, y_n)$ (H8): FIRKODE records whether the current point has one.

#### Reasons

1. RADAU5's policy makes step counts comparable with RADAU5.
2. CVODE's policy suits expensive block setups such as algebraic multigrid.
3. Everything is linearized at $(t_n, y_n)$, so a CVODE preconditioner for $I - \gamma J$ is reused
   unchanged.

#### Cost

Two policies to test under every failure path (Table 8.3). EXPERT mode (v0.4) adds a user policy.

Revisit after the M4 measurements with band, PCG and AMG block solvers (Section 11.4): if the
better of the two policies depends on the block solver, a cost-based policy on the step record, the
data the v0.4 hook already receives, becomes a third built-in policy and the BEGINNER default is
reconsidered. A third lever, recycling part of the FGMRES subspace across Newton iterations and steps
[2006P+], since consecutive stage systems differ only through $Z$ and $h$, is Open question 11.

### 5.10 Tolerances of the inner solves

A block solve gets its tolerance from its purpose. $\epsilon_L$ is set by `FIRKodeSetEpsLin`,
default 0.05. Tolerances are stated in the WRMS norm; Table 8.6 gives the conversion to the 2-norm
stopping tests of the SUNDIALS Krylov solvers.

| Purpose                                                    | Tolerance (WRMS)                                                                                                                                                                                                                                                                            | On non-convergence                                            |
| ---------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------- |
| Newton system, $s = 1$                                     | $\epsilon_L \epsilon_{\text{nls}}$, as in CVODE; tested on $x_0$ (see below)                                                                                                                                                                                                                | recoverable failure                                           |
| error estimate                                             | $\epsilon_L \min(0.1, \lVert x_0 \rVert)$: the fixed accuracy 0.1 of the estimate, which decides the error test against 1, tightened for a small estimate so that it is resolved rather than returned as zero; independent of $\epsilon_{\text{nls}}$, which follows the relative tolerance | recoverable failure                                           |
| preconditioner application                                 | direct block solver: exact; iterative block solver: one `psolve` per block with $\epsilon_L \lVert b \rVert$, relative, since $b$ is a Krylov vector of arbitrary scale, and no inner Krylov iteration                                                                                      | accepted; FGMRES absorbs it                                   |
| mass solve ($\dot y$ for the initial step and rootfinding) | $\epsilon_{L,M} \lVert b \rVert$, relative (Section 4.8)                                                                                                                                                                                                                                    | `SUNLS_RES_REDUCED` accepted; otherwise `FIRK_MASSSOLVE_FAIL` |
| mass solve ($x_0$, the units conversion below)             | $\epsilon_{L,M} \lVert b \rVert$, relative to the whole right-hand side $b$ of the enclosing solve                                                                                                                                                                                          | failure of the enclosing solve, with its sign                 |

*Table 5.7. Tolerances of the inner solves.*

The stacked FGMRES stops when the residual norm falls below $\epsilon_{stk} \epsilon_{\text{nls}}$, default
$\epsilon_{stk} = 0.05$. With right preconditioning this is the true residual. Krylov dimension
$\max(5, 3s)$, one restart, modified Gram–Schmidt; these defaults are provisional until Milestone M4
measures iteration counts.

The tolerances bound the WRMS norm of a correction, which has the units of $y$, with the weights of
$y_n$. With a mass matrix the right-hand side $b = -G$ of a Newton or estimate solve has the units of
$M y$, so for entries of $M$ of size $m$ a test on $\lVert b \rVert$ accepts a correction of size
$\frac{\epsilon}{m}$: with P1 finite-element masses, $m \sim \Delta x$, the stacked solve then accepts the
predictor on fine meshes and the step size collapses `[derived]`, which `firk_test_mass_fem` guards
against `[fact: test/unit_tests/firkode/CXX_serial/firk_test_mass_fem.cpp]`. FIRKODE therefore
converts $b$ to the units of $y$ by mass solves, $x_0 = M^{-1} b$, which is the right-hand side of
the same solve for the equivalent problem $\dot y = M^{-1} f$; with $M = I$ nothing changes:

- stacked solve, $s \ge 2$: $x_0 = (I_s \otimes M)^{-1} b$, one mass solve per stage, whatever the
  block solver; the FGMRES tolerance is multiplied by $\frac{\lVert b \rVert}{\lVert x_0 \rVert}$, which
  is consistent because FGMRES tests the true residual (right preconditioning only);
- $s = 1$ Newton solve and estimate with an iterative block solver: $x_0 = M^{-1} b$ replaces $b$ in
  the small-residual test of Listing 8.4 and in the estimate tolerance, and the Krylov residual is
  converted with $\sigma = \frac{\lVert x_0 \rVert}{\lVert b \rVert}$: a matrix-free solver solves the
  block system multiplied by $\sigma$ (right-hand side, products and preconditioner), which leaves
  its iterates and a left-preconditioned residual unchanged and puts a right-preconditioned or
  unpreconditioned residual in the units of $y$; a matrix-based iterative solver, whose test FIRKODE
  cannot see, gets the tolerance divided by $\sigma$;
- direct block solvers solve the $s = 1$ Newton system and the estimate exactly and need no
  conversion.

The cost is $s$ mass solves per Newton iteration for $s \ge 2$ and one per iterative block solve for
the $s = 1$ Newton solve and the estimate (Table 7.1); a conversion that fails makes the enclosing
solve fail, recoverably or not, with the sign of the mass-solve failure.

The tolerances above are fixed forcing terms. An adaptive forcing term in the manner of [1996EW],
which tightens $\epsilon_{stk}$ as the Newton iteration contracts and relaxes it when the previous
linear solve was more accurate than the nonlinear progress required, is compatible with the inexact
simplified Newton theory of [2000J] and is standard practice in Newton–Krylov solvers [2004KK]. It is
an EXPERT option of Milestone M6, bounded above by the fixed value and below by the $\sqrt u$ floor of
difference-quotient products (Section 4.4); the default stays fixed until the FGMRES counts of Section
11.4 say otherwise `[to measure]` (Open question 10).

### 5.11 Conventions

- Errors: `int` return codes `FIRK_*` and `FIRKLS_*` with the explicit numbering of Table 8.2 (D10).
  Messages go through `firkProcessError`, as in the other packages; warnings through the logger.
  `SUNCheckCall` and `SUNAssert` are used only in the optional `SUNStepper` glue, as ARKODE does
  `[fact: src/arkode/arkode_sunstepper.c]`.
- Logging: `SUNLogInfo` and `SUNLogDebug` records per step, Newton iteration and Krylov iteration
  with ARKODE's labels (Table 8.4), from the first milestone on.
- Profiling: `SUNDIALS_MARK_FUNCTION_BEGIN` and `SUNDIALS_MARK_FUNCTION_END` regions (Table 8.5).
- Host access: only the difference-quotient dense and band Jacobians touch raw vector data, on the
  user's subvectors.
- Naming: array parameters of public functions carry `_1d` or `_2d` suffixes and pointer outputs
  `_ptr` `[fact: doc/superbuild/source/developers/source_code/Naming.rst]`; no variadic public
  functions; callbacks carry `void* user_data`; C99 without variable-length arrays, so every per-stage
  scalar array has the fixed size `FIRK_MAX_STAGES`.

### 5.12 Smaller decisions

| ID  | Decision                                | Verdict                                                                                                                                                                                                                                                                                                                 |
| --- | --------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| D10 | return codes                            | explicit table (8.2): CVODE's values where the meaning is shared; the range $-13 \dots -19$ is left unused because CVODE's nonlinear-solver codes and ARKODE's mass codes collide there; FIRKODE-only codes in $-40 \dots -49$; `FIRKLS_*` equals ARKLS                                                                 |
| D11 | residual form                           | $G_i = M Z_i - h \sum_j a_{ij} F_j$, the form the preconditioner analysis assumes, not RADAU5's $A^{-1}$ form                                                                                                                                                                                                           |
| D12 | when $f(t_{n+1}, y_{n+1})$ is evaluated | at the start of the next step attempt under an `fn_current` flag; every evaluation of $f$ by the integrator goes through one internal wrapper, `firkRhs`, that counts calls in `nfe` and maps recoverable failures; the FIRKLS difference quotients count in `nfeDQ` and may use the `jtimes` RHS function, as in CVODE |
| D13 | rootfinding                             | ARKODE's `arkode_root.c` adapted with FIRKODE names on the collocation polynomial; a shared SUNDIALS rootfinding module is proposed for after 1.0                                                                                                                                                                       |
| D14 | reaction to a convergence failure       | reduce $h$ by the Newton driver's factor; re-evaluate $J$ only if it was stale; no retry at the same $h$                                                                                                                                                                                                                |

*Table 5.8. Smaller decisions.*

### 5.13 Defaults

| Parameter                                          | Default                                                                                                                                                                     |
| -------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| mode                                               | BEGINNER                                                                                                                                                                    |
| stages                                             | v0.1.0: $s = 3$; v0.2 EXPERT: any $1 \le s \le 9$, default 3; v0.3 BEGINNER: in situ, $3 \le s \le 7$                                                                       |
| $\epsilon_{\text{nls}}$, $m$, divergence threshold | RADAU5's rule $\max(10u/\mathrm{rtol}, \min(0.03, \sqrt{\mathrm{rtol}}))$, 0.03 without a relative tolerance (Section 4.4); 7; 0.99                                         |
| rate estimate, first-iteration exponent            | geometric mean of the last two ratios, 0.8                                                                                                                                  |
| safety                                             | $0.9$ on the controller's proposal; $\kappa = 0.9 \min\Bigl(1, \frac{1 + 2m}{k_{\text{nls}} + 2m}\Bigr)$ on the I-controller minimum (BEGINNER) or on the proposal (EXPERT) |
| step ratio limits                                  | $[0.2, 8]$; $10^4$ on the first step after `Init`, `ReInit` or `Reset`; 0.1 when it is rejected; 1 after a rejection; 0.3 after 2 error failures; 0.5 after a divergence    |
| dead band                                          | $[1, 1.2]$ when $\theta \le 10^{-3}$                                                                                                                                        |
| dsm floor                                          | $10^{-10}$                                                                                                                                                                  |
| controller                                         | `SUNAdaptController_ImpGus`, $k_1 = k_2 = 1$, with the I-controller minimum                                                                                                 |
| reuse policy                                       | `RADAU5`; `PERIODIC` uses 20 steps or $\Delta\gamma > 0.3$ for the setup and 51 steps for $J$                                                                               |
| $\epsilon_L$, $\epsilon_{stk}$, $\epsilon_{L,M}$   | 0.05, 0.05, 0.05                                                                                                                                                            |
| FGMRES                                             | $\mathrm{maxl} = \max(5, 3s)$, one restart, modified Gram–Schmidt                                                                                                           |
| refilter, predictor                                | on, extrapolation                                                                                                                                                           |
| linear-solution scaling                            | on for a matrix-based block solver; it acts only in the $s = 1$ Newton solve                                                                                                |
| failure limits                                     | 7 error-test failures, 10 convergence failures, 500 steps per `Evolve`, 10 warnings for $t + h = t$                                                                         |
| order rule (v0.3)                                  | raise at $\tilde\theta \le 0.002$, lower at $\tilde\theta \ge 0.8$, window $0.8 < \eta < 1.2$, hold 10 steps, veto $k_{lin} > 2s$                                           |

*Table 5.9. Default parameters.*

## 6. Comparison with other codes

### 6.1 Open-source Radau IIA and fully implicit Runge–Kutta codes

| Code                                | Methods                                                                              | Problem class                                                        | Stage linear algebra                                                               | Krylov, matrix-free                      | Parallelism                                | Step and order adaptivity                   | Licence     |
| ----------------------------------- | ------------------------------------------------------------------------------------ | -------------------------------------------------------------------- | ---------------------------------------------------------------------------------- | ---------------------------------------- | ------------------------------------------ | ------------------------------------------- | ----------- |
| RADAU5 [1996HW, 1999HW]             | Radau IIA, $s = 3$                                                                   | $M\dot y = f$, constant and possibly singular $M$, DAE index $\le 3$ | 1 real + 1 complex LU, dense or banded                                             | no                                       | serial                                     | step; fixed order                           | BSD-2-style |
| RADAU [1999HW]                      | $s \in \{1, 3, 5, 7\}$, orders 5, 9, 13 by default                                   | as RADAU5                                                            | 1 real + $\frac{s-1}{2}$ complex LU                                                | no                                       | serial                                     | step and order                              | BSD-2-style |
| SciPy `Radau` [2020V+]              | $s = 3$                                                                              | $\dot y = f$                                                         | 1 real + 1 complex LU                                                              | no                                       | serial                                     | step                                        | BSD-3       |
| OrdinaryDiffEq.jl [2017RN, 2025ESR] | RadauIIA3/5/9; `AdaptiveRadau`, odd $s$, arbitrary order                             | constant $M$, singular allowed                                       | 1 real + $\frac{s-1}{2}$ complex; any LinearSolve.jl solver                        | yes                                      | threads over the complex solves            | step and order                              | MIT         |
| PETSc `TSIRK`                       | Gauss–Legendre; user tableaus                                                        | $\dot y = f$                                                         | Newton on the monolithic $sN$ system                                               | any KSP/PC                               | MPI                                        | step through TSAdapt; no embedded estimator | BSD-2       |
| Irksome [2021FKM, 2025KM]           | Radau IIA (any $s$), Gauss, Lobatto, DIRKs                                           | finite-element forms with mass matrix; DAEs                          | monolithic stage-coupled system                                                    | GMRES with block or monolithic multigrid | MPI                                        | step                                        | LGPL-3      |
| PSIDE [1997HS]                      | Radau IIA, $s = 4$                                                                   | $g(t, y, \dot y) = 0$, index $\le 3$                                 | 4 independent real LUs                                                             | no                                       | 4-way shared memory                        | step                                        | not stated  |
| accelerInt [2017CNS]                | Radau IIA, $s = 3$                                                                   | $\dot y = f$, chemical kinetics                                      | dense LU per system                                                                | no                                       | CUDA, one ODE per thread                   | step                                        | MIT         |
| FIRKODE (proposed)                  | Radau IIA, $s = 3$ in v0.1.0; $1 \le s \le 9$ from v0.2; in situ odd order from v0.3 | $M\dot y = f$, constant nonsingular $M$                              | FGMRES with a single-shift block preconditioner; any SUNDIALS solver for the block | yes                                      | any `N_Vector`, MPI and GPU `[to measure]` | step; order from v0.3                       | BSD-3       |

*Table 6.1. Open-source Radau IIA and fully implicit RK codes, from their documentation and sources
as of September 2026 (see the software list after the references).*

deal.II, MFEM, Trilinos Tempus and diffrax offer only diagonally implicit RK methods; Boost.Odeint
offers implicit Euler and a Rosenbrock method; GSL offers implicit Gauss methods of orders 2 and 4;
Drake's `RadauIntegrator` has one or two stages. No general-purpose C library offers an adaptive
high-order Radau IIA method on user-supplied vectors and linear solvers.

### 6.2 Shared core and positioning

All adaptive Radau codes in Table 6.1 share RADAU5's core: simplified Newton with one $J$, a filtered
embedded estimate and a Gustafsson controller. FIRKODE keeps that core and replaces the complex
transformation. Among the codes surveyed:

- FIRKODE alone runs a stage-system Krylov solver on top of an arbitrary $N \times N$ user solver and
  any vector backend;
- FIRKODE and Irksome are the only MPI-capable codes with adaptive Radau IIA; PETSc `TSIRK` ships no
  embedded error estimator;
- RADAU5, RADAU and PSIDE handle a singular $M$ up to index 3; FIRKODE does not;
- spectral deferred correction frameworks iterate to the same collocation solution with
  stage-parallel sweeps [2000DGR, 2018S, 2025CLGR]; their sweep matrices are reachable in FIRKODE
  through the M6 hook (Section 5.2);
- no performance claim is made: iteration counts and run times are `[to measure]` (Section 11.4).

### 6.3 Deviations from RADAU5

The following differences from `radau5.f` are deliberate choices, listed because they change step
counts or make nominally equal runs different accuracy requests.

| Element                  | RADAU5                                                                                           | FIRKODE                                                                                                                  | Effect on comparisons                                                                                  |
| ------------------------ | ------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------ |
| tolerances               | replaces $\mathrm{rtol}$ by $0.1 \mathrm{rtol}^{\frac{2}{3}}$ and scales $\mathrm{atol}$ with it | uses the tolerances as given, as CVODE does                                                                              | equal nominal tolerances are different accuracy requests; the comparison harness of M2 transforms them |
| Newton tolerance         | $\max\Bigl(\frac{10u}{\mathrm{rtol}'}, \min(0.03, \sqrt{\mathrm{rtol}'})\Bigr)$                  | the same rule with the relative tolerance as given and $u$ = `SUN_UNIT_ROUNDOFF` (RADAU5: $10^{-16}$); EXPERT can fix it | none above $\mathrm{rtol} \approx 10^{-10}$                                                            |
| initial step             | fixed default $10^{-6}$                                                                          | CVODE's heuristic (Section 4.8)                                                                                          | the first steps differ                                                                                 |
| end of the integration   | lands on `XEND` with a shortened last step                                                       | interpolates at `tout` with the degree-$s$ collocation polynomial unless a stop time is set, which it then lands on      | the comparison harness sets the stop time                                                              |
| stage solve              | exact, complex LU                                                                                | inexact FGMRES to $\epsilon_{stk} \epsilon_{\text{nls}}$                                                                 | Newton counts may differ by the inexactness; see [2000J]                                               |
| error-test failures      | I-controller with safety                                                                         | same, plus the factor 0.3 after two consecutive failures                                                                 | rare                                                                                                   |
| norm                     | RMS with scale $\mathrm{atol} + \mathrm{rtol} \lvert y_n \rvert$                                 | WRMS with the same weights                                                                                               | identical                                                                                              |

*Table 6.2. Deviations from RADAU5. Everything else, including the rate estimate, the safety factor,
the dead band with* $\theta \le 10^{-3}$ *and the Jacobian policy, is adopted in BEGINNER mode.*

### 6.4 Choosing a SUNDIALS integrator

| Problem                                                                            | Package                                       | Reason                                                                                                                              |
| ---------------------------------------------------------------------------------- | --------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------- |
| non-stiff                                                                          | CVODE Adams, ERKStep                          | explicit or Adams methods are cheaper                                                                                               |
| stiff with a real spectrum and moderate accuracy                                   | CVODE BDF                                     | one implicit solve per step with a cheap factorization reuse                                                                        |
| stiff with an oscillatory spectrum, or tight tolerances                            | FIRKODE                                       | A- and L-stable at every order, stage order $s$, one-step                                                                           |
| stiff with a mass matrix and a stiffly accurate high-order method needed           | FIRKODE                                       | $M \dot y = f$ with the same block interface                                                                                        |
| stiff nonlinear dissipative, where contractivity of the numerical solution matters | FIRKODE                                       | algebraically stable, B-convergent of order $s$ (Section 4.2)                                                                       |
| many independent small stiff systems on a GPU                                      | CVODE, or FIRKODE with a batched block solver | one batched factorization per setup; FIRKODE adds order and L-stability for $s$ RHS evaluations per Newton iteration `[to measure]` |
| additively split stiff and non-stiff parts                                         | ARKStep IMEX                                  | implicit treatment of the stiff part only                                                                                           |
| index-1 DAE, or $F(t, y, \dot y) = 0$                                              | IDA                                           | consistent initialization and DAE error control                                                                                     |
| multiple time scales                                                               | MRIStep                                       | multirate infrastructure                                                                                                            |
| parabolic with explicit stabilization                                              | LSRKStep                                      | extended real stability without a solve                                                                                             |
| separable Hamiltonian                                                              | SPRKStep                                      | symplectic methods                                                                                                                  |

*Table 6.3. Choosing a SUNDIALS integrator.*

## 7. Costs, risks and maintenance

### 7.1 Cost model

Table 7.1 counts the work of one step in units that the user controls: RHS evaluations, block solves
with $M - \gamma J$, products with $J$ and $M$, setups, and vectors of length $N$ stored. The FIRKODE
entries are `[derived]` from the algorithms of Chapter 4; $k_{\text{nls}}$ Newton iterations and $k_{lin}$
FGMRES iterations per Newton iteration are parameters until measured.

| Quantity per step                      | FIRKODE, $s$ stages                                                                                                                                                                                     | CVODE BDF                    | ARKStep ESDIRK, $s_D$ implicit stages | RADAU5, $s = 3$                                                                                          |
| -------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------- | ------------------------------------- | -------------------------------------------------------------------------------------------------------- |
| RHS evaluations                        | $s k_{\text{nls}} + 1$ (+1 with refilter)                                                                                                                                                               | $k_{\text{nls}} + 1$         | $\approx s_D k_{\text{nls}} + 1$      | $3 k_{\text{nls}} + 1$ (+1 with refilter)                                                                |
| block solves with $M - \gamma J$       | $s k_{lin} k_{\text{nls}} + 1$ (+1 with refilter)                                                                                                                                                       | $k_{\text{nls}}$             | $s_D k_{\text{nls}}$                  | $k_{\text{nls}}$ real and $k_{\text{nls}}$ complex (about $5 k_{\text{nls}}$ real-solve equivalents) + 1 |
| products $J v$                         | $s k_{lin} k_{\text{nls}}$ (or RHS evaluations with difference quotients)                                                                                                                               | 0                            | 0                                     | 0                                                                                                        |
| products $M v$ ($M \ne I$)             | $s k_{\text{nls}} (1 + k_{lin})$ + 1                                                                                                                                                                    | not applicable               | $s_D k_{\text{nls}}$                  | folded into the transformation                                                                           |
| mass solves ($M \ne I$)                | $s k_{\text{nls}} + 1$ with an iterative block solver; $s k_{\text{nls}}$ ($s \ge 2$) or 0 ($s = 1$) with a direct one (Section 5.10)                                                                   | not applicable               | not counted here                      | none                                                                                                     |
| setups                                 | 1 real factorization of $M - \gamma J$ per setup                                                                                                                                                        | 1 per setup                  | 1 per setup                           | 1 real + 1 complex factorization (about 5 real LU equivalents) per setup                                 |
| vectors of length $N$ stored           | $11 + s(2\mathrm{maxl} + 10)$: 95 for $s = 3$, $\mathrm{maxl} = 9$                                                                                                                                      | $q + 1$ Nordsieck + about 10 | $s_D$ stage RHS + about 10            | about 12 plus the dense work array $4N^2$                                                                |
| global reductions per Krylov iteration | $j + 2$ dot products at Arnoldi step $j$, each one allreduce with MPIManyVector, $s$ with a plain ManyVector over MPI subvectors, which FIRKODE builds only when `y0` has no communicator (Section 8.3) | solver's                     | solver's                              | none (direct)                                                                                            |

*Table 7.1. Work per step. The complex-arithmetic equivalents count a complex multiply-add as four
real ones. A block solve is an exact solve with the factored block for a direct block solver and one
application of the user's preconditioner for an iterative one (Section 8.6).*

The break-even against RADAU5's dense path follows from flop counts `[derived]`: with a dense LU at
$\frac{2}{3} N^3$ flops, a triangular solve pair and a dense matrix-vector product at $2N^2$ each,
$s = 3$ and $k_{\text{nls}} = 2$:

| Setup frequency        | FIRKODE cheaper in flops when         | $k_{lin} = 4$   | $k_{lin} = 12$ (bound of Table 4.1)    |
| ---------------------- | ------------------------------------- | --------------- | -------------------------------------- |
| every step             | $N > 9 k_{lin} - 7.5$               | $N > 29$        | $N > 101$                              |
| every 20 steps         | $N > 180 k_{lin} - 150$             | $N > 570$       | $N > 2010$                             |

*Table 7.2. Break-even against a RADAU5-style complex LU for dense blocks. Memory: FIRKODE's*
$2N^2 + 95N$ *(two matrices and 95 vectors) is below RADAU5's* $4N^2 + 12N$ *for* $N > 41$.

The message is plain: for small dense systems under a CVODE-style reuse policy the complex
factorization is cheaper. FIRKODE's design pays off when setups are frequent (the RADAU5 policy),
when $N$ is in the hundreds or more, and whenever the block solver is iterative, matrix-free, batched
or on a GPU, where complex arithmetic is not available in SUNDIALS at all.

### 7.2 Risk register

| Risk                                                                                                | Likelihood        | Impact | Evidence today                                                                                                | Mitigation                                                                                                                                                                       | Gate                           |
| --------------------------------------------------------------------------------------------------- | ----------------- | ------ | ------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------ |
| FGMRES iteration counts for $s \ge 4$ are high                                                      | medium            | high   | Table 4.1 is a pessimistic bound; no measurement `[to measure]`                                               | measurement protocol of Section 11.4; $k_{lin}$ veto in the order rule; EXPERT stage preconditioners                                                                             | M4 gated on the measured table |
| Krylov defaults below the worst-case bound for $s \ge 6$                                            | high              | medium | $6s < k$ in Table 4.1                                                                                         | defaults revisited after M4; restarts configurable                                                                                                                               | M4                             |
| rank-local max norm and min, and one reduction per stage, on a plain ManyVector over MPI subvectors | high if unguarded | high   | `[fact]`, Table 5.5                                                                                           | automatic MPIManyVector; `Init` rejects a communicator without local reductions (a cost guard); documented `nvgetcommunicator` requirement; operation whitelist; mock tests (H9) | M1, M3                         |
| collective `MPI_Comm_dup` on every clone                                                            | medium            | medium | `[fact]`                                                                                                      | clones only at `Init` and order changes; a clone counter test (H14)                                                                                                              | M3                             |
| memory at high $s$: 587 vectors for $s = 9$, $\mathrm{maxl} = 27$                                   | medium            | medium | Table 8.7                                                                                                     | default range $3 \le s \le 7$; allocation at order changes only; `firk_test_workspace`                                                                                           | M4                             |
| single precision cannot verify $s \ge 5$, nor `CheckOrder` for $s = 4$                              | certain           | low    | $B(18)$ residual $9.4 \cdot 10^{-11}$ `[proved]`; $B(8)$ residual $1.0 \cdot 10^{-4}$ for $s = 4$ `[derived]` | `exclude-single` for $s \ge 5$; table tests for $s \le 4$, `CheckOrder` for $s \le 3$                                                                                            | M0                             |
| extended-precision literals and `-Wdouble-promotion`, `-Wconversion`                                | medium            | low    | CI matrix `[fact: .github/workflows]`                                                                         | `SUN_RCONST` literals with 40 digits; `sunindextype` throughout                                                                                                                  | M0                             |
| $\gamma_0$ for even $s$ is a heuristic                                                              | certain           | low    | Lean items R6, R7 open                                                                                        | EXPERT only; measured in M4                                                                                                                                                      | M4                             |
| Krylov-path answer files differ across platforms                                                    | high              | medium | iteration counts are platform-dependent                                                                       | exit-code tests for Krylov paths; answer files only for deterministic output                                                                                                     | M3                             |
| the RADAU5 comparison is confounded by Table 6.2                                                    | medium            | medium | deviations are enumerated                                                                                     | the harness transforms tolerances and uses the RADAU5 policy                                                                                                                     | M2                             |
| duplicated rootfinder diverges from ARKODE's                                                        | low               | low    | Table 5.2                                                                                                     | adapted verbatim; shared module proposed after 1.0                                                                                                                               | M2                             |
| bindings absent at first release                                                                    | certain           | medium | the developer checklist asks for SWIG and litgen runs                                                         | bindable API (Section 9.4); dated follow-up                                                                                                                                      | after v0.1.0                   |
| observed order on parabolic PDEs with time-dependent boundary data is about $s + 1$, not $2s - 1$   | certain           | low    | Table 4.3 [1993LO] `[cited]`                                                                                  | the guide states Table 4.3; the parabolic-order test of Table 11.2                                                                                                               | M2, M4                         |
| one external contributor                                                                            | certain           | high   | unfunded work (P1)                                                                                            | small milestones, complete documentation, independent versioning, explicit maintenance table                                                                                     | every release                  |

*Table 7.3. Risk register.*

### 7.3 Maintenance surfaces

| Path or process                                                                                                                                                                             | One-time                                                     | Per release                             |
| ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------ | --------------------------------------- |
| `src/firkode`, `include/firkode`                                                                                                                                                            | the package                                                  | fixes tracking CVODE and ARKODE         |
| `doc/firkode/guide` and the superbuild registration (`doc/superbuild/source/index.rst`, `conf.py`, `*_package_links.rst`, `doc/shared/sundials_vars.py`, `doc/shared/sundials/Install.rst`) | guide and links                                              | version strings                         |
| `examples/firkode`, `test/unit_tests/firkode`                                                                                                                                               | examples and tests                                           | answer-file refresh when output changes |
| `test/answers` (the `sundials-codes/answers` repository)                                                                                                                                    | answer files for three precisions                            | same                                    |
| `scripts/firkode.sh`, `scripts/tarscript.sh`, `scripts/updateVersion.sh`, `scripts/startReleaseCycle.sh`                                                                                    | tarball list and version plumbing                            | version bump                            |
| `CHANGELOG.md`, `doc/shared/RecentChanges.rst`, `CITATIONS.md`, `README.md`                                                                                                                 | announcement                                                 | entries                                 |
| `scripts/spack/packages/sundials/package.py`                                                                                                                                                | a `firkode` variant                                          | none                                    |
| `swig/`, `bindings/sundials4py`                                                                                                                                                             | after the API settles                                        | regeneration                            |
| `cmake/SundialsBuildOptionsPre.cmake`, `SundialsBuildOptionsPost.cmake`                                                                                                                     | `SUNDIALS_ENABLE_FIRKODE`; the ManyVector dependency in Post | none                                    |
| CI                                                                                                                                                                                          | none, if the conventions of Chapter 10 are followed          | none                                    |
| `verification/lean` (optional)                                                                                                                                                              | generator target, pinned toolchain, copyright handler        | toolchain pin                           |

*Table 7.4. Maintenance surfaces.*

### 7.4 Effort by analogy

CVODE's core is about 9300 lines (`cvode.c` 5003, `cvode_ls.c` 1976, `cvode_io.c` 1841,
`cvode_nls.c` 490 `[fact]`), plus ARKODE's rootfinder (961) and option parser (328). FIRKODE's
v0.1.0 is of the same order of magnitude, plus tests, examples and a guide. The milestones of Chapter
10 are sized so that each is reviewable on its own (P1).

### 7.5 What is asked of the maintainers, and what is offered

#### What is asked

Review time for four milestones, CI minutes within the existing matrix, access to the answers
repository, and the release-checklist items of Table 7.4. The package is labelled experimental in the
changelog, is versioned independently of the suite, and stays behind `SUNDIALS_ENABLE_FIRKODE`.

#### What is offered

A stiffly accurate, L-stable, high-order integrator on SUNDIALS' vectors and solvers; a
stiff validation suite with stated pass criteria that other packages can reuse; the table generator
and its certificates as tooling; and a design document that is kept current with every milestone.
This is unfunded work by one contributor (P1); the plan is therefore conservative about scope and
explicit about what remains unverified.

## 8. Architecture

This chapter is the implementer's specification. Pseudocode is normative for control flow and for the
order of operations; names of internal functions are suggestions, names of public functions and
constants are binding (Chapter 9).

### 8.1 Package shape

- `libsundials_firkode` is built from `sundials_core` and the object libraries
  `sundials_nvecmanyvector_obj`, `sundials_sunlinsolspfgmr_obj` and
  `sundials_sunadaptcontrollersoderlind_obj`, as ARKODE is `[fact: src/arkode/CMakeLists.txt]`, plus
  `sundials_nvecmpimanyvector_obj` under `SUNDIALS_ENABLE_MPI`, and, as the other packages, the
  convenience objects of the serial vector, the dense, band and sparse matrices, and the dense,
  band and Krylov linear solvers (R14); no nonlinear-solver or ImExGus objects. The CMake option
  `SUNDIALS_ENABLE_FIRKODE` is defined with the other packages in `cmake/SundialsBuildOptionsPre.cmake`;
  the ManyVector options are defined later, in `cmake/SundialsBuildOptionsPost.cmake` `[fact]`, so the
  dependency is enforced there: with FIRKODE on, `SUNDIALS_ENABLE_NVECTOR_MANYVECTOR` (and
  `SUNDIALS_ENABLE_NVECTOR_MPIMANYVECTOR` under MPI) is forced on with a status message.
- No header from `include/arkode` or `include/cvode` is included. Callback typedefs are redeclared
  with identical signatures (Table 9.2); C function-pointer compatibility lets users pass the same
  functions they pass to CVODE or ARKODE.
- Vocabulary: a *block* is an $N \times N$ object or a vector of length $N$; a *stacked* object has
  $sN$ rows. The user's `SUNLinearSolver` and `SUNMatrix` are block objects. The internal SPFGMR is
  the stacked solver.

### 8.2 Memory objects

| Object (header)                                  | Field groups                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| ------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `FIRKodeMem` (`firkode_impl.h`)                  | problem: `f`, `user_data`, `sunctx`, `mode`; time and step: `tn`, `h`, `hprime`, `hu`, `h0u`, `hmin`, `hmax_inv`, `tstop`, `tstopset`, `tstopinterp`, `fixedstep`; tolerances: `itol`, `reltol`, `Sabstol`, `Vabstol`, `efun`, `e_data`, `ewt`; vectors: `yn`, `fn`, `fn_current`, `tempv1..3`, `err`; method: `s`, `tbl` (pointer to a `FIRKodeTableMem`), `smin`, `smax`; stack: `Z`, `F`, `G`, `delta`, `Jv`, `W` (weight view), `stack_kind`; Newton: `maxcor`, `nlscoef_user`, `nlscoef`, `div_rate`, `rate_pow`, `theta_jac`, `jac_bad_rate`, `eta_prev`, `nls_status`, `nls_iters`, `nls_theta`, `nls_hfac`, hook pointers; adaptivity: `C`, `own_C`, `safety`, `newton_safety`, `etamx1`, `growth`, `etamin`, `etamxf`, `small_nef`, `etacf`, `lbound`, `ubound`, `refilter`, `dsm_floor`, `dsm`, `eta`, `attempt_rejected`, `predictor`; reuse: `policy`, `policy_fn`, `policy_data`, `msbp`, `dgmax_lsetup`, `nstlp`, `gamma`, `gammap`, `gamrat`, `jcur`, `keep_setup`, `record`; history: `hist = {t, h, tbl, Z, valid}`; limits and counters: `mxstep`, `mxhnil`, `maxnef`, `maxncf`, `nst`, `nst_attempts`, `nfe`, `ncfn`, `nni`, `nnf`, `netf`, `nsetups`, `nhnil`, `tolsf`, `lrw`, `liw`; interfaces: `lmem`, `mass_mem`, `mass_iterative`, `root_mem`, `order_mem` (v0.3) |
| `FIRKLsMem` (`firkode_ls_impl.h`)                | block solver: `LS`, `A`, `savedJ`, `iterative`, `matrixbased`, `scalesol`, `x`, `ytemp`; Jacobian: `jac`, `user_jac`, `J_data`, `linsys`, `user_linsys`, `jtsetup`, `jtimes`, `jt_data`, `jtimesRhs`, `jbad`, `msbj`, `dgmax_jbad`, `nstlj`, `tnlj`; preconditioner: `pset`, `psolve`, `P_data`; tolerance: `eplifac`, `nrmfac`, `sysfac`; stacked solver: `SPFGMR`, `maxl`, `maxrs`, `gstype`, `eps_stk`, `x_stk`; counters: `nje`, `nfeDQ`, `npe`, `nps`, `nli`, `ncfl`, `njtsetup`, `njtimes`, `nsli`, `nscf`, `nblocksolves`, `nblockfails`, `last_flag`; hooks (v0.4): `stage_pset`, `stage_psolve`, `stage_LS`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| `FIRKLsMassMem` (`firkode_ls_impl.h`)            | `LS`, `M`, `M_lu`, `mass`, `M_data`, `time_dependent` (must be false), `mtsetup`, `mtimes`, `mt_data`, `pset`, `psolve`, `P_data`, `eplifac`, `nrmfac`, `x`, `msetuptime`; counters: `nmsetups`, `nmsolves`, `nmmults`, `nmtsetup`, `nmpe`, `nmps`, `nmli`, `nmcfl`, `last_flag`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| `FIRKodeRootMem` (`firkode_root_impl.h`)         | `gfun`, `nrtfn`, `iroots`, `rootdir`, `gactive`, `glo`, `ghi`, `grout`, `tlo`, `thi`, `trout`, `toutc`, `ttol`, `taskc`, `irfnd`, `nge`, `mxgnull`; the user data pointer is read from `FIRKodeMem` at call time                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| `FIRKodeTableMem` (public, `firkode_tables.h`)   | `s`, `q`, `p`, `c[]`, `A[]` (row-major, $a_{ij}$ = `A[i * s + j]`), `b[]`, `Ainv[]`, `d[]`, `e[]`, `gamma0`, `P[]` (dense-output basis coefficients, `P[k * s + i]`), `stiffly_accurate`, `owned`; arrays of the fixed size `FIRK_MAX_STAGES` or its square; the built-in tables are `static const` objects, returned through the `FIRKodeTable` handle read-only and not owned (`_Free` ignores them, `_Copy` gives a modifiable owned table), and allocate nothing                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| `FIRKodeOrderMem` (v0.3, `firkode_order_impl.h`) | `theta_filt`, `knls_last`, `klin_avg`, `eta_last`, `flag_newton_fail`, `flag_unexpected_reject`, `steps_since_change`, `steps_since_eval`, `nchanges`, `nst_per_order[FIRK_MAX_STAGES + 1]`, `W[FIRK_MAX_STAGES + 1]`, cost samples `C_f`, `C_sol`, `C_Jv`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |

*Table 8.1. Memory objects and their field groups. Every per-stage scalar array has the fixed size
`FIRK_MAX_STAGES` (C99, no variable-length arrays).*

### 8.3 The stage stack

#### Kind

`N_VNew_ManyVector(s, vecs, sunctx)` over $s$ clones of `y0`, or, when `SUNDIALS_ENABLE_MPI`
is defined and `N_VGetCommunicator(y0) != SUN_COMM_NULL`, `N_VMake_MPIManyVector(comm, s, vecs,
sunctx)` `[fact: include/nvector/nvector_mpimanyvector.h]`. A NULL return is a `FIRK_MEM_FAIL`.
`Init` decides the kind (`firkStackSetKind` in `firkode_stack.c`). Before choosing the MPI kind it
checks that `y0->ops->nvdotprodlocal` and `nvwsqrsumlocal` exist and returns `FIRK_ILL_INPUT`
otherwise: without them the MPI kind calls each subvector's global reduction and then reduces
across ranks, correct but no cheaper than the plain kind (D4). A distributed `y0` without
`nvgetcommunicator` cannot be detected and is excluded by documentation. For $s = 1$ there is no
stack; the accessor returns the vector itself.

#### Accessor

One internal function `firkStage(mem, stk, i)` dispatches on the kind recorded in `Init` and
returns the $i$-th subvector; the vector ID cannot tell the kinds, or a stack from a ManyVector
`y0`, apart. All stage-wise work (predictor, residual, estimate, operator products) runs on
subvectors of length $N$; only SPFGMR and the stacked norm see the stacked object.

#### Weight view

The stacked weight vector `W` is a ManyVector whose $s$ subvector pointers all point
at `ewt`. It is created without data ownership, never cloned, and destroyed without freeing `ewt`
(H16). It is rebuilt at an order change.

#### Clones

Clones happen only in `FIRKodeInit` and at an order change (H14), since every clone of an
MPIManyVector performs a collective `MPI_Comm_dup` `[fact: src/nvector/manyvector/nvector_manyvector.c]`.

#### Permitted operations and vector budget

| Operation on a stacked vector                                                            | Status      | Reason                                                                                                                                                                                                                         |
| ---------------------------------------------------------------------------------------- | ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `N_VConst`, `N_VScale`, `N_VLinearSum`, `N_VProd`, `N_VDiv`, `N_VDotProd`, `N_VWrmsNorm` | permitted   | correct on both kinds: one allreduce on the MPI kind, one per stage on the plain kind (Table 5.5)                                                                                                                              |
| `N_VDotProdMulti`, `N_VLinearCombination`                                                | SPFGMR only | the plain ManyVector allocates a host pointer array per call; FIRKODE's own code uses subvector loops                                                                                                                          |
| `N_VMaxNorm`, `N_VMin`, `N_VInvTest`, `N_VConstrMask`, `N_VMinQuotient`                  | forbidden   | not needed by the stage solver; the whitelist stays minimal so that the H9 mock traps everything else, and on the plain kind they are rank-local over distributed stages `[fact: src/nvector/manyvector/nvector_manyvector.c]` |
| `N_VGetLocalLength`, `N_VGetArrayPointer`                                                | forbidden   | no such operation on a ManyVector; `N_VGetLength` is the sum of the subvector lengths                                                                                                                                          |

*Table 8.8. Permitted and forbidden operations on the stack. A unit test with a mock subvector whose
forbidden operations trap enforces the table (H9).*

| Owner            | Vectors of length $N$                                                    | Count                                                                      |
| ---------------- | ------------------------------------------------------------------------ | -------------------------------------------------------------------------- |
| integrator       | `yn`, `ycur`, `fn`, `ewt`, `ele`, `ftemp`, `tempv1`, `tempv2`, `tempv3`  | 9, and `Vabstol` with vector absolute tolerances                           |
| history          | `Zprev`                                                                  | $s$                                                                        |
| Newton           | `Z`, `F`, `G` ($b$ formed in place), `delta`                             | $4s$                                                                       |
| stacked operator | `Jv`                                                                     | $s$                                                                        |
| SPFGMR           | $2\mathrm{maxl} + 4$ stacked vectors                                     | $s(2\mathrm{maxl} + 4)$                                                    |
| FIRKLS block     | `x`, `ytemp`                                                             | 2                                                                          |
| mass             | `x`, `ydn` (`ydn` aliases `fn` without a mass matrix)                    | 2                                                                          |
| weight view      | pointers only                                                            | 0                                                                          |
| total            | $11 + s(2\mathrm{maxl} + 10)$ for $s \ge 2$, 16 for $s = 1$; +2 with $M$ | 95 for $s = 3$, $\mathrm{maxl} = 9$; 587 for $s = 9$, $\mathrm{maxl} = 27$ |

*Table 8.7. Vector budget, checked by `firk_test_workspace`, which counts the live clones of
`y0` (H27). SUNDIALS deprecates the work-space getters (removal in 8.0.0), so FIRKODE has no
`GetWorkSpace` (R20). For* $s = 1$ *the stacked rows are absent.*

### 8.4 The evolve loop and one step

```
FIRKodeEvolve(mem, tout, yout, &tret, itask):
  check: mem, Init done, yout and tret not NULL, itask in {FIRK_NORMAL, FIRK_ONE_STEP}
  if first:                                     # first call since Init, ReInit or Reset
    if |tout - tn| < 2u max(|tn|, |tout|): return FIRK_TOO_CLOSE   # also with an initial step, as CVODE
    check: tolerances set, linear solver attached
    ewt <- weights(y0); a zero weight -> FIRK_ILL_INPUT
    FIRKLS linit; mass linit and one setup of the constant M
    fn <- firkRhs(t0, y0); unrecoverable -> FIRK_RHSFUNC_FAIL, recoverable -> FIRK_FIRST_RHSFUNC_ERR
    if tstop is set and not beyond t0 towards tout: return FIRK_ILL_INPUT
    h <- fixed step if set, else SetInitStep value if set, else firkHin(tout, or tstop if nearer)   # ydot = M^-1 f by one mass solve
    h pointing away from tout -> FIRK_ILL_INPUT; h <- within [hmin, hmax]
    if the step would pass tstop: h <- (tstop - tn) (1 - 4u)
    rootfinding: RootCheck1
  else:
    rootfinding: RootCheck2, then RootCheck3 over the rest of the last step -> FIRK_ROOT_RETURN
    if tn is at tstop and tout is not before it: tret <- tstop; yout <- interpolated or copied;
                                                  clear tstop; return FIRK_TSTOP_RETURN
    if the next step would pass tstop: h <- (tstop - tn) (1 - 4u)
    if itask == FIRK_NORMAL and (tn - tout) h >= 0: tret <- tout; yout <- GetDky(tout, 0); return FIRK_SUCCESS
    if itask == FIRK_ONE_STEP and tn was not returned yet: tret <- tn; yout <- yn; return FIRK_SUCCESS
  loop:
    if steps in this call >= mxstep: return FIRK_TOO_MUCH_WORK
    ewt <- weights(yn), except on the first step; a zero weight -> FIRK_ILL_INPUT
    if u * ||yn||_W > 1: tolsf <- 2 u ||yn||_W; return FIRK_TOO_MUCH_ACC   # CVODE's factor 2
    if tn + h == tn: warn, at most mxhnil times
    if not fn_current: fn <- firkRhs(tn, yn); unrecoverable -> FIRK_RHSFUNC_FAIL, recoverable -> FIRK_UNREC_RHSFUNC_ERR   # D12
    (v0.3) order selection if due
    kflag <- firkStep()                                    # Listing 8.2
    if kflag < 0: return firkHandleFailure(kflag)
    rootfinding: RootCheck3 -> FIRK_ROOT_RETURN
    if (tn - tstop) h >= 0: tret <- tstop; yout <- interpolated or copied per SetInterpolateStopTime;
                            clear tstop; return FIRK_TSTOP_RETURN
    if itask == FIRK_ONE_STEP: tret <- tn; yout <- yn; return FIRK_SUCCESS
    if (tn - tout) h >= 0: tret <- tout; yout <- GetDky(tout, 0); return FIRK_SUCCESS
    if the next step would pass tstop: h <- (tstop - tn) (1 - 4u)
```

*Listing 8.1. The evolve loop, after `CVode` `[fact: src/cvode/cvode.c]`. The stop-time semantics
are CVODE's: `SetStopTime` after the time has passed is `FIRK_ILL_INPUT`; a `FIRK_TSTOP_RETURN`
clears the stop time; in fixed-step mode the last step is shortened and the fixed step resumes
afterwards (H21). `first` is set by `Init`, `ReInit` and `Reset` and cleared when the first step is
completed, as ARKODE's `initsetup` `[fact: src/arkode/arkode.c]`; CVODE tests `nst == 0`
`[fact: src/cvode/cvode.c]`, which `Reset` does not restore because it keeps the counters. "The
first step" means this step everywhere in this document. The `FIRK_TOO_CLOSE` test and the factor
2 on `tolsf` are CVODE's `[fact: src/cvode/cvode.c, cvInitialSetup and CVode]`.*

```
firkStep():
  nflag <- FIRST_CALL; ncf <- 0; nef <- 0; attempt_rejected <- FALSE
  loop:
    nst_attempts += 1; log begin-step-attempt               # one end-step-attempt per attempt, on every exit (H31)
    gamma <- h * gamma0
    (jbad, need_setup) <- policy(nflag, record)            # D9
    firkPredict()                                           # Listing 8.7
    if need_setup:
      flag <- firkLsSetup(jbad, &jcur)                      # Listing 8.4
      if flag < 0: return FIRK_LSETUP_FAIL
      if flag > 0: nls <- FIRK_NLS_LSETUP_RECVR; nls_hfac <- etacf; goto handle   # Table 8.3, R5
    nls <- firkNewton()                                     # Listing 8.3
  handle:
    kflag <- firkHandleNFlag(nls, &ncf)
    if kflag == PREDICT_AGAIN: nflag <- PREV_CONV_FAIL; attempt_rejected <- TRUE; continue   # RADAU5's REJECT
    if kflag < 0: return kflag
    dsm <- 0
    if not fixedstep:
      flag <- firkErrEstimate(&dsm)                         # Listing 8.5
      if flag > 0: nls <- FIRK_NLS_LSOLVE_RECVR (FIRK_NLS_RHS_RECVR for the refilter RHS); nls_hfac <- etacf; goto handle
      if flag < 0: return FIRK_LSOLVE_FAIL
      kflag <- firkCheckTemporalError(&nflag, &nef, dsm)    # Listing 8.6; stores eta and keep_setup
      if kflag == TRY_AGAIN: attempt_rejected <- TRUE; continue
      if kflag < 0: return kflag
    firkCompleteStep(dsm)
    return FIRK_SUCCESS

firkHandleNFlag(nls, &ncf):
  if nls == FIRK_NLS_SUCCESS: return FIRK_SUCCESS
  ncfn += 1
  if nls < 0: return nls                                    # FIRK_RHSFUNC_FAIL, FIRK_LSOLVE_FAIL, FIRK_MASSMULT_FAIL, FIRK_VECTOROP_ERR
  ncf += 1; etamax <- 1                                     # no growth on the next accepted step
  too_many <- fixedstep ? (ncf > 1 or jcur) : (ncf == maxncf or |h| <= hmin_eff (1 + u))   # hmin_eff = max(hmin, 10 u |tn|) (R18)
  if too_many: return (nls == FIRK_NLS_RHS_RECVR) ? FIRK_REPTD_RHSFUNC_ERR : FIRK_CONV_FAILURE
  if not fixedstep: eta <- max(nls_hfac, hmin_eff / |h|); h <- eta * h   # Table 8.3, D14
                                                            # fixed step (R3): retried once, unreduced, when J was stale
  record.flag_newton_fail <- TRUE
  return PREDICT_AGAIN

firkCompleteStep(dsm):
  if not fixedstep: SUNAdaptController_UpdateH(C, h, max(1e-2, dsm) if default_controller else dsm)  # eta computed before UpdateH (H36)
                    failure -> FIRK_CONTROLLER_ERR, the state of the last step intact
  yn <- yn + sum_j d_j Z_j                                  # in place, one linear combination; Radau IIA: d = e_s
  tn <- tn + h, snapped to tstop within roundoff; hu <- h; nst += 1
  history: swap Z and hist.Z when s is unchanged, else replace; hist <- {tn - h, h, tbl, valid}
  fn_current <- FALSE                                       # D12
  h <- fixedstep ? the fixed step : eta * h                 # the fixed step resumes after a step shortened to tstop (H21)
  record <- {nls_theta, nls_iters, k_lin of this step, eta, flags}; clear the flags
```

*Listing 8.2. One step attempt. `nflag`, `ncf` and `nef` are local to the step; `nst_attempts`
counts every attempt, at its start.*

| Outcome of the attempt                                                                             | $h$                                                | Setup at the next attempt            | $J$ re-evaluated                                           |
| -------------------------------------------------------------------------------------------------- | -------------------------------------------------- | ------------------------------------ | ---------------------------------------------------------- |
| `FIRK_NLS_DIVERGED`, `FIRK_NLS_MAXITER`                                                            | $\times 0.5$                                       | yes ($\gamma$ changed)               | only if it was not current at the failed setup             |
| `FIRK_NLS_PREDICTED_FAIL`                                                                          | $\times 0.8 q^{-\frac{1}{4+m-1-k}}$                | yes                                  | only if not current                                        |
| `FIRK_NLS_RHS_RECVR`, `FIRK_NLS_LSOLVE_RECVR`, `FIRK_NLS_LSETUP_RECVR`, recoverable estimate solve | $\times$ `etacf` (0.25, `SetMaxCFailGrowth`)       | yes                                  | only if not current                                        |
| any Newton or setup failure in fixed-step mode                                                     | $\times 1$, retried once if $J$ was not current    | yes                                  | yes; a second failure is `FIRK_CONV_FAILURE`               |
| error test failed, fewer than 2 consecutive failures                                               | $\times \eta$ from the controller                  | yes                                  | RADAU5: only if not current; PERIODIC: if due              |
| error test failed, 2 or more consecutive failures                                                  | $\times \min(\eta, 0.3)$; first step: $\times 0.1$ | yes                                  | RADAU5: only if not current; PERIODIC: if due              |
| accepted after a rejected attempt of the same step                                                 | $\times \min(\eta, 1)$                             | per policy                           | per policy                                                 |
| accepted, inside the dead band                                                                     | $\times 1$                                         | no                                   | no                                                         |
| accepted, outside the dead band                                                                    | $\times \eta$                                      | RADAU5 policy: yes; PERIODIC: if due | RADAU5: unless $\theta \le \theta_{jac}$; PERIODIC: if due |

*Table 8.3. Step-driver outcomes and actions; the convergence-failure rows implement D14.*

### 8.5 The Newton driver

```
firkNewton():            # Z holds the predictor on entry; returns FIRK_NLS_*; sets nls_iters, nls_theta, nls_hfac
  eta <- max(eta_prev, u)^p (p = 0.8, SetNonlinRateExponent); theta <- theta_jac (1e-3); dynold <- 0; ratio_prev <- 0
  status <- FIRK_NLS_MAXITER; hfac <- 0.5
  for k = 1 .. m:
    for i = 1..s: F_i <- firkRhs(tn + c_i h, yn + Z_i)       # nfe += s; recoverable -> status FIRK_NLS_RHS_RECVR, hfac etacf (R5), stop
    for i = 1..s: b_i <- -M Z_i + h sum_j a_ij F_j            # b = -G(Z), residual form D11; mass product or Z_i; one linear combination of s + 1 subvectors into G_i
    flag <- (s == 1) ? firkLsSolveBlock(b_1, delta_1, NEWTON, k) : firkLsSolveStacked(b, delta, k)
    nni += 1                                                   # every linear solve counts (R4)
    if flag > 0: status <- FIRK_NLS_LSOLVE_RECVR; hfac <- etacf; stop
    if flag < 0: return FIRK_LSOLVE_FAIL
    dyno <- ||delta||_W                                        # WRMS over all sN components
    if dyno is not finite: status <- FIRK_NLS_DIVERGED; hfac <- 0.5; stop   # R4
    if 2 <= k < m:
      ratio <- dyno / dynold
      theta <- (k == 2) ? ratio : sqrt(ratio * ratio_prev); ratio_prev <- ratio
      if theta >= 0.99: status <- FIRK_NLS_DIVERGED; hfac <- 0.5; break
      eta <- theta / (1 - theta)
      dyth <- eta * dyno * theta^(m - 1 - k) / eps_nls
      if dyth >= 1: q <- clamp(dyth, 1e-4, 20); hfac <- 0.8 * q^(-1 / (4 + m - 1 - k));
                    status <- FIRK_NLS_PREDICTED_FAIL; break
    dynold <- max(dyno, u)
    Z <- Z + delta
    log nonlinear-iterate (k, nni, dyno, theta, eta)
    if eta * dyno <= eps_nls: status <- FIRK_NLS_SUCCESS; break
  eta_prev <- eta
  nls_iters <- min(k, m); nls_theta <- theta; nls_hfac <- hfac
  if status == FIRK_NLS_SUCCESS: jcur <- FALSE else: nnf += 1
  return status
```

*Listing 8.3. The Newton driver with RADAU5's test (Section 4.4). The v0.4 hooks `SetNlsConvTestFn`
and `SetNlsNormFn` replace the test block and the norm; a custom test must set `hfac` on failure.*

#### Contract

The Newton driver never changes `h`. It returns the status and writes `nls_hfac` in
$(0, 1]$; `firkHandleNFlag` is the only consumer. The two failure paths that bypass the driver, a
recoverable setup failure and a recoverable estimate solve, write `nls_hfac` $=$ `etacf` themselves
before entering the handler (Listing 8.2, Table 8.3). `etacf` is the factor of every recoverable
failure of a function the integrator calls (R5); the factors 0.5 and the predicted-failure formula
are RADAU5's. The `FIRK_NLS_*` values are a private positive
enum in `firkode_impl.h`; negative values are public `FIRK_*` codes passed through unchanged.

### 8.6 FIRKLS: the linear-solver interface

FIRKLS owns the user's block solver and matrix, the saved Jacobian, the stacked SPFGMR, the mass
interface and the difference-quotient routines. It does not decide when to set up; the reuse policy
does (D9).

```
firkLsSetup(jbad, &jcur):
  gamma <- h * gamma0
  if matrixbased:
    if user linsys: linsys(tn, yn, fn, A, M, jok = not jbad, &jcur, gamma, user_data, tmp1, tmp2, tmp3)
    else if jbad: SUNMatZero(A) for direct solvers; jac(tn, yn, fn, A, user_data, tmp1, tmp2, tmp3);
                  savedJ <- A; nje += 1; nstlj <- nst; jcur <- TRUE
    else:         A <- savedJ; jcur <- FALSE
    A <- (M attached) ? SUNMatScaleAdd(-gamma, A, M) : SUNMatScaleAddI(-gamma, A)
  else:
    jcur <- jbad                                           # the user psetup receives jok = not jbad and may change jcur
  flag <- SUNLinSolSetup(LS, A); npe += 1 if a preconditioner setup ran
  gammap <- gamma; nsetups += 1; nstlp <- nst
  return 0 | recoverable (FIRKLS_JACFUNC_RECVR, SUNLS_PSET_FAIL_REC) | unrecoverable (FIRKLS_SUNMAT_FAIL, FIRKLS_JACFUNC_UNRECVR, FIRKLS_SUNLS_FAIL)

firkLsSolveBlock(b, x, purpose, k):                        # one N x N solve with the current M - gamma J
  x0 <- b; sigma <- 1
  if iterative and M attached and purpose != PRECOND:      # units of y (Section 5.10)
    x0 <- M^{-1} b (mass solve, tol eps_LM ||b||_W); a failure fails this solve with its sign; sigma <- ||x0||_W / ||b||_W
  tol_wrms by purpose (Table 5.7, with ||x0||); if iterative: SUNLinSolSetScalingVectors(LS, ewt, ewt); delta <- tol_wrms * nrmfac
  gamma_op <- (purpose == NEWTON or no setup ran) ? gamma : gammap                 # the gamma of firkLsATimes and of the user's psolve (Table 8.6)
  if iterative and purpose == NEWTON and k >= 2 and ||x0||_W <= tol_wrms: x <- 0; return 0   # CVLS's shortcut, iterative solvers only; never for ESTIMATE or PRECOND (H17)
  if sigma != 1: matrix-free solver: solve sigma (M - gamma J) x = sigma b (ATimes, PSolve and b scaled); matrix-based: delta <- delta / sigma
  if a jtsetup is attached and (purpose == NEWTON or (purpose == ESTIMATE and not matrixbased)):
    jtsetup(tn, yn, fn, user_data); njtsetup += 1                                   # after the shortcut, as in CVLS (H6)
  x <- 0; SUNLinSolSetZeroGuess(LS, TRUE)                                           # Krylov solvers return early without writing x; the flag resets after every solve (H11)
  flag <- SUNLinSolSolve(LS, A, x, b, delta); nli, nps, ncfl updated from the solver
  if purpose == NEWTON and matrixbased and linear-solution scaling on: x <- x * 2 / (1 + gamma / gammap)
  map: 0 -> 0;  SUNLS_RES_REDUCED -> (purpose == NEWTON and k == 1) ? 0 : +1 (H18);
       SUNLS_RECOV_FAILURE, SUNLS_CONV_FAIL, SUNLS_ATIMES_FAIL_REC, SUNLS_PSOLVE_FAIL_REC, SUNLS_PSET_FAIL_REC,
       SUNLS_PACKAGE_FAIL_REC, SUNLS_QRFACT_FAIL, SUNLS_LUFACT_FAIL -> +1;
       SUNLS_GS_FAIL, SUNLS_QRSOL_FAIL, SUNLS_ATIMES_FAIL_UNREC, SUNLS_PSOLVE_FAIL_UNREC, SUNLS_PSET_FAIL_UNREC,
       SUNLS_ATIMES_NULL, SUNLS_PSOLVE_NULL, SUN_ERR_EXT_FAIL (R8: there is no SUNLS_PACKAGE_FAIL_UNREC),
       SUN_ERR_ARG_*, SUN_ERR_MEM_FAIL and any other value -> -1
  for purpose PRECOND a +1 result is accepted and counted in nblockfails ("Block prec fails" in PrintAllStats)

firkLsSolveStacked(b, x, k):                                # s >= 2: the Newton system
  if a jtsetup is attached: jtsetup(tn, yn, fn, user_data) exactly once per Newton iteration; njtsetup += 1 (H6)
  SUNLinSolSetScalingVectors(SPFGMR, W, W); delta <- eps_stk * eps_nls * sqrt(s) * nrmfac
  if M attached: x0_i <- M^{-1} b_i for every stage (mass solves, tol eps_LM ||b||_W; a failure fails this solve with its sign)
                 delta <- delta * ||b||_W / ||x0||_W    # units of y (Section 5.10); exact for the true residual of right preconditioning
  x <- 0; SUNLinSolSetZeroGuess(SPFGMR, TRUE)              # no small right-hand side shortcut: SPFGMR's zero-guess return gives x = 0
  flag <- SUNLinSolSolve(SPFGMR, NULL, x, b, delta); nsli, nscf updated
  same mapping as above; SUNLS_RES_REDUCED accepted only when k == 1

firkLsATimes(v, z):                                         # z = (I_s (x) M - h A (x) J) v
  if savedJ and the matrix provides it: SUNMatMatvecSetup(savedJ), once per Jacobian evaluation, not per call
  for j: Jv_j <- savedJ ? SUNMatMatvec(savedJ, v_j) : jtimes attached ? jtimes(v_j, Jv_j, tn, yn, fn, user_data, tmp) :
                 difference quotient: sig <- ||v_j||_W; if sig == 0: Jv_j <- 0 (H7) else sig <- 1/sig;
                 Jv_j <- (firkRhs(tn, yn + sig v_j) - fn) / sig, up to 3 tries with sig <- sig/4 on a recoverable failure; nfeDQ += 1
  for i: z_i <- M v_i - h sum_j a_ij Jv_j                     # one linear combination of s + 1 subvectors
  return 0 | SUNLS_ATIMES_FAIL_REC | SUNLS_ATIMES_FAIL_UNREC

firkLsPSolve(r, z, tol, lr):                                # lr is SUN_PREC_RIGHT; tol is ignored (relative per block)
  for i:
    if matrixbased:          firkLsSolveBlock(r_i, z_i, PRECOND, -)          # direct: exact solve with the factored block;
                                                            # matrix-iterative: one inner solve to eps_L ||r_i||_W (R6)
    else if psolve attached and lr_user != NONE:
                             psolve(tn, yn, fn, r_i, z_i, gamma_p, eps_L ||r_i||_W nrmfac, lr_user, user_data); nps += 1
                             nblocksolves += 1              # the user's preconditioner only; the user's Krylov iteration never runs here (D2, H37);
                                                            # lr_user is the user's pretype; SUN_PREC_BOTH applies left then right (2 psolves)
    else:                    z_i <- r_i                     # no preconditioner: unpreconditioned FGMRES; not a block solve
  a recoverable failure is accepted and counted in nblockfails; an unrecoverable one returns nonzero,
  which SPFGMR reports as SUNLS_PSOLVE_FAIL_UNREC
  lr_user (R19): the pretype of a SUNDIALS Krylov solver, recognized by its solve routine so that a
  wrapper reporting the same ID is never read as one (PCG: LEFT unless NONE); LEFT for any other solver
```

*Listing 8.4. The FIRKLS protocol. The stacked operator is the block product*
$\mathrm{vec}(MV - h J V A^T)$ *for* $V = [v_1, \dots, v_s]$ *[2000VL], which a future multivector
layout could exploit.*

| Solve                 | Operator and solver                                                                               | Stopping test of the solver                                                                                      | `delta` passed                                                                                                                                                            | On non-convergence                                            |
| --------------------- | ------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------- |
| Newton, $s \ge 2$     | $I_s \otimes M - hA \otimes J$; internal SPFGMR, right preconditioned, scaling $(W, W)$           | $\lVert W \circ r \rVert_2 = \sqrt{sN} \lVert r \rVert_{WRMS}$ `[fact: src/sunlinsol/spfgmr/sunlinsol_spfgmr.c]` | $\epsilon_{stk} \epsilon_{\text{nls}} \sqrt{s}\; \mathrm{nrmfac}$, times $\frac{\lVert b \rVert}{\lVert x_0 \rVert}$ with a mass matrix                                   | recoverable; `SUNLS_RES_REDUCED` accepted iff $k = 1$         |
| Newton, $s = 1$       | $M - \gamma J$; the user's solver, scaling $(w, w)$                                               | the solver's                                                                                                     | $\epsilon_L \epsilon_{\text{nls}}\; \mathrm{nrmfac}$ (CVLS); with a mass matrix and an iterative solver, $\sigma$-scaled system or $\frac{\delta}{\sigma}$ (Section 5.10) | recoverable                                                   |
| error-estimate filter | $M - \gamma_p J$ ($\gamma_p = \gamma$ under the RADAU5 policy); the user's solver                 | the solver's                                                                                                     | $\epsilon_L \min(0.1, \lVert x_0 \rVert_{WRMS})\; \mathrm{nrmfac}$ (Table 5.7), converted as for $s = 1$; no small-residual shortcut                                      | recoverable, then $h \times 0.25$                             |
| preconditioner block  | $M - \gamma_p J$; a direct user solver solves exactly, an iterative one applies its `psolve` once | none; the `psolve` receives the tolerance as `delta`                                                             | $\epsilon_L \lVert r_i \rVert_{WRMS}\; \mathrm{nrmfac}$                                                                                                                   | accepted and counted                                          |
| mass solve            | $M$; the user's mass solver                                                                       | the solver's                                                                                                     | $\epsilon_{L,M} \lVert b \rVert_{WRMS}\; \mathrm{nrmfac}$, relative (Section 4.8)                                                                                         | `SUNLS_RES_REDUCED` accepted; otherwise `FIRK_MASSSOLVE_FAIL` |
| units conversion      | $M$; the user's mass solver, per stage block for $s \ge 2$                                        | the solver's                                                                                                     | $\epsilon_{L,M} \lVert b \rVert_{WRMS}\; \mathrm{nrmfac}$, $b$ the whole right-hand side of the enclosing solve                                                           | failure of the enclosing solve, with its sign                 |

*Table 8.6. Tolerance conversions. The SUNDIALS Krylov solvers test the 2-norm of the scaled
residual, so the WRMS tolerances of Table 5.7 are multiplied by* $\mathrm{nrmfac} = \sqrt{N}$, *as
CVLS and ARKLS do `[fact: src/cvode/cvode_ls.c, src/arkode/arkode_ls.c]`; `SetLSNormFactor`
overrides it with ARKLS semantics.* $N$ *comes from `N_VGetLength`, never from `N_VGetLocalLength`
(H10). Direct solvers ignore `delta`.*

Further rules of the interface:

- Iterative block solvers: only the user's preconditioner is applied inside `firkLsPSolve`, once per
  block and per FGMRES iteration, with `delta` $= \epsilon_L \lVert r_i \rVert$; the user's Krylov
  iteration never runs inside the stacked solve and serves the $s = 1$ Newton solve and the estimate
  (H37). Without a user preconditioner the stacked FGMRES is unpreconditioned. The $Jv$ products of
  the stacked operator use the saved matrix, the user's `jtimes`, or difference quotients [1989BH].
- The saved matrix is applied with `SUNMatMatvec`, after `SUNMatMatvecSetup` once per Jacobian
  evaluation when the matrix provides it `[fact: include/sundials/sundials_matrix.h]`. A
  matrix-based block solver uses a `SUNMatrix` that implements `SUNMatMatvec`, which the dense,
  band, sparse and GPU matrices of the suite do; for a matrix without it the stacked operator takes
  $Jv$ from an attached `jtimes`, else from difference quotients with one warning at initialization
  (R7). With `SetLinSysFn` no Jacobian is saved, so $Jv$ comes from `jtimes` or from difference
  quotients; setting a user `linsys` drops a saved Jacobian, and setting `NULL` saves one again from
  the next evaluation; the guide says so.
- `jok` and `jcur`: the user's `linsys` and `psetup` receive `jok = !jbad`, and whatever they report
  in `jcur` is recorded as the currency of the Jacobian (H29).
- Difference-quotient dense and band Jacobians follow CVLS: increments
  $\max\Bigl(\sqrt u \lvert y_j \rvert, \frac{\mathrm{minInc}}{w_j}\Bigr)$ with $\mathrm{minInc} = 1000 \lvert h \rvert u N \lVert f \rVert_{WRMS}$ or 1, and the perturbed component is restored even when the
  RHS reports a recoverable failure (H2). They are the only routines that read raw vector data.
- Mass matrix: $J$ and $M$ are both matrix-based or both matrix-free (ARKLS rule);
  `time_dep = SUNTRUE` is rejected with `FIRKLS_ILL_INPUT` until M9; the constant $M$ is set up once
  per problem, at the first step after `Init` or `ReInit` or after a mass routine is attached, and
  reused, also after `Reset` (H28). Mass products appear in the residual ($s$ per Newton iteration),
  in the stacked operator ($s$ per FGMRES iteration), and once in the estimate. Mass solves appear
  in the initial step, in rootfinding, and in the units conversion of Section 5.10: $s$ per Newton
  iteration for $s \ge 2$, and one per iterative block solve for the $s = 1$ Newton solve and the
  estimate.
- Linear-solution scaling by $\frac{2}{1 + \frac{\gamma}{\gamma_p}}$ applies to matrix-based block
  solvers (direct or matrix-iterative: they solve with the matrix of the last setup) in the $s = 1$
  Newton solve only, as in CVODE; never inside the preconditioner. It is on by default and
  `SetLinearSolutionScaling` rejects a matrix-free solver with `FIRKLS_ILL_INPUT`. The estimate
  solve uses the current factorization without correction: under the RADAU5 policy
  $\gamma_p = \gamma$ whenever the estimate is formed, because the factorization is redone when $h$
  changes; under the PERIODIC policy $\gamma_p$ may differ from $\gamma$ by up to 30%, and the
  filter is then $(M - \gamma_p J)^{-1}$, which changes neither the order of the estimate nor its
  stiff-limit bound.

### 8.7 The error estimate

```
firkErrEstimate(&dsm):
  v <- sum_i e_i Z_i; if M attached: v <- M v
  b <- gamma * fn + v
  flag <- firkLsSolveBlock(b, err, ESTIMATE, -)             # Table 8.6; no shortcut
  if flag > 0: return +1;  if flag < 0: return FIRK_LSOLVE_FAIL
  dsm <- max(||err||_w, dsm_floor)
  if dsm >= 1 and refilter and (first or attempt_rejected):   # attempt_rejected: any rejected attempt of this step
    ytmp <- yn + err; flag <- firkRhs(tn, ytmp, ftmp); nfe += 1
    if flag < 0: return FIRK_RHSFUNC_FAIL
    if flag == 0: b <- gamma * ftmp + v; solve into a temporary
                  if the solve fails unrecoverably: return FIRK_LSOLVE_FAIL
                  if it succeeds: err <- the temporary; dsm <- max(||err||_w, dsm_floor)
    a recoverable failure of the RHS or of the solve keeps the first err and dsm (H20)
  log error-estimate (dsm, refiltered)
  keep err for GetEstLocalErrors
  return 0
```

*Listing 8.5. The filtered estimate with RADAU5's refilter (Section 4.6). A recoverable failure of
the first solve fails the attempt with the factor `etacf` (Listing 8.2).*

### 8.8 Step control

```
firkAdapt(dsm, k_nls, accepted) -> eta, keep_setup:          # once per attempt, before UpdateH (H36)
  dsm <- max(dsm, dsm_floor)                                 # pow(0, negative) guard (H13)
  kappa <- safety * (newton_safety ? min(1, (1 + 2m) / (k_nls + 2m)) : 1)
  h_I <- kappa * |h| * dsm^(-1/(s+1))                         # RADAU5's I-controller proposal with FAC
  if mode == BEGINNER and not accepted: |hnew| <- h_I         # RADAU5 computes FACGUS on accepted steps only
  else:
    SUNAdaptController_EstimateStep(C, h, p = s, dsm, &hnew); failure -> FIRK_CONTROLLER_ERR
    if mode == BEGINNER: |hnew| <- min(safety * |hnew|, h_I)  # RADAU5: SAFE on FACGUS; QUOT = MAX(QUOT, FACGUS)
    else if accepted:    |hnew| <- kappa * |hnew|             # EXPERT: no I-controller minimum
    else:                |hnew| <- min(kappa * |hnew|, h_I)   # a rejected step always shrinks
  cap <- etamax    # etamx1 on the first step, growth after an accepted step, 1 after a rejection (RADAU5's REJECT)
  |hnew| <- min(|hnew|, cap * |h|);  |hnew| <- max(|hnew|, etamin * |h|)
  keep_setup <- FALSE
  if accepted and nls_theta <= theta_jac and lbound * |h| < |hnew| < ubound * |h| (up to u): hnew <- h; keep_setup <- TRUE
  eta <- hnew / h; eta <- max(eta, hmin_eff / |h|); eta <- eta / max(1, |h| * hmax_inv * eta)   # hmin, hmax last (R2, R18)
  log new-step-before-bounds, new-step-after-max-min-bounds, new-step-eta

firkCheckTemporalError(&nflag, &nef, dsm):
  (eta, keep_setup) <- firkAdapt(dsm, nls_iters, dsm <= 1)    # eta for firkCompleteStep; keep_setup records the dead band for the tests (the policies see gamma == gammap)
  if dsm <= 1: return FIRK_SUCCESS
  nef += 1; netf += 1; nflag <- PREV_ERR_FAIL
  if nef == maxnef or |h| <= hmin_eff (1 + u): return FIRK_ERR_FAILURE
  if first: eta <- min(eta, 0.1)                              # RADAU5: first step rejected
  else if nef >= small_nef: eta <- min(eta, etamxf)
  eta <- max(eta, hmin_eff / |h|); etamax <- 1                 # no growth for the rest of the step
  h <- eta * h; return TRY_AGAIN                               # gamma changed: the policy will set up
```

*Listing 8.6. Step control (Section 4.7). `EstimateStep` is called at most once per attempt, in
`firkCheckTemporalError`, before the history update in `firkCompleteStep`, as ARKODE does
`[fact: src/arkode/arkode_adapt.c, src/arkode/arkode.c]`; a second call after `UpdateH` would see the
current step as its own history and reduce `ImpGus` to the I-controller (H36).*

### 8.9 Dense output, predictor and order change

```
FIRKodeGetDky(mem, t, k, dky):
  require dky != NULL (FIRK_BAD_DKY), 0 <= k <= s_hist (FIRK_BAD_K)
  if not hist.valid (no step since Init, ReInit or Reset): k == 0 and t == tn up to 100u(|tn| + |h|) -> dky <- yn; else FIRK_BAD_T
  require t within [tn - hu, tn] up to 100u(|tn| + |hu|) (FIRK_BAD_T)
  theta <- (t - hist.t) / hist.h
  coefficients c_j <- L_j^(k)(theta) / hist.h^k - [k == 0] [j == s]      # Horner on hist.tbl->P
  dky <- [k == 0] yn + sum_j c_j hist.Z_j

firkPredict():
  if not hist.valid or hu == 0 or predictor == FIRK_PREDICT_TRIVIAL: Z <- 0; return
  for i: theta_i <- 1 + c_i h / hist.h; Z_i <- sum_j (L_j^hist(theta_i) - [j == s_hist]) hist.Z_j    # valid across an order change

firkOrderSelect() [v0.3; when a setup is due, or after 10 kept steps]:
  if s + 2 <= smax and knls_last > 1 and theta_filt <= 0.002 and 0.8 < eta_last < 1.2 and klin_avg <= 2s and steps_since_change >= 10: change to s + 2
  else if s - 2 >= smin and (theta_filt >= 0.8 or flag_newton_fail or flag_unexpected_reject): change to s - 2

firkChangeOrder(s_new):
  tbl <- FIRKodeTable_RadauIIA(s_new)                        # a copy of the built-in table, freed by FIRKodeFree
  allocate Z, F, G, delta, Jv for s_new; rebuild W; destroy and recreate SPFGMR with maxl = max(5, 3 s_new), restarts, GS type, ATimes, PSolve, Initialize
  on any failure: restore the old objects; return FIRK_ORDER_CHANGE_FAIL (H33)
  need_setup <- TRUE (gamma0 changed, J kept); SUNAdaptController_Reset(C); h unchanged; steps_since_change <- 0; nchanges += 1
  log order-change; hist keeps the old table and increments until the next completion
```

*Listing 8.7. Dense output, predictor and the order change of v0.3.*

### 8.10 Rootfinding

`firkode_root.c` adapts `src/arkode/arkode_root.c` with FIRKODE names: `firkRootCheck1` (initial
time), `firkRootCheck2` (the last returned point), `firkRootCheck3` (the new step), `firkRootfind`
(Illinois iteration on `FIRKodeGetDky`), `ttol = 100 u (\lvert t \rvert + \lvert h \rvert)`, root
directions, the inactive-root warning and `mxgnull` (D13). The public functions are `RootInit`,
`SetRootDirection`, `SetNoInactiveRootWarn`, `GetRootInfo` and `GetNumGEvals`. `firkGetYdot`, which
needs one mass solve when $M \ne I$, to the relative tolerance of Section 4.8, serves only the
initial-step heuristic and this module.

### 8.11 Error handling and return codes

| Range                               | Codes                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| ----------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| success and warnings                | `FIRK_SUCCESS` 0, `FIRK_TSTOP_RETURN` 1, `FIRK_ROOT_RETURN` 2, `FIRK_WARNING` 99                                                                                                                                                                                                                                                                                                                                                                                                                              |
| shared with CVODE (same meaning)    | `FIRK_TOO_MUCH_WORK` −1, `FIRK_TOO_MUCH_ACC` −2, `FIRK_ERR_FAILURE` −3, `FIRK_CONV_FAILURE` −4, `FIRK_LINIT_FAIL` −5, `FIRK_LSETUP_FAIL` −6, `FIRK_LSOLVE_FAIL` −7, `FIRK_RHSFUNC_FAIL` −8, `FIRK_FIRST_RHSFUNC_ERR` −9, `FIRK_REPTD_RHSFUNC_ERR` −10, `FIRK_UNREC_RHSFUNC_ERR` −11, `FIRK_RTFUNC_FAIL` −12                                                                                                                                                                                                   |
| unused                              | −13 to −19: CVODE's nonlinear-solver codes and ARKODE's mass codes collide here `[fact: include/cvode/cvode.h, include/arkode/arkode.h]`                                                                                                                                                                                                                                                                                                                                                                      |
| memory and input, shared with CVODE | `FIRK_MEM_FAIL` −20, `FIRK_MEM_NULL` −21, `FIRK_ILL_INPUT` −22, `FIRK_NO_MALLOC` −23, `FIRK_BAD_K` −24, `FIRK_BAD_T` −25, `FIRK_BAD_DKY` −26, `FIRK_TOO_CLOSE` −27, `FIRK_VECTOROP_ERR` −28                                                                                                                                                                                                                                                                                                                   |
| FIRKODE only                        | `FIRK_MASSINIT_FAIL` −14, `FIRK_MASSSETUP_FAIL` −15, `FIRK_MASSSOLVE_FAIL` −16, `FIRK_MASSMULT_FAIL` −18 (the ARKODE values), `FIRK_TABLE_FAIL` −44, `FIRK_CONTROLLER_ERR` −50 (the values of `ARK_INVALID_TABLE` and `ARK_CONTROLLER_ERR`), `FIRK_UNRECOGNIZED_ERR` −99                                                                                                                                                                                                                                      |
| FIRKLS (ARKLS values)               | `FIRKLS_SUCCESS` 0, `FIRKLS_MEM_NULL` −1, `FIRKLS_LMEM_NULL` −2, `FIRKLS_ILL_INPUT` −3, `FIRKLS_MEM_FAIL` −4, −5 unused (ARKLS's preconditioner-module code), `FIRKLS_MASSMEM_NULL` −6, `FIRKLS_JACFUNC_UNRECVR` −7, `FIRKLS_JACFUNC_RECVR` −8, `FIRKLS_MASSFUNC_UNRECVR` −9, `FIRKLS_MASSFUNC_RECVR` −10, `FIRKLS_SUNMAT_FAIL` −11, `FIRKLS_SUNLS_FAIL` −12                                                                                                                                                                           |
| private (`firkode_impl.h`)          | step flags `PREDICT_AGAIN` 3, `TRY_AGAIN` 5, `FIRST_CALL` 6, `PREV_CONV_FAIL` 7, `PREV_ERR_FAIL` 8; Newton statuses `FIRK_NLS_SUCCESS`, `FIRK_NLS_MAXITER`, `FIRK_NLS_DIVERGED`, `FIRK_NLS_PREDICTED_FAIL`, `FIRK_NLS_RHS_RECVR`, `FIRK_NLS_LSOLVE_RECVR`, `FIRK_NLS_LSETUP_RECVR` (positive); from v0.4 the hook return values `FIRK_NLS_CONTINUE` 901 and `FIRK_NLS_CONV_RECVR` 902 are public, mirroring `SUN_NLS_CONTINUE` and `SUN_NLS_CONV_RECVR` `[fact: include/sundials/sundials_nonlinearsolver.h]` |

*Table 8.2. Return codes (D10). `FIRKodeGetReturnFlagName` and `FIRKodeGetLinReturnFlagName` map
every value to its name; the FIRKLS setters of EXPERT options return the integrator code
`FIRK_ILL_INPUT` in BEGINNER mode, which `FIRKodeGetLinReturnFlagName` names, as it names every
negative integrator code that is not a FIRKLS code.*

Errors are reported through `firkProcessError(firk, code, __LINE__, __func__, __FILE__, fmt, ...)`,
which forwards to the context's error handlers as `cvProcessError` does; warnings go to the logger at
the warning level. `SUNCheckCall` and `SUNAssert` appear only in `firkode_sunstepper.c` (M8). Every
example program returns a nonzero exit code when the solver fails (H34).

### 8.12 Logging and profiling

| Level       | Label                                                                     | Content                                                                                                                                                                                                                                                |
| ----------- | ------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| info        | `begin-step-attempt`                                                      | `step = %li, tn = %g, h = %g, s = %i`                                                                                                                                                                                                                  |
| info        | `end-step-attempt`                                                        | `status = success, dsm = %g`, or `status = failed error test, dsm = %g`, or `status = failed solve, nflag = %i`                                                                                                                                        |
| info        | `begin-nonlinear-solve`, `end-nonlinear-solve`                            | `tol = %g`; `status = success, iters = %li` or `status = failed, retval = %i, iters = %li, hfac = %g`                                                                                                                                                  |
| info        | `begin-iterations-list`, `end-iterations-list`                            | around each Newton iteration, as `SUNNonlinSol_Newton` logs them: `status = continue`, `success` or `failure, retval = %i`                                                                                                                             |
| info        | `begin-linear-solve`, `end-linear-solve`                                  | per Newton or estimate solve: `purpose = newton or estimate, [stacked = 1,] iterative = %i, b-norm = %g, x-norm = %g, b-tol = %g, res-tol = %g`; `status = %s, retval = %i, iters = %li, p-solves = %li` (stacked: `block-solves = %li, resnorm = %g`) |
| info        | `begin-estimate-list`, `end-estimate-list`                                | around each estimate solve, of which a refiltered attempt has two                                                                                                                                                                                      |
| info        | `begin-mass-linear-solve-list`, `end-mass-linear-solve-list`              | `iterative = %i, b-tol = %g, res-tol = %g`; `status = %s, retval = %i, last-flag = %li, iters = %li`; a list, since one region can hold several mass solves                                                                                            |
| info        | `order-change` (v0.3)                                                     | `s-old = %i, s-new = %i, reason = raise or lower`                                                                                                                                                                                                      |
| debug       | `nonlinear-iterate`                                                       | `cur-iter = %i, total-iters = %li, update-norm = %g, theta = %g, eta = %g`                                                                                                                                                                             |
| debug       | `newton-exit`, `lsetup`, `error-estimate`                                 | `nls-status = %i, nls-hfac = %g`; `jbad = %i, jcur = %i, gamma = %g, retval = %i`; `est-norm = %g, refiltered = %i`                                                                                                                                        |
| debug       | `new-step-before-bounds`, `new-step-after-max-min-bounds`, `new-step-eta` | as ARKODE (`h_acc_controller`, `h_acc`, `eta`), with `h_I` before the bounds                                                                                                                                                                           |
| extra debug | `block-solve`, `predictor`                                                | one summary per stacked preconditioner application (`blocks = %i, accepted-failures = %li`); the predicted stages as `Z_i(:)` vectors                                                                                                                  |

*Table 8.4. Log records. The ARKODE labels are kept so that `suntools/src/suntools/logs.py` parses
FIRKODE logs; the parser merges the records of a region into one dictionary, so keys are unique
within a region and a region that can repeat within its parent is a list. Every exit path emits the
matching `end-*` record (H31): a thin wrapper closes the region a solve opened. `linear-iterate` records come
from SPFGMR. Logging compiles to nothing below `SUNDIALS_LOGGING_LEVEL` 3 `[fact:
cmake/SundialsBuildOptionsPre.cmake]`.*

| Region                                                                                                                                                                                                                                                                                                                       | Where                                                                                                                                            |
| ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------ |
| `FIRKodeInit`, `FIRKodeReInit`, `FIRKodeReset`, `FIRKodeEvolve`, `FIRKodeGetDky` (function regions, as CVODE's), `firkStep`, `firkNewton`, `firkLsSetup`, `firkLsSolveStacked`, `firkLsSolveBlock`, `firkRhs`, `firkLsMassMult`, `firkLsMassSolve`, `firkErrorEstimate`, `firkRootfind`, `firkHin`, `firkChangeOrder` (v0.3) | `SUNDIALS_MARK_FUNCTION_BEGIN` and `SUNDIALS_MARK_FUNCTION_END` with `firk_mem->sunctx->profiler`; active only under `SUNDIALS_ENABLE_PROFILING` |

*Table 8.5. Profiling regions.*

### 8.13 Files, build and wiring

| File                                                                     | Contents                                                                               | Milestone  |
| ------------------------------------------------------------------------ | -------------------------------------------------------------------------------------- | ---------- |
| `include/firkode/firkode.h`                                              | return codes, constants, callback typedefs, integrator prototypes                      | M0         |
| `include/firkode/firkode_ls.h`                                           | FIRKLS and mass-matrix prototypes and typedefs                                         | M0         |
| `include/firkode/firkode_tables.h`                                       | `FIRKodeTableMem`, `FIRKodeTable_*`                                                    | M0         |
| `src/firkode/firkode_tables.c` (generated region)                        | literal tables, generated between markers; `FIRKodeTable_*` written by hand            | M0         |
| `src/firkode/firkode_impl.h`, `firkode_ls_impl.h`, `firkode_root_impl.h` | the memory objects of Table 8.1, private flags and defaults                            | M0, M1, M2 |
| `src/firkode/firkode.c`                                                  | create, init, free, evolve loop, initial step, failure handling, `GetDky` entry        | M1, M2     |
| `src/firkode/firkode_step.c`                                             | step attempt, Newton-flag handling, estimate, complete step                            | M1, M2     |
| `src/firkode/firkode_nls.c`                                              | the Newton driver                                                                      | M1         |
| `src/firkode/firkode_dense.c`                                            | collocation polynomial: predictor, dense output                                        | M1         |
| `src/firkode/firkode_ls.c`                                               | FIRKLS: block and stacked solves, operator, preconditioner, difference quotients, mass | M1         |
| `src/firkode/firkode_stack.c`                                            | stages' stack: kind, creation, weight view, stage accessor                             | M1         |
| `src/firkode/firkode_adapt.c`                                            | controller glue, reuse policies                                                        | M2         |
| `src/firkode/firkode_root.c`                                             | rootfinding                                                                            | M2         |
| `src/firkode/firkode_io.c`, `firkode_cli.c`                              | setters, getters, statistics, `WriteParameters`, `PrintAllStats`; `SetOptions`         | M1, M2     |
| `src/firkode/firkode_order.c` (+ `_impl.h`)                              | order selection                                                                        | M5         |
| `src/firkode/firkode_sunstepper.c`                                       | `SUNStepper` glue                                                                      | M8         |
| `scripts/firkode_radau_tables.py`                                        | table generator: `--emit-c`, `--check`, `--table1`, `--emit-lean`                      | M0         |
| `src/firkode/CMakeLists.txt`, `README.md`, `LICENSE`                     | library definition and package metadata                                                | M0         |

*Table 8.9. Planned file layout.*

Wiring that the first milestone performs, following the developer checklist
`[fact: doc/superbuild/source/developers/getting_started/Checklist.rst]`:

- CMake: `SUNDIALS_ENABLE_FIRKODE` in `cmake/SundialsBuildOptionsPre.cmake`, with the ManyVector and
  MPIManyVector options forced on in `cmake/SundialsBuildOptionsPost.cmake` when FIRKODE is on
  (Section 8.1); `firkodelib_VERSION` and `SOVERSION` in the root `CMakeLists.txt`;
  `sundials_add_library(sundials_firkode SOURCES ... HEADERS ... INCLUDE_SUBDIR firkode
  LINK_LIBRARIES PUBLIC sundials_core OBJECT_LIBRARIES ... OUTPUT_NAME sundials_firkode VERSION ...
  SOVERSION ...)`; subdirectory entries in `src/CMakeLists.txt`, `examples/CMakeLists.txt` and
  `test/unit_tests/CMakeLists.txt`.
- Examples: `examples/firkode/C_serial/CMakeLists.txt` with `sundials_add_executable`,
  `sundials_add_test(... ANSWER_DIR ... ANSWER_FILE ... EXAMPLE_TYPE ...)` and
  `sundials_install_examples`; `EXAMPLE_TYPE develop` for tests that need
  `SUNDIALS_TEST_ENABLE_DEV_TESTS`, `exclude-single` for tests that single precision cannot pass.
- Unit tests: `test/unit_tests/firkode/C_serial` (exit-code tests, linking
  `sundials_firkode_obj` for internal symbols) and `CXX_serial` (diffed developer tests).
- Documentation: `doc/firkode/guide`, registration in `doc/superbuild/source/index.rst` and
  `conf.py`, the `*_package_links.rst` files, `doc/shared/sundials_vars.py`,
  `doc/shared/sundials/Install.rst`, `doc/Makefile`, `CITATIONS.md`, `README.md`.
- Release: `scripts/firkode.sh`, `scripts/tarscript.sh`, `scripts/updateVersion.sh` (independent
  version), `scripts/startReleaseCycle.sh`, a `firkode` variant in
  `scripts/spack/packages/sundials/package.py`, `CHANGELOG.md` and `doc/shared/RecentChanges.rst`.
- CI: no new workflow; the generator's `--check` runs in the format-check job, and `ruff` formats
  the Python.

## 9. Public API

Function names omit the `FIRKode` prefix except where stated. The mode column says whether BEGINNER
mode accepts the call (B) or rejects it with `FIRK_ILL_INPUT` (E). The version column names the
release in which the function first does something; functions of later versions are absent from the
headers of earlier ones, except `SetNumStages` and `SetOrder`, which exist from v0.1.0 and accept
only $s = 3$ until v0.2.

### 9.1 Functions by group

| Group                          | Functions                                                                                                                                                                                                                                                                                                                                                                                                                                                     | Version                     | Mode                 |
| ------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------- | -------------------- |
| lifecycle                      | `Create(sunctx)`, `Init(mem, f, t0, y0)`, `ReInit(mem, t0, y0)` (zeroes the counters), `Reset(mem, tR, yR)` (keeps them), `Free(&mem)`, `SetMode(mem, FIRK_BEGINNER or FIRK_EXPERT)`, `GetMode`                                                                                                                                                                                                                                                               | v0.1.0                      | B                    |
| tolerances                     | `SStolerances(mem, rtol, atol)`, `SVtolerances(mem, rtol, atol_vec)`, `WFtolerances(mem, efun)`                                                                                                                                                                                                                                                                                                                                                               | v0.1.0                      | B                    |
| integration                    | `Evolve(mem, tout, yout, &tret, itask)` with `FIRK_NORMAL` or `FIRK_ONE_STEP`; `GetDky(mem, t, k, dky)`, $0 \le k \le s$                                                                                                                                                                                                                                                                                                                                      | v0.1.0                      | B                    |
| integrator options             | `SetUserData`, `SetFixedStep`, `SetStopTime`, `SetInterpolateStopTime`, `ClearStopTime`, `SetMaxNumSteps`, `SetMaxHnilWarns`, `SetInitStep`, `SetMinStep`, `SetMaxStep`                                                                                                                                                                                                                                                                                       | v0.1.0                      | B                    |
|                                | `SetMaxErrTestFails`, `SetMaxConvFails`, `SetPredictorMethod(FIRK_PREDICT_EXTRAPOLATE or FIRK_PREDICT_TRIVIAL)`, `SetErrorRefilter(on)`                                                                                                                                                                                                                                                                                                                       | v0.1.0                      | E                    |
| step adaptivity (ARKODE names) | `SetSafetyFactor`, `SetErrorBias`, `SetMaxGrowth`, `SetMinReduction`, `SetFixedStepBounds`, `SetMaxFirstGrowth`, `SetMaxEFailGrowth`, `SetSmallNumEFails`, `SetMaxCFailGrowth`, `SetNewtonCountSafety(on)`, `SetAdaptController(C)` (of type `SUN_ADAPTCONTROLLER_H`), `SetAdaptControllerByName(name)` (the Soderlind-family names of ARKODE, without `ImExGus`)                                                                                             | v0.1.0                      | E                    |
| Newton                         | `SetMaxNonlinIters`, `SetNonlinConvCoef`, `SetNonlinDivergenceRate`, `SetNonlinRateExponent`, `SetJacReuseConvRate` ($\theta_{jac}$, shared by the dead band and the RADAU5 policy), `SetJacBadConvRate`                                                                                                                                                                                                                                                      | v0.1.0                      | E                    |
| reuse                          | `SetReusePolicy(FIRK_POLICY_RADAU5 or FIRK_POLICY_PERIODIC)`, `SetLSetupFrequency`, `SetDeltaGammaMaxLSetup`; in FIRKLS: `SetJacEvalFrequency`, `SetDeltaGammaMaxBadJac`                                                                                                                                                                                                                                                                                      | v0.1.0                      | E                    |
| method                         | `SetNumStages(s)`, `SetOrder(q)` (the smallest $s$ with $2s - 1 \ge q$, the default for $q \le 0$; $q > 17$ is `FIRK_ILL_INPUT`; v0.1.0 accepts $q \le 0$, 4 and 5), `GetNumStages`, `GetCurrentTable(mem, &tbl_ptr)` (not owned by the caller)                                                                                                                                                                                                               | v0.1.0 ($s = 3$ only), v0.2 | setters E, getters B |
| order adaptivity               | `SetOrderRange(smin, smax)` (odd), `GetNumOrderChanges`, `GetNumStepsPerOrder(mem, nst_1d)`                                                                                                                                                                                                                                                                                                                                                                   | v0.3                        | setter E, getters B  |
| tables                         | `SetTable(mem, tbl)`                                                                                                                                                                                                                                                                                                                                                                                                                                          | v0.4                        | E                    |
| tables (`firkode_tables.h`)    | `FIRKodeTable_RadauIIA(s)` returns a caller-owned copy of a built-in table; `FIRKodeTable_Alloc(s)`, `_Copy(tbl)`, `_Free(tbl)`, `_Write(tbl, FILE*)`, `_CheckOrder(tbl, &q, &p, &stage_order, FILE*)`                                                                                                                                                                                                                                                                    | v0.1.0                      | both                 |
| FIRKLS                         | `SetLinearSolver(mem, LS, A)` (`A` implements `SUNMatMatvec` unless `jtimes` is attached, Section 8.6), `SetJacFn`, `SetLinSysFn`, `SetPreconditioner(mem, psetup, psolve)`, `SetJacTimes(mem, jtsetup, jtimes)`, `SetJacTimesRhsFn`, `SetEpsLin`, `SetLSNormFactor`, `SetLinearSolutionScaling`                                                                                                                                                              | v0.1.0                      | B                    |
| stage solver                   | `SetStageSolverMaxl`, `SetStageSolverMaxRestarts`, `SetStageSolverGSType`, `SetStageEpsLin`                                                                                                                                                                                                                                                                                                                                                                   | v0.1.0                      | E                    |
|                                | `SetStagePreconditioner(mem, psetup, psolve)`, `SetStageLinearSolver(mem, LS)`                                                                                                                                                                                                                                                                                                                                                                                | v0.4                        | E                    |
| mass matrix                    | `SetMassLinearSolver(mem, LS, M, time_dep)`, `SetMassFn`, `SetMassTimes(mem, mtsetup, mtimes, mtimes_data)`, `SetMassPreconditioner`, `SetMassEpsLin`, `SetMassLSNormFactor`                                                                                                                                                                                                                                                                                  | v0.1.0                      | B                    |
| rootfinding                    | `RootInit(mem, nrtfn, g)`, `SetRootDirection(mem, rootdir_1d)`, `SetNoInactiveRootWarn`, `GetRootInfo(mem, rootsfound_1d)`, `GetNumGEvals`                                                                                                                                                                                                                                                                                                                    | v0.1.0                      | B                    |
| hooks                          | `SetReusePolicyFn(mem, fn, data)`, `SetNlsConvTestFn`, `SetNlsNormFn`, `SetPredictorFn`, `SetDomEigEstimator(mem, DEE)`, `SetDomEigFrequency`, `SetSpectralRadiusFn`, `SetPrecShift`                                                                                                                                                                                                                                                                          | v0.4                        | E                    |
| options and output             | `SetOptions(mem, "firkode", NULL, argc, argv)`, `WriteParameters(mem, FILE*)`, `PrintAllStats(mem, FILE*, SUNOutputFormat)`, `GetReturnFlagName(flag)`, `GetLinReturnFlagName(flag)`                                                                                                                                                                                                                                                                          | v0.1.0                      | B                    |
| integrator statistics          | `GetNumSteps`, `GetNumStepAttempts`, `GetNumRhsEvals`, `GetNumLinSolvSetups`, `GetNumErrTestFails`, `GetNumStepSolveFails`, `GetNumNonlinSolvIters`, `GetNumNonlinSolvConvFails`, `GetNonlinSolvStats`, `GetLastNonlinConvRate`, `GetActualInitStep`, `GetLastStep`, `GetCurrentStep`, `GetCurrentTime`, `GetCurrentState`, `GetCurrentGamma`, `GetTolScaleFactor`, `GetErrWeights`, `GetEstLocalErrors`, `GetIntegratorStats`, `GetUserData`                 | v0.1.0                      | B                    |
| FIRKLS statistics              | `GetJac`, `GetJacTime`, `GetJacNumSteps`, `GetNumJacEvals`, `GetNumPrecEvals`, `GetNumPrecSolves`, `GetNumLinIters`, `GetNumLinConvFails`, `GetNumJTSetupEvals`, `GetNumJtimesEvals`, `GetNumLinRhsEvals`, `GetLastLinFlag`, `GetNumStageLinIters`, `GetNumStageLinConvFails`, `GetNumBlockSolves`                                                                                                                                                            | v0.1.0                      | B                    |
| mass statistics                | `GetCurrentMassMatrix`, `GetNumMassSetups`, `GetNumMassMultSetups`, `GetNumMassMult`, `GetNumMassSolves`, `GetNumMassPrecEvals`, `GetNumMassPrecSolves`, `GetNumMassIters`, `GetNumMassConvFails`, `GetLastMassFlag`                                                                                                                                                                                                                                          | v0.1.0                      | B                    |
| SUNStepper                     | `CreateSUNStepper(mem, &stepper)`, `SetStepDirection`                                                                                                                                                                                                                                                                                                                                                                                                         | optional M8                 | B                    |

*Table 9.1. Public functions by group.*

#### Constants

- `FIRK_NORMAL = 1`
- `FIRK_ONE_STEP = 2`
- `FIRK_BEGINNER = 0`
- `FIRK_EXPERT = 1`
- `FIRK_PREDICT_EXTRAPOLATE = 0`
- `FIRK_PREDICT_TRIVIAL = 1`
- `FIRK_POLICY_RADAU5 = 0`
- `FIRK_POLICY_PERIODIC = 1`
- `FIRK_MAX_STAGES = 9`
- `FIRK_DEFAULT_STAGES = 3`

#### Mode rules

The default is `FIRK_BEGINNER`; `SetMode` is accepted until the first step attempt and again after
`ReInit`, and a return to BEGINNER restores every EXPERT default; a rejected call returns
`FIRK_ILL_INPUT` and reports "`<function>` requires `FIRKodeSetMode(mem, FIRK_EXPERT)`";
`WriteParameters` prints the mode first. A FIRKLS setter of an EXPERT option checks the memory
(`FIRKLS_MEM_NULL`), then the mode (`FIRK_ILL_INPUT`), then the linear solver memory
(`FIRKLS_LMEM_NULL`) and its argument (`FIRKLS_ILL_INPUT`).

### 9.2 User-supplied functions

| Typedef                        | Signature                                                                                                                                  | Mirrors                     |
| ------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------ | --------------------------- |
| `FIRKRhsFn`                    | `int (t, N_Vector y, N_Vector ydot, void* user_data)`                                                                                      | `CVRhsFn`, `ARKRhsFn`       |
| `FIRKRootFn`                   | `int (t, N_Vector y, sunrealtype* gout_1d, void* user_data)`                                                                               | `ARKRootFn`                 |
| `FIRKEwtFn`                    | `int (N_Vector y, N_Vector ewt, void* user_data)`                                                                                          | `CVEwtFn`                   |
| `FIRKLsJacFn`                  | `int (t, y, fy, SUNMatrix Jac, void* user_data, N_Vector tmp1, tmp2, tmp3)`                                                                | `ARKLsJacFn`                |
| `FIRKLsLinSysFn`               | `int (t, y, fy, SUNMatrix A, SUNMatrix M, sunbooleantype jok, sunbooleantype* jcur, sunrealtype gamma, void* user_data, tmp1, tmp2, tmp3)` | `ARKLsLinSysFn`             |
| `FIRKLsPrecSetupFn`            | `int (t, y, fy, sunbooleantype jok, sunbooleantype* jcurPtr, sunrealtype gamma, void* user_data)`                                          | `ARKLsPrecSetupFn`          |
| `FIRKLsPrecSolveFn`            | `int (t, y, fy, N_Vector r, N_Vector z, sunrealtype gamma, sunrealtype delta, int lr, void* user_data)`                                    | `ARKLsPrecSolveFn`          |
| `FIRKLsJacTimesSetupFn`        | `int (t, y, fy, void* user_data)`                                                                                                          | `ARKLsJacTimesSetupFn`      |
| `FIRKLsJacTimesVecFn`          | `int (N_Vector v, N_Vector Jv, t, y, fy, void* user_data, N_Vector tmp)`                                                                   | `ARKLsJacTimesVecFn`        |
| `FIRKLsMassFn`                 | `int (t, SUNMatrix M, void* user_data, tmp1, tmp2, tmp3)`                                                                                  | `ARKLsMassFn`               |
| `FIRKLsMassTimesSetupFn`       | `int (t, void* mtimes_data)`                                                                                                               | `ARKLsMassTimesSetupFn`     |
| `FIRKLsMassTimesVecFn`         | `int (N_Vector v, N_Vector Mv, t, void* mtimes_data)`                                                                                      | `ARKLsMassTimesVecFn`       |
| `FIRKLsMassPrecSetupFn`        | `int (t, void* user_data)`                                                                                                                 | `ARKLsMassPrecSetupFn`      |
| `FIRKLsMassPrecSolveFn`        | `int (t, N_Vector r, N_Vector z, sunrealtype delta, int lr, void* user_data)`                                                              | `ARKLsMassPrecSolveFn`      |

*Table 9.2. User-supplied functions `[fact: include/arkode/arkode_ls.h, include/cvode/cvode.h]`.
Every Jacobian-related callback is evaluated at one point per step,* $(t_n, y_n, f(t_n, y_n))$,
*with `gamma`* $= h\gamma_0$. *Return conventions are CVODE's: 0 success, positive recoverable,
negative unrecoverable.*

### 9.3 Hooks (v0.4)

| Typedef                                              | Signature                                                                                                                                                                                                                                                                                        | Replaces                                        |
| ---------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ----------------------------------------------- |
| `FIRKodeReusePolicyFn`                               | `int (void* firkode_mem, const FIRKodeStepRecord* rec, sunbooleantype* jbad_ptr, sunbooleantype* setup_ptr, void* data)`                                                                                                                                                                         | the built-in policies                           |
| `FIRKodeNlsConvTestFn`                               | `int (void* firkode_mem, int iter, N_Vector del, sunrealtype delnrm, sunrealtype tol, N_Vector ewt_stk, sunrealtype* hfac_ptr, void* data)` returning `FIRK_NLS_CONTINUE`, `FIRK_SUCCESS` or `FIRK_NLS_CONV_RECVR`                                                                               | the test block of Listing 8.3                   |
| `FIRKodeNlsNormFn`                                   | `int (N_Vector del, N_Vector w, sunrealtype* nrm_ptr, void* data)`                                                                                                                                                                                                                               | the stacked WRMS norm                           |
| `FIRKodePredictorFn`                                 | `int (void* firkode_mem, sunrealtype t, sunrealtype h, N_Vector* Z_1d, void* data)`                                                                                                                                                                                                              | the extrapolation predictor                     |
| `FIRKodeStagePrecSetupFn`, `FIRKodeStagePrecSolveFn` | `int (t, y, fy, sunbooleantype jok, sunbooleantype* jcurPtr, sunrealtype gamma, sunrealtype h, const FIRKodeTableMem* tbl, void* user_data)`; `int (t, y, fy, N_Vector r_stk, N_Vector z_stk, sunrealtype gamma, sunrealtype h, const FIRKodeTableMem* tbl, sunrealtype delta, void* user_data)` | the single-shift block preconditioner           |
| `FIRKodeSpectralRadiusFn`                            | `int (t, y, fy, sunrealtype* rho_ptr, void* user_data)`                                                                                                                                                                                                                                          | the dominant-eigenvalue estimator               |

*Table 9.3. Hook typedefs of v0.4. `FIRKodeStepRecord` is the public per-step record of Section
4.11* ($\theta$, $k_{\text{nls}}$, $k_{lin}$, *dsm,* $\eta$, $\frac{\Delta\gamma}{\gamma}$, *steps since
the last* $J$ *and setup, failure flags). The default path must be bit-identical when a hook is
unset.*

### 9.4 Options and bindability

- `SetOptions(mem, "firkode", NULL, argc, argv)` reads `firkode.<key> <value>` pairs with the key and
  the value as separate arguments, through the helpers of `src/sundials/sundials_cli.h` as
  `src/arkode/arkode_cli.c` does `[fact]`. There is one key per scalar setter, named after the setter
  in snake case: `firkode.max_nonlin_iters 7`, `firkode.reuse_policy radau5`,
  `firkode.scalar_tolerances 1e-6 1e-10`, `firkode.mode expert`. A key that maps to an EXPERT setter
  is rejected in BEGINNER mode like the setter; `firkode.mode` is applied before every other key of
  the same call, whatever its position in `argv`, so that `firkode.mode expert` unlocks the EXPERT
  keys beside it. File-based options, a key without its value and an unknown name return
  `FIRK_ILL_INPUT`; an unknown key warns and is skipped.
- Bindability rules, so that SWIG-Fortran and `sundials4py` can follow later without API changes:
  array parameters carry `_1d` or `_2d` suffixes and the other pointer parameters `_ptr`
  `[fact: doc/superbuild/source/developers/source_code/Naming.rst]`, outputs of callbacks included
  (`jcur_ptr`). The generator reads the array suffixes: it returns every other pointer to an integer
  or real as an output in a tuple
  `[fact: bindings/sundials4py/sundials4py-generate/generate.py:81–85]`.
  Exceptions: the opaque `void* firkode_mem` and `void* user_data`, `FILE*`, strings, and
  `void** firkode_mem` of `Free(&mem)`, as in every package. No `**` or `***` parameters in public
  functions, except `void**` in `Free` and `GetUserData`, pointers to handles returned without
  ownership (`N_Vector*`, `SUNMatrix*`, `FIRKodeTable*`) as in ARKODE
  `[fact: include/arkode/arkode.h]`, and `char* argv[]` in `SetOptions` as in `ARKodeSetOptions`;
  no variadic public functions; every callback carries a `void* user_data` or `mtimes_data` pointer;
  pointers returned without ownership carry the `// nb::rv_policy::reference` hint
  `[fact: include/nvector/nvector_manyvector.h]`.
- A C++ view `firkode.hpp` with a deleter and a `sundials::experimental::FIRKodeView` follows the
  other packages.

## 10. Implementation milestones

Each milestone lists its goal, deliverables, tests, the hazard rows of Table 11.3 it must close, the
principle that sizes it, and the criteria that close it. Every criterion is a test in the tree or a
CI job. Every milestone also ships the user-guide pages for its API and an entry in `CHANGELOG.md`
and `doc/shared/RecentChanges.rst`.

| Milestone | Version  | Goal                                                                 | Depends on | Principle |
| --------- | -------- | -------------------------------------------------------------------- | ---------- | --------- |
| M0        | v0.1.0   | tables, generator, package skeleton, build and release wiring        | –          | P1        |
| M1        | v0.1.0   | fixed-step $s = 3$ with the final stage solver                       | M0         | P4        |
| M2        | v0.1.0   | adaptive $s = 3$: estimate, control, reuse, failures, roots, options | M1         | P2        |
| M3        | v0.1.0   | hardening, CI matrix, answer files, MPI example                      | M2         | P1        |
| M4        | v0.2     | the family $1 \le s \le 9$ at fixed order, with measurements         | M3         | P4        |
| M5        | v0.3     | in situ order selection                                              | M4         | P4        |
| M6        | v0.4     | EXPERT stage-solver hooks and spectral information                   | M4         | P4        |
| M7        | optional | Lean certification of the generated literals                         | M0         | P1        |
| M8        | optional | `SUNStepper` glue                                                    | M2         | P3 (last) |
| M9        | optional | time-dependent $M$                                                   | M4         | –         |
| M10       | optional | singly implicit tables                                               | M6         | –         |
| M11       | optional | resize for adaptive meshes                                           | M2         | –         |
| M12       | optional | relaxation                                                           | M2         | –         |
| M13       | optional | discrete adjoint                                                     | M4         | –         |

*Table 10.1. Milestones. M7 runs in parallel with everything after M0; M8–M13 are not scheduled.*

### 10.1 M0: tables and skeleton

#### Goal

every coefficient the integrator will use exists as a checked literal, and the empty package
builds in every CI configuration.

#### Deliverables


- `scripts/firkode_radau_tables.py`: nodes by bisection of $P_s(2c-1) - P_{s-1}(2c-1)$ with `mpmath`
  at 50 digits; $A$ from the collocation conditions; $A^{-1}$, $d$, $e$, $\gamma_0$ (Section 4.5) and
  the dense-output coefficients $P$; `--emit-c` writes `src/firkode/firkode_tables.c` with 40-digit
  `SUN_RCONST` literals for $1 \le s \le 9$ and a header comment with the generator version;
  `--check` exits nonzero when the file is stale; `--table1` prints Table 4.1; `--emit-lean` writes
  the exact rational data of M7. Formatted with `ruff`.
- `include/firkode/firkode.h`, `firkode_ls.h`, `firkode_tables.h`: all codes, constants and
  typedefs of Chapter 9; prototypes of the v0.1.0 functions.
- `FIRKodeTable_RadauIIA`, `_Alloc`, `_Copy`, `_Free`, `_Write`, `_CheckOrder` (order conditions
  $B$, $C$, $D$ in double precision with tolerance $10^3 s^2 u$, returning $q$, $p$ and the stage
  order; $p$ is the largest $k$ for which the embedded quadrature of Section 4.6, the weights
  $b^T + e^T A$ on the nodes and $\gamma_0$ on the node 0, is exact for degree $k - 1$).
- `src/firkode/CMakeLists.txt`, `firkode_impl.h` skeleton, `README.md`, `LICENSE`; the CMake, docs
  and release wiring of Section 8.13; the documentation skeleton `doc/firkode/guide/source`.
- `test/unit_tests/firkode/C_serial/firk_test_tables.c`.

#### Done when


- the library builds in every CI configuration: gcc int32 and int64 with double (Debug and
  Release), single and extended precision, clang with logging levels 0–5 and profiling on and off,
  ASan, oneAPI, macOS, MSVC, Intel and MinGW on Windows, all with `CMAKE_COMPILE_WARNING_AS_ERROR`
  and `SUNDIALS_ENABLE_ALL_WARNINGS` `[fact: .github/workflows]` (H24, H25);
- `firk_test_tables` passes for $s \le 9$ in double and in extended precision, both 80-bit x87 and
  binary128 `long double`, and for $s \le 4$ in single (the `CheckOrder` result for $s \le 3$):
  nodes and $c_s = 1$, row sums, $B(2s-1)$, $C(s)$, $A^{-1} A = I$, $d = \mathbf e_s$,
  $e_s = (-1)^s \frac{\gamma_0}{s}$, $L_i(c_j) = \delta_{ij}$, $\det(A^{-1} - \gamma_0^{-1} I) = 0$ for odd
  $s$, RADAU5's constants for $s = 3$, and `CheckOrder` returning $(2s - 1, s, s)$ with negative
  controls (a perturbed table fails the right condition);
- `--check` is green in CI; the documentation builds.

### 10.2 M1: fixed-step RADAU5 with the final stage solver

#### Goal

a fixed-step $s = 3$ integrator whose stage solve is the final architecture (P4), so that
every later milestone changes control logic, never the solver.

#### Deliverables

`firkode.c` (`Create`, `Init`, `Free`, `Evolve` in both tasks with collocation output
at `tout`, `GetDky`, flag names, tolerances including `WFtolerances`, `SetFixedStep`, `SetStopTime`,
`SetUserData`), `firkode_step.c` without the estimate, `firkode_nls.c` (Listing 8.3),
`firkode_dense.c`, `firkode_ls.c` complete (Listing 8.4: block and stacked solves,
difference-quotient dense and band Jacobians, difference-quotient $Jv$, preconditioner, `linsys`,
constant mass matrix), `firkode_io.c` with the M1 setters, every counter and getter, the mode gate,
logging and profiling, the stage stack with its accessor and weight view, the example
`firk_analytic` with its answer file.

Tests (Chapter 11): `firk_test_dahlquist`, `firk_test_order` on Prothero–Robinson and
Kværnø–Prothero–Robinson in fixed-step mode, `firk_test_krylov`, `firk_test_mass`,
`firk_test_dense_output`, `firk_test_workspace`, `firk_test_stack_ops` (the trapping mock of H9),
`firk_test_mode`, `firk_test_firkls` (the recording mock of Listing 8.4), `firk_test_newton`,
`firk_test_ewt`. Hazard rows H2, H4, H6, H7, H9–H11, H15–H18, H24, H25, H27–H29, H31, H35, H37.

#### Done when


- Dahlquist: in EXPERT mode with $\epsilon_{\text{nls}} = 10^{-12}$, $\epsilon_{stk} = 1$ and unit
  error weights, so that the Newton and Krylov tolerances lie below the pass threshold, one fixed
  step reproduces $R(hL) y_0$ within $10^5 u (1 + \lVert y_{\mathrm{ref}} \rVert)$ for real $h\lambda$
  down to $-10^7$ and for complex rotation blocks; a linear problem takes at most 2 Newton iterations;
- fixed-step order 5 on Kværnø–Prothero–Robinson and non-stiff Prothero–Robinson, at least 3 on
  stiff Prothero–Robinson ($\lambda = -10^4$); last observed order $\ge$ expected $- 0.3$;
- the matrix-free path (unpreconditioned FGMRES with difference-quotient $Jv$) agrees with the dense
  path to $10^{-8}$; $M \dot y = f$ agrees with $\dot y = M^{-1} f$ with dense and PCG mass solvers,
  and the mass counters equal $s(1 + k_{lin})$ products and $s$ solves per Newton iteration;
- `GetDky` is exact for polynomial solutions of degree $\le s$ and fails for degree $s + 1$;
- the vector budget equals Table 8.7 (`firk_test_workspace`);
- ASan and UBSan clean; builds and passes in single, double and extended precision.

### 10.3 M2: adaptive RADAU5

#### Goal

RADAU5's adaptive behavior for $s = 3$ plus the SUNDIALS conveniences (P2: this is the
functional content of v0.1.0).

#### Deliverables

the estimate with filter and refilter (Listing 8.5), `firkode_adapt.c` (Listing 8.6,
the reuse policies, `SetAdaptController` and the ARKODE-named adaptivity setters), the initial-step
heuristic, failure handling (Table 8.3, recoverable RHS and callback failures, the failure limits),
`ReInit`, `Reset`, backward integration, `SetInterpolateStopTime`, `ClearStopTime`, `firkode_root.c`,
`firkode_cli.c`, `WriteParameters`, `PrintAllStats`, `GetIntegratorStats`, `GetEstLocalErrors`, the
examples `firk_robertson`, `firk_robertson_root`, `firk_diurnal_kry`, `firk_heat1D_fem_mass`, and the
RADAU5 comparison harness: a script that runs `radau5.f` offline on Robertson, van der Pol
($\mu = 1000$) and the Oregonator with transformed tolerances and records the step sequences, committed
with the reference files under `test/`, and the test `firk_test_radau5_steps` that reads them.

#### Tests

`firk_test_errest`, `firk_test_robertson`, `firk_test_vdp`, `firk_test_radau5_steps`,
`firk_test_tstop`, `firk_test_reset`, `firk_test_root`, `firk_test_rootfind`, `firk_test_recoverable`,
`firk_test_options`, `firk_test_logging`, `firk_test_adapt`, `firk_test_reuse`, `firk_test_userdata`,
`firk_test_hnil`. Hazard rows H1, H3, H5, H8, H13, H19–H22, H30, H32, H34, H36.

#### Done when


- the estimate order $p = 3$ is measured on a smooth problem (slope of $\mathrm{dsm}$ against $h$);
- on Robertson, van der Pol and the Oregonator the global error at the final time, against
  reference solutions computed with `radau5.f` at $\mathrm{rtol} = 10^{-14}$, is within a factor 2
  of `radau5.f`'s own error or below half the tolerance, at equal transformed tolerances for
  $\mathrm{rtol} \in \{10^{-4}, 10^{-6}, 10^{-8}\}$, with RADAU5's initial step and the final time as
  a stop time (Table 6.2); the floor is needed because both errors lie well below the tolerance,
  where the error is set by the path through the stiff transitions (Section 12, question 15); the
  estimate agrees within 50% between a dense and an iterative block solver;
- step counts on Robertson, van der Pol and the Oregonator are within 20% of `radau5.f` at equal
  transformed tolerances under the RADAU5 policy;
- stop-time semantics in adaptive and fixed-step mode, `ReInit` zeroing and `Reset` keeping every
  counter, backward integration, root locations and directions, every recoverable-failure path and
  every option key pass their tests;
- a level-4 log of a run parses with `suntools.logs` and every `begin-*` has its `end-*`.

### 10.4 M3: hardening and the validation matrix

#### Goal

the suite that later milestones run unchanged (P1), green in the complete CI matrix.

#### Deliverables

`examples/firkode/C_parallel` with `NVECTOR_PARALLEL` subvectors and
`test/unit_tests/firkode/C_parallel` (H9, H14, H15 on 1 and 4 ranks); a gated CUDA or HIP
smoke test; a gated batched-kinetics test with many independent Robertson systems in one `N_Vector`
and a batched block solver (`SUNLinSol_KokkosDense` or `SUNLinSol_MagmaDense`)
`[fact: include/sunlinsol]`; `test/unit_tests/firkode/CXX_serial` developer tests with answer files;
answer files for
every diffed test in the `sundials-codes/answers` repository for single, double and extended
precision `[fact: test/env/docker.sh sets SUNDIALS_TEST_ANSWER_DIR]`; `sundials_install_examples`;
Jenkins and GitLab registration; the user guide complete for the v0.1.0 API with the mode of every
function; a developer page `doc/superbuild/source/developers/packages/firkode` with the step-attempt
flowchart.

#### Done when

`ctest -R firk` is green in every CI configuration including the answers path; in fixed-step mode
with a direct block solver the 1-rank and 4-rank MPI runs agree to $10^{-12}$ relative, the only
order-dependent operations being the allreduce sums; in adaptive mode they are compared by step count
(within 2%) and by error against the reference solution, not bitwise, because the reduction order
perturbs `dsm` at roundoff and can flip an accept or reject decision; neither run performs a clone
per step; every test header states its pass criterion; every hazard row of Table 11.3 tagged M0–M3
has a test.

### 10.5 Release criteria for v0.1.0

- M0–M3 closed.
- `FIRKodeSetNumStages` and `SetOrder` accept $s = 3$ only and reject the rest with
  `FIRK_ILL_INPUT`; no EXPERT hook of v0.4 is declared.
- The guide documents every function with its mode; `CHANGELOG.md` and `RecentChanges.rst`
  announce the package as experimental with independent versioning; `CITATIONS.md` has an entry.
- The measurement table of Section 11.4 for $s = 3$ (FGMRES iterations per Newton iteration on the
  problem set) is recorded in the guide: the mathematics chapter, from
  `benchmarks/firkode/firk_klin` (2.5–4.2 with a direct block solver, 12–14 with a Jacobi-
  preconditioned matrix-free solver on the heat equation) `[fact]`.

### 10.6 M4 (v0.2): the Radau IIA family

#### Goal

every $1 \le s \le 9$ on the same code path, in EXPERT mode (P4: the Krylov counts per $s$
are measured before the family is exposed).

#### Deliverables

the paths for general $s$ exercised ($s = 1$ with no stack and the block solver on the
Newton system; even $s$ with the heuristic $\gamma_0$), `SetNumStages`, `SetOrder`, `GetNumStages`,
`GetCurrentTable` unlocked, every M1–M2 test parametrized over $s$, the measurement script of Section
11.4, extended-precision table checks.

#### Done when

all tests pass for $s = 1, \dots, 9$ in double precision and $s \le 4$ in single;
fixed-step order $2s - 1$ on Kværnø–Prothero–Robinson and non-stiff Prothero–Robinson and at least $s$
on stiff Prothero–Robinson for every $s$; estimate order $p = s$ measured for every $s$;
`firk_analytic` runs for $s \in \{1, 2, 3, 5, 7, 9\}$ with answer files; the $\gamma_0$ literals
for even $s$ agree with an independent `mpmath` root of the Padé denominator; the $k_{lin}$ table
is recorded in this document and in the guide, and the Krylov defaults are revised if it says so.

### 10.7 M5 (v0.3): in situ order selection

#### Deliverables

`firkode_order.c` (Section 4.11, Listing 8.7), `SetOrderRange`, `GetNumOrderChanges`,
`GetNumStepsPerOrder`, the BEGINNER default $3 \le s \le 7$, the current $s$ in every log record.

#### Done when

on Robertson, van der Pol and the Oregonator the variable-order run needs at most 10% more
RHS evaluations than the best fixed odd $s$ and never fails where fixed $s = 3$ succeeds; on a linear
problem with a prescribed contractivity profile the order sequence equals RADAU's with `NSMIN = 3`,
`NSMAX = 7`; a stress test that forces a change every 10 steps is ASan-clean and keeps `GetDky` and
the root locations consistent across changes (H23, H33).

### 10.8 M6 (v0.4): EXPERT stage-solver hooks

#### Deliverables

`SetTable`; `SetStagePreconditioner` (replaces the single-shift block preconditioner
on the stacked vector) and `SetStageLinearSolver` (replaces the internal SPFGMR by a user solver for
the whole $sN$ system with the stacked operator as its `ATimes`), which admit monolithic multigrid
[2024K] and stage-parallel solvers [2024M+]; `SetReusePolicyFn`; `SetNlsConvTestFn`, `SetNlsNormFn`,
`SetPredictorFn`; `SetDomEigEstimator`, `SetDomEigFrequency` (as LSRKStep
`[fact: include/arkode/arkode_lsrkstep.h]`), `SetSpectralRadiusFn`, used for a bound on the initial
step, the Krylov dimension, a cap on the order range, and an optional separate preconditioner shift
`SetPrecShift`, which costs a second setup.

#### Done when

each hook has a unit test showing it is called with the documented arguments and the
default path is bit-identical when the hook is unset; BEGINNER mode rejects every hook with
`FIRK_ILL_INPUT`; the guide lists the mode of every setter.

### 10.9 Optional milestones

**M7: formal certification (optional).** `verification/lean/gen_firkode_tables.py` reads the
generated `firkode_tables.c`; a core-Lean package certifies each literal table against $C(s)$,
$B(2s-1)$, $e_s = (-1)^s \frac{\gamma_0}{s}$, the $\gamma_0$ definition and the dense-output coefficients
in exact rational arithmetic to half a unit in the last printed digit; a Mathlib package holds the
general theorems (Section 11.5); the axiom audit admits only `propext`, `Classical.choice` and
`Quot.sound`; a CI job runs the core package and `--check`. The pilot in `verification/lean`
demonstrates feasibility: it certifies the same claims for tables transcribed by hand. Done when both
packages build on the generated sources and the claim table of Section 11.5 is current.

**M8: `SUNStepper` glue (optional, last; P3).** `firkode_sunstepper.c` sets every `SUNStepper_Set*Fn`
operation (`Evolve`, `OneStep`, `FullRhs`, `ReInit`, `Reset`, `ResetCheckpointIndex`,
`SetStopTime`, `SetStepDirection`, `SetForcing`, `GetNumSteps`, `Destroy`), because
`SUNStepper_Create` leaves some operation pointers uninitialized `[fact:
src/sundials/sundials_stepper.c]` (H12); forcing terms are added inside the RHS wrapper; only
negative codes map to errors. Done when every operation responds in a unit test.

**M9: time-dependent** $M$ **(optional).** Unknowns become the stage derivatives $K_j$ with
$M(t_n + c_j h) K_j = f(t_n + c_j h, y_n + h \sum_l a_{jl} K_l)$ and $Z = hAK$; the Newton matrix
$I_s \otimes M(t_n) - hA \otimes J$ keeps the single-shift preconditioner; the estimate filter uses
$M(t_n)$ and its order must be re-derived; `time_dep = SUNTRUE` as in ARKODE, one mass setup per step
and $s$ mass evaluations per Newton iteration. A singular $M$ stays out: it needs consistent
initialization and index-aware error scaling [1989HLR]. Test: a time-dependent change of variables of
an explicit ODE against the explicit form.

**M10: singly implicit tables (optional).** A SIRK method of stage order $s$ is the collocation
method on $c_i = \lambda \xi_i$ with $\xi_i$ the zeros of the Laguerre polynomial $L_s$ [1976N,
1978B]; $A$ has the single eigenvalue $\lambda$ [1979B] and the order is at most $s + 1$ [1977NW];
L-stable choices with $\frac{1}{\lambda}$ a zero of $L_s$ exist for $s \le 8$ except $s = 7$ [2009B], with
abscissae beyond 1 for $s \ge 3$. With $\gamma_0 = \lambda$ the preconditioned operator is identity
plus nilpotent, so FGMRES with exact block solves terminates in at most $s$ iterations: tables
through `SetTable` and an embedded estimator suffice, with a `tstop` guard for abscissae beyond the
step. The same structure gives a preconditioner $I_s \otimes M - hB \otimes J$ with a singly implicit
$B \approx A$ for Radau IIA [1983CB], to be measured against Table 4.1; methods with distinct real
eigenvalues [1989K, 1997B] would decouple the stages exactly at a further loss of order.

**M11: resize (optional).** `FIRKodeResize(mem, ynew, hscale, t0, resize_fn, resize_data)` with the
semantics of `ARKodeResize` `[fact: include/arkode/arkode.h]`: the stage stack, the weight view, the
Newton and FGMRES work vectors and the history are re-cloned from `ynew`, or transformed in place by
the user's `resize_fn` when one is given; the history is invalidated, so the next step uses the
trivial predictor; the user reattaches the linear and mass solvers and the vector tolerances, as
ARKODE requires. A one-step method keeps no multistep history, so nothing else changes. Test: a 1-D
heat equation on a mesh refined every 10 steps against the run on the fine mesh throughout.

**M12: relaxation (optional).** `SetRelaxFn(mem, rfn, rjac)` with ARKODE's semantics and setters
`[fact: include/arkode/arkode.h, ARKodeSetRelax*]`. The relaxed update
$y_{n+1} = y_n + r \sum_j d_j Z_j$ needs the direction $\sum_j d_j Z_j$, which the step forms anyway
(Section 4.3), and a scalar root solve for the relaxation factor $r$ [2019K, 2020R+]; the accepted
step then has length $r h$, with $r$ confined to a band around 1 as in ARKODE
`[fact: doc/arkode/guide/source/Mathematics.rst]`. Algebraic stability already makes the unrelaxed
method contractive (Section 4.2); relaxation adds exact conservation of one convex functional. Test: a
quadratic invariant of a stiff nonlinear oscillator held to roundoff.

**M13: discrete adjoint (optional).** `SUNAdjointStepper` and `SUNAdjointCheckpointScheme` exist for
explicit Runge–Kutta methods `[fact: include/sundials/sundials_adjointstepper.h]`. The discrete
adjoint of a Runge–Kutta method is a Runge–Kutta-type method with the transposed coupling [2006S]; for
Radau IIA it needs the stage increments of every forward step from the checkpoint scheme and a stage
solver for $I_s \otimes M^T - h A^T \otimes J^T$, that is a transposed block solve from the user. Not
before the API freeze (Chapter 12).

## 11. Verification plan

### 11.1 Test taxonomy

| Kind                           | Where                                              | Registration                                                                        | Passes on                              | Precision                                      |
| ------------------------------ | -------------------------------------------------- | ----------------------------------------------------------------------------------- | -------------------------------------- | ---------------------------------------------- |
| C unit tests                   | `test/unit_tests/firkode/C_serial`                 | `add_test`, linking `sundials_firkode_obj` for internal symbols                     | exit code                              | every precision; $s \le 4$ in single           |
| C++ developer tests            | `test/unit_tests/firkode/CXX_serial`               | `sundials_add_test(... EXAMPLE_TYPE develop ...)` with answer files                 | exit code and answer-file diff         | `exclude-single` where the output cannot match |
| examples with self-checks      | `examples/firkode/C_serial`, `C_parallel`          | `sundials_add_test` with answer files; `exclude-single` where needed                | exit code and diff                     | double and extended                            |
| Krylov-path tests              | both trees                                         | `NODIFF` or exit-code only                                                          | exit code                              | all                                            |
| MPI and GPU tests              | `examples/firkode/C_parallel`, gated GPU directory | `MPI_NPROCS 4`; GPU behind the backend option                                       | exit code                              | double                                         |
| measurement scripts            | `test/` or `scripts/`                              | not run by CTest; produce the tables of Section 11.4                                | –                                      | double                                         |

*Table 11.1. Test taxonomy. Answer files for diffed tests must exist in the `sundials-codes/answers`
repository for every precision CI runs, because CI overrides the per-test answer directory
`[fact: cmake/macros/SundialsAddTest.cmake, test/env/docker.sh]`. Krylov paths are not diffed because
iteration counts differ across platforms.*

### 11.2 Problem set

| Problem                                                                                               | Purpose                                                                       | Used in |
| ----------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------- | ------- |
| Dahlquist $\dot y = \lambda y$, real and $2 \times 2$ rotation blocks                                 | one-step reproduction of $R(z)$, L-stability, Newton count on linear problems | M1      |
| Prothero–Robinson, $\lambda = -1$ and $-10^4$                                                         | fixed-step order; stage-order reduction                                       | M1, M4  |
| Kværnø–Prothero–Robinson                                                                              | fixed-step order on a coupled stiff and non-stiff pair                        | M1, M4  |
| polynomial solutions of degree $\le s + 1$                                                            | dense output and its derivatives, negative control                            | M1      |
| analytic problem $\dot y = \lambda y + \frac{1}{1+t^2} - \lambda \arctan t$                           | example, answer files per $s$                                                 | M1, M4  |
| Robertson kinetics, to $t = 10^{11}$                                                                  | stiff adaptive accuracy and cost; rootfinding; RADAU5 step counts             | M2, M5  |
| van der Pol, $\mu = 1000$                                                                             | stiff adaptive accuracy; RADAU5 step counts                                   | M2, M5  |
| Oregonator                                                                                            | RADAU5 step counts; order selection                                           | M2, M5  |
| 2-species diurnal kinetics PDE, $10 \times 10$                                                        | matrix-free with user $Jv$ and block-diagonal preconditioner                  | M2      |
| 1-D heat equation, P1 finite elements with mass matrix                                                | banded $M - \gamma J$ and $M$                                                 | M2      |
| constant $M = [[2, 1], [1, 3]]$ against $\dot y = M^{-1} f$                                           | mass interface, dense and PCG                                                 | M1      |
| 2-D heat and advection–diffusion, $N = 10^4$ to $10^6$, band, PCG or AMG blocks                       | $k_{lin}$ measurements; MPI weak scaling                                      | M3, M4  |
| linear problem with a prescribed contractivity profile                                                | the order sequence of RADAU                                                   | M5      |
| 1-D heat equation with time-dependent Dirichlet data                                                  | temporal order on a parabolic PDE against Table 4.3                           | M2, M4  |
| batched Robertson, $10^3$ to $10^5$ independent copies in one vector                                  | batched block solver on a GPU; kernel-launch cost of the stack (D4)           | M3      |
| 2-D Brusselator reaction–diffusion (`ark_brusselator1D` precedent `[fact: examples/arkode/C_serial]`) | $k_{lin}$ on a nonlinear reaction–diffusion PDE                               | M4      |
| stiff problems of the Bari IVP test set [2012MCS]                                                     | cross-check of the stiff validation beyond the three RADAU5 problems          | M3, M5  |

*Table 11.2. Problem set.*

### 11.3 Hazard checklist

Each row is a requirement and the test that enforces it; the tag names the earliest milestone that
must close it. `test/unit_tests/firkode/README.md` maps each row to the test files.

| ID       | Requirement                                                                                                                                                                                                                                                                        | Test                                                                                                                                                                    |
| -------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| H1 (M2)  | When the Newton rate of an accepted step exceeds the bad-convergence rate, the next setup re-evaluates $J$ and refactors, under both policies                                                                                                                                      | a counting Jacobian mock: `nje` increments after a slow step                                                                                                            |
| H2 (M1)  | The dense and band difference-quotient Jacobians restore the perturbed component after a recoverable RHS failure and return `FIRKLS_JACFUNC_RECVR`                                                                                                                                 | RHS returns 1 when component $j$ is perturbed; `yn` unchanged bitwise; `Evolve` continues                                                                               |
| H3 (M2)  | `Reset` keeps every counter (integrator, FIRKLS, stage, mass, root); `ReInit` zeroes them all                                                                                                                                                                                      | read every getter before and after                                                                                                                                      |
| H4 (M1)  | Stage-solver setters called after the first step never free a live object; they take effect at the next reallocation or return `FIRK_ILL_INPUT`                                                                                                                                    | `SetStageSolverMaxl` between two `Evolve` calls, then `Evolve`; ASan                                                                                                    |
| H5 (M2)  | `SetUserData` propagates to the root, preconditioner, `jtimes`, mass and weight callbacks whenever it is called                                                                                                                                                                    | set after `RootInit` and `SetPreconditioner`; callbacks assert pointer identity                                                                                         |
| H6 (M1)  | `jtsetup` runs exactly once per Newton iteration that reaches a linear solver (for $s = 1$ the $k \ge 2$ shortcut of Listing 8.4 skips both, as in CVLS), and once per estimate solve only with a matrix-free block solver                                                         | `njtsetup` equals `nni` plus the estimate solves for $s = 3$, and the block solver calls of a recording mock for $s = 1$                                                |
| H7 (M1)  | The difference-quotient $Jv$ with $\lVert v \rVert = 0$ returns $Jv = 0$ without a division                                                                                                                                                                                        | a unit test of the stacked operator with a zero stage block                                                                                                             |
| H8 (M2)  | $J$ is never re-evaluated at the same $(t_n, y_n)$: after a convergence failure with a current $J$ only the factorization is redone                                                                                                                                                | `nje == 1` across a step with two forced convergence failures                                                                                                           |
| H9 (M1)  | No forbidden operation of Table 8.8 is ever applied to a stacked vector; `Init` rejects a `y0` with a communicator but without `nvdotprodlocal` and `nvwsqrsumlocal`                                                                                                               | a mock subvector whose forbidden operations trap; a mock with a communicator and no local reductions; 1-rank and 4-rank fixed-step MPI agreement                        |
| H10 (M1) | Lengths come from `N_VGetLength` only                                                                                                                                                                                                                                              | covered by the H9 mock                                                                                                                                                  |
| H11 (M1) | `SUNLinSolSetZeroGuess(.., SUNTRUE)` precedes every `SUNLinSolSolve` (stacked, block, mass)                                                                                                                                                                                        | a mock linear solver asserting the flag on each call                                                                                                                    |
| H12 (M8) | Every `SUNStepper` operation is set; positive FIRKODE codes are not errors in the glue                                                                                                                                                                                             | call every `SUNStepper_*` operation; `FIRK_TSTOP_RETURN` and `FIRK_ROOT_RETURN` map to success                                                                          |
| H13 (M2) | `dsm` is clamped to at least $10^{-10}$ before `EstimateStep` (the estimate also floors it, as RADAU5's `ESTRAD`)                                                                                                                                                                  | a polynomial exact solution: the step after the first grows by $\kappa (10^{-10})^{-1/4}$, the I-controller proposal at the floor, and the later ones by the growth cap |
| H14 (M3) | A NULL return of the MPIManyVector constructor is handled; no clone during the steps (the stacks are cloned in `Init`, the stacked SPFGMR when the linear solver is initialized at the first `Evolve` or after a stage-solver setter, H4); identical allocation order on all ranks | a clone counter equal to 0 per step; a 4-rank run completes                                                                                                             |
| H15 (M1) | No `N_VLinearCombination` or `N_VScaleAddMulti` on a stacked vector outside SPFGMR; modified Gram–Schmidt by default                                                                                                                                                               | an operation-counting mock on the stack                                                                                                                                 |
| H16 (M1) | The weight view is never cloned, is destroyed without ownership, and is rebuilt at an order change                                                                                                                                                                                 | ASan in M1 and in the M5 stress test                                                                                                                                    |
| H17 (M1) | The estimate solve never takes the small-residual shortcut; only the `dsm` floor applies                                                                                                                                                                                           | dense against iterative estimate within 50%                                                                                                                             |
| H18 (M1) | `SUNLS_RES_REDUCED` is accepted only at Newton iteration 1                                                                                                                                                                                                                         | a mock solver returning it at $k = 2$: recoverable, step reduced                                                                                                        |
| H19 (M2) | A recoverable RHS failure in a stage reduces $h$ by 0.25 and retries; `FIRK_REPTD_RHSFUNC_ERR` after `maxncf`; an unrecoverable failure returns `FIRK_RHSFUNC_FAIL`                                                                                                                | a domain-restricted RHS                                                                                                                                                 |
| H20 (M2) | A RHS failure in the refilter keeps the first estimate                                                                                                                                                                                                                             | a RHS that fails only at $y_n + \mathrm{err}$                                                                                                                           |
| H21 (M2) | Stop time: clamp $(t_{stop} - t_n)(1 - 4u)$; `SetStopTime` after passing it is `FIRK_ILL_INPUT`; auto-clear on `FIRK_TSTOP_RETURN`; the fixed step resumes after a shortened step                                                                                                  | `firk_test_tstop`, adaptive and fixed                                                                                                                                   |
| H22 (M2) | Rootfinding on the collocation polynomial with `ttol = 100u(\lvert t \rvert + \lvert h \rvert)`; directions; inactive-root warning; a root at $t_0$                                                                                                                                | `firk_test_root`, `firk_test_rootfind`                                                                                                                                  |
| H23 (M5) | After an order change `GetDky` uses the history's own table; roots are located across the change                                                                                                                                                                                   | the forced-change stress test                                                                                                                                           |
| H24 (M0) | Literals compile in single, double and extended precision under `-Wconversion` and `-Wdouble-promotion`; every constant through `SUN_RCONST`                                                                                                                                       | the CI matrix                                                                                                                                                           |
| H25 (M0) | C99, no variable-length arrays, fixed `FIRK_MAX_STAGES` arrays; `-Wvla -Walloca` clean                                                                                                                                                                                             | the all-warnings CI job                                                                                                                                                 |
| H26 (M3) | Every diffed test has answer files for three precisions in the answers repository                                                                                                                                                                                                  | the answers-diff CI path                                                                                                                                                |
| H27 (M1) | The vector budget equals Table 8.7                                                                                                                                                                                                                                                 | `firk_test_workspace` (live clones of `y0`)                                                                                                                             |
| H28 (M1) | Mass: $J$ and $M$ both matrix-based or both matrix-free; `time_dep` rejected; the constant $M$ set up once; product and solve counters as specified                                                                                                                                | `firk_test_mass`, `firk_test_reset`                                                                                                                                     |
| H29 (M1) | `jok = !jbad` and `jcurPtr` honored for the user's `psetup` and `linsys` under both policies                                                                                                                                                                                       | a recording preconditioner mock                                                                                                                                         |
| H30 (M2) | The dead band keeps $\gamma$ and the factorization (`nsetups` unchanged)                                                                                                                                                                                                           | a smooth problem with setups counted                                                                                                                                    |
| H31 (M1) | Log records are balanced on every exit path                                                                                                                                                                                                                                        | a level-4 run parsed by `suntools.logs`                                                                                                                                 |
| H32 (M2) | Warnings for $t + h = t$ are capped at `mxhnil` and go through the logger                                                                                                                                                                                                          | a tiny-step test                                                                                                                                                        |
| H33 (M5) | An allocation failure during an order change leaves the state intact and returns `FIRK_ORDER_CHANGE_FAIL`                                                                                                                                                                          | a fault-injecting `N_VClone` wrapper                                                                                                                                    |
| H34 (M2) | Every example returns a nonzero exit code when the solver fails                                                                                                                                                                                                                    | review plus CI exit codes                                                                                                                                               |
| H35 (M1) | A zero error weight is `FIRK_ILL_INPUT`; weights are evaluated at $y_n$ once per step                                                                                                                                                                                              | `SStolerances(.., 0)` with $y = 0$                                                                                                                                      |
| H36 (M2) | `SUNAdaptController_EstimateStep` is called at most once per step attempt (not for a rejected attempt in BEGINNER mode) and before `UpdateH`, so that the BEGINNER default is Gustafsson's predictive controller and not the I-controller                                          | a recording mock controller checks the call order; on a problem with varying error the accepted-step proposal differs from $\kappa h \, \mathrm{dsm}^{-1/(s+1)}$        |
| H37 (M1) | With an iterative block solver the stacked solve applies only the user's `psolve`; the user's Krylov iteration runs only for the $s = 1$ Newton solve and the estimate                                                                                                             | the user solver's `nli` is unchanged by stacked solves; `nps` grows by $s$ per FGMRES iteration                                                                         |

*Table 11.3. Hazard checklist.*

### 11.4 Measurement protocols

Sentences that assert a number this document does not have may be written only after the
corresponding measurement exists.

| Blocked claim                                                                         | Protocol                                                                                                                                                                                                                                                                                                                        | Gate                                          |
| ------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------- |
| "FGMRES needs $k_{lin}$ iterations per Newton iteration at stage count $s$"           | 2-D heat and advection–diffusion with $N = 10^4$ to $10^6$ and band, PCG and AMG block solvers; Robertson, van der Pol, Oregonator; one advection-dominated case; record mean and maximum $k_{lin}$, restarts and failures for $s = 1, \dots, 9$                                                                                | M4; $s = 3$ recorded at v0.1.0                |
| "step counts match RADAU5 within 20%"                                                 | equal effective tolerances after $0.1 \mathrm{rtol}^{\frac{2}{3}}$, $k_1 = k_2 = 1$, the RADAU5 policy, RADAU5's initial step $10^{-6}$ and the final time as a stop time; reference sequences and solutions from `radau5.f` with the script committed; errors within twice RADAU5's or below half the tolerance (Section 10.3) | M2                                            |
| "estimate order $p = s$" and "order at least $s$ on stiff Prothero–Robinson"          | slopes on smooth and stiff Prothero–Robinson for every $s$                                                                                                                                                                                                                                                                      | M2, M4                                        |
| the memory formula                                                                    | count `N_VClone` calls against $11 + s(2\mathrm{maxl} + 10)$                                                                                                                                                                                                                                                                    | M1                                            |
| "one allreduce per reduction with MPIManyVector"                                      | count collectives per step with MPIManyVector over `NVECTOR_PARALLEL` (a plain ManyVector over `NVECTOR_PARALLEL` would reduce once per stage, and its max norm rank-locally); weak scaling on 2-D heat with hypre                                                                                                              | M3                                            |
| "order about $s + 1$ on parabolic PDEs with time-dependent boundary data" (Table 4.3) | slopes on the 1-D heat equation with time-dependent and with homogeneous Dirichlet data, discrete $L^2$ norm in space, $s = 3$ first and then every $s$                                                                                                                                                                         | M2, M4                                        |
| "an adaptive forcing term reduces the FGMRES work" (Section 5.10)                     | FGMRES iterations per step and RHS evaluations with fixed and with adaptive $\epsilon_{stk}$ on the problem set of this table                                                                                                                                                                                                   | M6                                            |
| work–precision against CVODE BDF and ARKStep ESDIRK                                   | the stiff set of [1996HW, §IV.10] plus one PDE at scale                                                                                                                                                                                                                                                                         | a later performance report, not this document |

*Table 11.4. Measurement protocols.*

### 11.5 Formal verification (optional)

The mathematics behind the tables is small enough to be machine-checked. The pilot in
`verification/lean` has two packages: a core-Lean package of 39 kernel-checked certificates
(`decide +kernel`, no Mathlib, builds in seconds) and a Mathlib package of 337 theorems; its axiom
audit fails the build on any axiom beyond `propext`, `Classical.choice` and `Quot.sound` [2021MU,
2020M]. Table 11.5 lists what it establishes for the Radau IIA tables with $s \le 9$ and what remains
open. Milestone M7 re-targets the certificates at the generated literals of M0.

| Claim                                                                                                                   | Status for $s \le 9$                                                  |
| ----------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------- |
| the nodes are the roots of $P_s(2c-1) - P_{s-1}(2c-1)$, real, simple, increasing, $c_s = 1$                             | proved                                                                |
| $C(s)$, $B(2s-1)$, $D(s-1)$; `CheckOrder` returns $q = 2s - 1$ in exact arithmetic                                      | proved; the order itself cites Butcher's theorem                      |
| the smallest violated condition for $s = 9$ is $B(18)$ with residual $9.4 \cdot 10^{-11}$                               | proved                                                                |
| $d = \mathbf e_s$ and $R(\infty) = 0$                                                                                   | proved                                                                |
| the stability function is the $(s-1, s)$ Padé approximant; A- and L-stability                                           | proved                                                                |
| $\gamma_0 = \frac{1}{U_1}$ with $U_1$ the unique real eigenvalue of $A^{-1}$ (odd $s$); no real eigenvalue for even $s$ | proved; the numeric value for even $s$ is not certified               |
| $e_j = -\frac{\gamma_0 l_j(0)}{c_j}$, $e_s = (-1)^s \frac{\gamma_0}{s}$, RADAU5's weights for $s = 3$                   | proved                                                                |
| the estimate is exact on polynomials of degree $\le s$ (order $p = s$)                                                  | proved; the $O(h^{s+1})$ bound for smooth $f$ is cited                |
| dense output interpolates the stages and ends at $y_{n+1}$; the derivative identities behind `GetDky`                   | proved                                                                |
| the statement that $\gamma_0$ minimizes the stiff-limit deviation                                                       | false in the worst case for $s = 3, 4$; not claimed                   |
| algebraic stability $BA + A^T B - b b^T \succeq 0$ (Section 4.2)                                                        | open; a finite check on each literal table, proposed as Lean item R24 |
| the Newton test, the step and order rules, rootfinding, the initial step, floating-point rounding, the C code           | not covered                                                           |

*Table 11.5. Formal verification: the claims of the pilot, to be re-established on the generated
literals in M7. Trust model: kernel-checked theorems; generated data; Butcher's theorem cited; C
semantics and rounding outside the scope.*

## 12. Open questions and decision log

1. **Choice of** $\gamma_0$. The rule of Section 4.5 is exact for the one real mode of odd $s$ but not optimal
   for the bound of Table 4.1. Minimizing the bound over the shift gives a larger value, about 0.37
   instead of 0.275 for $s = 3$ (bound 0.751 instead of 0.766) and about 0.20 instead of 0.086 for
   $s = 9$ (0.941 instead of 0.969) `[derived]`. The minimum is flat, any other shift loses the exact
   mode, and $\gamma_0$ also defines the estimate and its filter, so a separate preconditioner shift
   costs a second setup (`SetPrecShift`, v0.4). Decide from the measured iteration counts of M4.
2. **Krylov defaults.** $\mathrm{maxl} = \max(5, 3s)$ with one restart is below the worst-case bound
   for $s \ge 6$ (Section 4.5). Revise after M4.
3. **Alternative stage solvers** (D2) for the v0.4 hooks: block-triangular preconditioners
   [2021R+, 2024ADN]; one real preconditioner per conjugate pair [2022SKPD]; stage-parallel solves
   [2017PP, 2024M+], which need an MPIManyVector over sub-communicators, one block solver per group
   and an exchange for the coupling through $A$; monolithic multigrid [2005VV, 2024K]; low-
   synchronization Gram–Schmidt and pipelined GMRES for the stacked solve [2021S+, 2013G+];
   parallelism across the method in general [1988G, 1995JN].
4. **Even** $s$ (D7). A heuristic shift and no exactly preconditioned mode; EXPERT only;
   `[to measure]`.
5. **Order rule.** RADAU's thresholds plus an uncalibrated $k_{lin}$ veto; calibrate on the M4 data.
6. **Rootfinding** (D13). A shared SUNDIALS rootfinding module would remove the duplication of
   Table 5.2; proposed for after 1.0.
7. **Package or stepper** (D1). Re-evaluate before 1.0 against its revisit condition.
8. **GPU layout** (D4). A ManyVector launches $s$ kernels per operation; a contiguous multivector
   with a block product $JV$ would raise arithmetic intensity. Decide from measurements on device
   vectors [2021B+].
9. **API freeze.** Settle names and semantics before 1.0, in particular whether the stage-solver
   options become a stage-solver object.
10. **Adaptive forcing terms** (Section 5.10). An Eisenstat–Walker forcing term [1996EW] for
    $\epsilon_{stk}$ is an EXPERT option of M6; decide the default from the M6 measurement.
11. **Krylov recycling.** Consecutive stage systems differ only through $Z$ and $h$, so recycling
    part of the FGMRES subspace across Newton iterations and steps [2006P+] could lower $k_{lin}$ for
    $s \ge 5$; it needs a recycling variant of SPFGMR or a user stage solver (M6) `[to measure]`.
12. **Mixed precision.** FGMRES tolerates an inexact preconditioner, so the block factorization may
    be computed and applied in lower precision: [2026LBMC] does this for Radau IIA with the block
    preconditioner of [2021R+], and [2025DKGK] measures mixed-precision implicit solves on GPUs.
    SUNDIALS has one `sunrealtype`, so the path is a user `SUNLinearSolver` `[to measure]`.
13. **Stiffness-aware control.** `SUNDomEigEstimator` `[fact: include/sundomeigest]` estimates
    $\rho(J)$; with the stiffness characterizations of [2015SJC], $h \rho(J)$ could cap the order
    range of Section 4.11 and bound the initial step, as the v0.4 hooks foresee. Calibrate on the M4
    data.
14. **Positivity.** Constraints stay out (Section 3.5); if users ask, step rejection on constraint
    violation [2021NRK] is the route, since no Runge–Kutta method of order above one preserves
    positivity unconditionally [2021NRK].
15. **The M2 error gate.** With the rule of Section 4.4 the step sequences of FIRKODE and `radau5.f`
    agree almost step for step on the M2 problems (van der Pol at $10^{-6}$: 486 against 489 steps,
    211 against 212 through the jump at $t = 1.614$), but the ratios of the final-time errors range
    from 0.75 to 23 and swing a hundredfold under small changes of the settings, since both errors
    lie well below the tolerance and are set by the path through the stiff transitions; hence the
    floor of half the tolerance in Section 10.3
    `[fact: test/unit_tests/firkode/CXX_serial/firk_test_radau5_steps.cpp]`. A sharper comparison
    would average over perturbed tolerances or compare the errors at the end of the transition.
    The rule also tightens the stacked FGMRES tolerance $\epsilon_{stk} \epsilon_{\text{nls}}$ at
    tight tolerances (the finite-element heat problem of `firk_test_mass_fem_krylov` at
    $\mathrm{rtol} = 10^{-6}$: about 3.5 FGMRES iterations per Newton iteration with banded
    blocks, against 2.5 with $\epsilon_{\text{nls}} = 0.1$) `[fact]`, which strengthens the case
    of question 10.
16. **Known issues of the general-$s$ path** (reached through the internal stage setter only in
    v0.1.0), from the benchmark report of branch `work/firkode-benchmark`, Section 8: at
    $\mathrm{rtol} = 10^{-2}$ the $s = 5$ and $s = 7$ methods end in `FIRK_CONV_FAILURE` at a tiny
    step on HIRES, Robertson and the Oregonator (the Newton iteration of the general-$s$ path at
    loose tolerances, v0.2); at $s = 9$ the stacked FGMRES reaches its cap of $2 \cdot 3s$
    iterations in a few solves (question 2); Robertson at $\mathrm{atol} = \mathrm{rtol} = 10^{-2}$
    ends in `FIRK_ERR_FAILURE` once $v$ turns negative (an absolute tolerance far above the
    component, not a defect). The benchmark's step-size collapse at loose tolerances for $s = 3$,
    the error-test cascade on van der Pol and the step collapse on fine finite-element meshes were
    fixed by the minimum step of Section 4.7, the refilter and the units conversion of Section
    5.10, and `firk_test_vdp` and `firk_test_mass_fem` guard them `[fact]`.

## References

- [1963D] G. G. Dahlquist. A special stability problem for linear multistep methods. *BIT*
  3:27–43, 1963. [doi:10.1007/BF01963532](https://doi.org/10.1007/BF01963532)
- [1964Ba] J. C. Butcher. Implicit Runge–Kutta processes. *Math. Comp.* 18:50–64, 1964.
  [doi:10.1090/S0025-5718-1964-0159424-9](https://doi.org/10.1090/S0025-5718-1964-0159424-9)
- [1964Bb] J. C. Butcher. Integration processes based on Radau quadrature formulas. *Math. Comp.*
  18:233–244, 1964.
  [doi:10.1090/S0025-5718-1964-0165693-1](https://doi.org/10.1090/S0025-5718-1964-0165693-1)
- [1969A] O. Axelsson. A class of A-stable methods. *BIT* 9:185–199, 1969.
  [doi:10.1007/BF01946812](https://doi.org/10.1007/BF01946812)
- [1974LR] P. Lesaint, P.-A. Raviart. On a finite element method for solving the neutron transport
  equation. In C. de Boor (ed.), *Mathematical Aspects of Finite Elements in Partial Differential
  Equations*, pp. 89–123, Academic Press, 1974.
  [doi:10.1016/B978-0-12-208350-1.50008-X](https://doi.org/10.1016/B978-0-12-208350-1.50008-X)
- [1974PR] A. Prothero, A. Robinson. On the stability and accuracy of one-step methods for solving
  stiff systems of ordinary differential equations. *Math. Comp.* 28:145–162, 1974.
  [doi:10.1090/S0025-5718-1974-0331793-2](https://doi.org/10.1090/S0025-5718-1974-0331793-2)
- [1975B] J. C. Butcher. A stability property of implicit Runge–Kutta methods. *BIT* 15(4):358–361,
  1975. [doi:10.1007/BF01931672](https://doi.org/10.1007/BF01931672)
- [1976B] J. C. Butcher. On the implementation of implicit Runge–Kutta methods. *BIT*
  16:237–240, 1976. [doi:10.1007/BF01932265](https://doi.org/10.1007/BF01932265)
- [1976N] S. P. Nørsett. Runge–Kutta methods with a multiple real eigenvalue only. *BIT*
  16:388–393, 1976. [doi:10.1007/BF01932722](https://doi.org/10.1007/BF01932722)
- [1977B] T. A. Bickart. An efficient solution process for implicit Runge–Kutta methods. *SIAM J.
  Numer. Anal.* 14(6):1022–1027, 1977. [doi:10.1137/0714069](https://doi.org/10.1137/0714069)
- [1977NW] S. P. Nørsett, A. Wolfbrandt. Attainable order of rational approximations to the
  exponential function with only real poles. *BIT* 17:200–208, 1977.
  [doi:10.1007/BF01932291](https://doi.org/10.1007/BF01932291)
- [1978B] K. Burrage. A special family of Runge–Kutta methods for solving stiff differential
  equations. *BIT* 18:22–41, 1978. [doi:10.1007/BF01947741](https://doi.org/10.1007/BF01947741)
- [1978WHN] G. Wanner, E. Hairer, S. P. Nørsett. Order stars and stability theorems. *BIT*
  18:475–489, 1978. [doi:10.1007/BF01932026](https://doi.org/10.1007/BF01932026)
- [1979B] J. C. Butcher. A transformed implicit Runge–Kutta method. *J. ACM* 26(4):731–738, 1979.
  [doi:10.1145/322154.322163](https://doi.org/10.1145/322154.322163)
- [1979BB] K. Burrage, J. C. Butcher. Stability criteria for implicit Runge–Kutta methods. *SIAM J.
  Numer. Anal.* 16(1):46–57, 1979. [doi:10.1137/0716004](https://doi.org/10.1137/0716004)
- [1981FSU] R. Frank, J. Schneid, C. W. Ueberhuber. The concept of B-convergence. *SIAM J. Numer.
  Anal.* 18(5):753–780, 1981. [doi:10.1137/0718051](https://doi.org/10.1137/0718051)
- [1982DES] R. S. Dembo, S. C. Eisenstat, T. Steihaug. Inexact Newton methods. *SIAM J. Numer.
  Anal.* 19(2):400–408, 1982. [doi:10.1137/0719025](https://doi.org/10.1137/0719025)
- [1983CB] G. J. Cooper, J. C. Butcher. An iteration scheme for implicit Runge–Kutta methods. *IMA
  J. Numer. Anal.* 3(2):127–140, 1983.
  [doi:10.1093/imanum/3.2.127](https://doi.org/10.1093/imanum/3.2.127)
- [1984DV] K. Dekker, J. G. Verwer. *Stability of Runge–Kutta Methods for Stiff Nonlinear Differential
  Equations.* CWI Monographs 2, North-Holland, 1984.
- [1985FSU] R. Frank, J. Schneid, C. W. Ueberhuber. Order results for implicit Runge–Kutta methods
  applied to stiff systems. *SIAM J. Numer. Anal.* 22(3):515–534, 1985.
  [doi:10.1137/0722031](https://doi.org/10.1137/0722031)
- [1986BHV] K. Burrage, W. H. Hundsdorfer, J. G. Verwer. A study of B-convergence of Runge–Kutta
  methods. *Computing* 36(1–2):17–34, 1986.
  [doi:10.1007/BF02238189](https://doi.org/10.1007/BF02238189)
- [1986SS] Y. Saad, M. H. Schultz. GMRES: A generalized minimal residual algorithm for solving
  nonsymmetric linear systems. *SIAM J. Sci. Stat. Comput.* 7(3):856–869, 1986.
  [doi:10.1137/0907058](https://doi.org/10.1137/0907058)
- [1986SVH] J. M. Sanz-Serna, J. G. Verwer, W. H. Hundsdorfer. Convergence and order reduction of
  Runge–Kutta schemes applied to evolutionary problems in partial differential equations. *Numer.
  Math.* 50(4):405–418, 1986. [doi:10.1007/BF01396661](https://doi.org/10.1007/BF01396661)
- [1987BH] K. Burrage, W. H. Hundsdorfer. The order of B-convergence of algebraically stable
  Runge–Kutta methods. *BIT* 27:62–71, 1987.
  [doi:10.1007/BF01937355](https://doi.org/10.1007/BF01937355)
- [1988G] C. W. Gear. Parallel methods for ordinary differential equations. *Calcolo* 25:1–20,
  1988. [doi:10.1007/BF02575744](https://doi.org/10.1007/BF02575744)
- [1989BH] P. N. Brown, A. C. Hindmarsh. Reduced storage matrix methods in stiff ODE systems.
  *Appl. Math. Comput.* 31:40–91, 1989.
  [doi:10.1016/0096-3003(89)90110-0](https://doi.org/10.1016/0096-3003(89)90110-0)
- [1989HLR] E. Hairer, C. Lubich, M. Roche. *The Numerical Solution of Differential-Algebraic
  Systems by Runge–Kutta Methods.* Lecture Notes in Math. 1409, Springer, 1989.
  [doi:10.1007/BFb0093947](https://doi.org/10.1007/BFb0093947)
- [1989K] S. L. Keeling. On implicit Runge–Kutta methods with a stability function having distinct
  real poles. *BIT* 29:91–109, 1989. [doi:10.1007/BF01932708](https://doi.org/10.1007/BF01932708)
- [1991HS] P. J. van der Houwen, B. P. Sommeijer. Iterated Runge–Kutta methods on parallel
  computers. *SIAM J. Sci. Stat. Comput.* 12(5):1000–1028, 1991.
  [doi:10.1137/0912054](https://doi.org/10.1137/0912054)
- [1992OR] A. Ostermann, M. Roche. Runge–Kutta methods for partial differential equations and
  fractional orders of convergence. *Math. Comp.* 59(200):403–420, 1992.
  [doi:10.1090/S0025-5718-1992-1142285-6](https://doi.org/10.1090/S0025-5718-1992-1142285-6)
- [1993HNW] E. Hairer, S. P. Nørsett, G. Wanner. *Solving Ordinary Differential Equations I:
  Nonstiff Problems,* 2nd ed., Springer, 1993.
  [doi:10.1007/978-3-540-78862-1](https://doi.org/10.1007/978-3-540-78862-1)
- [1993LO] C. Lubich, A. Ostermann. Runge–Kutta methods for parabolic equations and convolution
  quadrature. *Math. Comp.* 60(201):105–131, 1993.
  [doi:10.1090/S0025-5718-1993-1153166-7](https://doi.org/10.1090/S0025-5718-1993-1153166-7)
- [1993S] Y. Saad. A flexible inner-outer preconditioned GMRES algorithm. *SIAM J. Sci. Comput.*
  14(2):461–469, 1993. [doi:10.1137/0914028](https://doi.org/10.1137/0914028)
- [1994G] K. Gustafsson. Control-theoretic techniques for stepsize selection in implicit
  Runge–Kutta methods. *ACM Trans. Math. Softw.* 20(4):496–517, 1994.
  [doi:10.1145/198429.198437](https://doi.org/10.1145/198429.198437)
- [1995JN] K. R. Jackson, S. P. Nørsett. The potential for parallelism in Runge–Kutta methods.
  Part 1: RK formulas in standard form. *SIAM J. Numer. Anal.* 32(1):49–82, 1995.
  [doi:10.1137/0732002](https://doi.org/10.1137/0732002)
- [1995LO] C. Lubich, A. Ostermann. Runge–Kutta approximation of quasi-linear parabolic equations.
  *Math. Comp.* 64(210):601–627, 1995.
  [doi:10.1090/S0025-5718-1995-1284670-0](https://doi.org/10.1090/S0025-5718-1995-1284670-0)
- [1996EW] S. C. Eisenstat, H. F. Walker. Choosing the forcing terms in an inexact Newton method.
  *SIAM J. Sci. Comput.* 17(1):16–32, 1996. [doi:10.1137/0917003](https://doi.org/10.1137/0917003)
- [1996HW] E. Hairer, G. Wanner. *Solving Ordinary Differential Equations II: Stiff and
  Differential-Algebraic Problems,* 2nd ed., Springer, 1996.
  [doi:10.1007/978-3-642-05221-7](https://doi.org/10.1007/978-3-642-05221-7)
- [1997B] C. Bendtsen. A parallel stiff ODE solver based on MIRKs. *Adv. Comput. Math.*
  7(1–2):27–36, 1997. [doi:10.1023/A:1018978331792](https://doi.org/10.1023/A:1018978331792)
- [1997HS] P. J. van der Houwen, J. J. B. de Swart. Parallel linear system solvers for Runge–Kutta
  methods. *Adv. Comput. Math.* 7(1–2):157–181, 1997.
  [doi:10.1023/A:1018990601750](https://doi.org/10.1023/A:1018990601750)
- [1997SS] J. J. B. de Swart, G. Söderlind. On the construction of error estimators for implicit
  Runge–Kutta methods. *J. Comput. Appl. Math.* 86(2):347–358, 1997.
  [doi:10.1016/S0377-0427(97)00166-0](https://doi.org/10.1016/S0377-0427(97)00166-0)
- [1999HW] E. Hairer, G. Wanner. Stiff differential equations solved by Radau methods. *J. Comput.
  Appl. Math.* 111(1–2):93–111, 1999.
  [doi:10.1016/S0377-0427(99)00134-X](https://doi.org/10.1016/S0377-0427(99)00134-X)
- [1999JB] L. O. Jay, T. Braconnier. A parallelizable preconditioner for the iterative solution of
  implicit Runge–Kutta-type methods. *J. Comput. Appl. Math.* 111(1–2):63–76, 1999.
  [doi:10.1016/S0377-0427(99)00132-6](https://doi.org/10.1016/S0377-0427(99)00132-6)
- [2000DGR] A. Dutt, L. Greengard, V. Rokhlin. Spectral deferred correction methods for ordinary
  differential equations. *BIT* 40(2):241–266, 2000.
  [doi:10.1023/A:1022338906936](https://doi.org/10.1023/A:1022338906936)
- [2000J] L. O. Jay. Inexact simplified Newton iterations for implicit Runge–Kutta methods. *SIAM
  J. Numer. Anal.* 38(4):1369–1388, 2000.
  [doi:10.1137/S0036142999360573](https://doi.org/10.1137/S0036142999360573)
- [2000VL] C. F. Van Loan. The ubiquitous Kronecker product. *J. Comput. Appl. Math.*
  123(1–2):85–100, 2000.
  [doi:10.1016/S0377-0427(00)00393-9](https://doi.org/10.1016/S0377-0427(00)00393-9)
- [2003HV] W. Hundsdorfer, J. G. Verwer. *Numerical Solution of Time-Dependent
  Advection-Diffusion-Reaction Equations.* Springer Series in Computational Mathematics 33, Springer,
  2003. [doi:10.1007/978-3-662-09017-6](https://doi.org/10.1007/978-3-662-09017-6)
- [2003S] G. Söderlind. Digital filters in adaptive time-stepping. *ACM Trans. Math. Softw.*
  29(1):1–26, 2003. [doi:10.1145/641876.641877](https://doi.org/10.1145/641876.641877)
- [2004KK] D. A. Knoll, D. E. Keyes. Jacobian-free Newton–Krylov methods: a survey of approaches and
  applications. *J. Comput. Phys.* 193(2):357–397, 2004.
  [doi:10.1016/j.jcp.2003.08.010](https://doi.org/10.1016/j.jcp.2003.08.010)
- [2005H+] A. C. Hindmarsh, P. N. Brown, K. E. Grant, S. L. Lee, R. Serban, D. E. Shumaker,
  C. S. Woodward. SUNDIALS: Suite of nonlinear and differential/algebraic equation solvers. *ACM
  Trans. Math. Softw.* 31(3):363–396, 2005.
  [doi:10.1145/1089014.1089020](https://doi.org/10.1145/1089014.1089020)
- [2005VV] J. Van Lent, S. Vandewalle. Multigrid methods for implicit Runge–Kutta and boundary
  value method discretizations of parabolic PDEs. *SIAM J. Sci. Comput.* 27(1):67–92, 2005.
  [doi:10.1137/030601144](https://doi.org/10.1137/030601144)
- [2006P+] M. L. Parks, E. de Sturler, G. Mackey, D. D. Johnson, S. Maiti. Recycling Krylov subspaces
  for sequences of linear systems. *SIAM J. Sci. Comput.* 28(5):1651–1674, 2006.
  [doi:10.1137/040607277](https://doi.org/10.1137/040607277)
- [2006S] A. Sandu. On the properties of Runge–Kutta discrete adjoints. In *Computational Science –
  ICCS 2006*, LNCS 3994, pp. 550–557, Springer, 2006.
  [doi:10.1007/11758549_76](https://doi.org/10.1007/11758549_76)
- [2006SMN] G. A. Staff, K.-A. Mardal, T. K. Nilssen. Preconditioning of fully implicit Runge–Kutta
  schemes for parabolic PDEs. *Model. Identif. Control* 27(2):109–123, 2006.
  [doi:10.4173/mic.2006.2.3](https://doi.org/10.4173/mic.2006.2.3)
- [2006SW] G. Söderlind, L. Wang. Adaptive time-stepping and computational stability. *J. Comput.
  Appl. Math.* 185(2):225–243, 2006.
  [doi:10.1016/j.cam.2005.03.008](https://doi.org/10.1016/j.cam.2005.03.008)
- [2007MNS] K.-A. Mardal, T. K. Nilssen, G. A. Staff. Order-optimal preconditioners for implicit
  Runge–Kutta schemes applied to parabolic PDEs. *SIAM J. Sci. Comput.* 29(1):361–375, 2007.
  [doi:10.1137/05064093X](https://doi.org/10.1137/05064093X)
- [2009B] J. C. Butcher. Practical Runge–Kutta methods for scientific computation. *ANZIAM J.*
  50(3):333–342, 2009.
  [doi:10.1017/S1446181109000030](https://doi.org/10.1017/S1446181109000030)
- [2010HO] M. Hochbruck, A. Ostermann. Exponential integrators. *Acta Numer.* 19:209–286, 2010.
  [doi:10.1017/S0962492910000048](https://doi.org/10.1017/S0962492910000048)
- [2012MCS] F. Mazzia, J. R. Cash, K. Soetaert. A test set for stiff initial value problem solvers in
  the open source software R: Package deTestSet. *J. Comput. Appl. Math.* 236(16):4119–4131, 2012.
  [doi:10.1016/j.cam.2012.03.014](https://doi.org/10.1016/j.cam.2012.03.014)
- [2013G+] P. Ghysels, T. J. Ashby, K. Meerbergen, W. Vanroose. Hiding global communication latency
  in the GMRES algorithm on massively parallel machines. *SIAM J. Sci. Comput.* 35(1):C48–C71,
  2013. [doi:10.1137/12086563X](https://doi.org/10.1137/12086563X)
- [2015SJC] G. Söderlind, L. Jay, M. Calvo. Stiffness 1952–2012: Sixty years in search of a
  definition. *BIT* 55(2):531–558, 2015.
  [doi:10.1007/s10543-014-0503-3](https://doi.org/10.1007/s10543-014-0503-3)
- [2017CNS] N. J. Curtis, K. E. Niemeyer, C.-J. Sung. An investigation of GPU-based stiff chemical
  kinetics integration methods. *Combust. Flame* 179:312–324, 2017.
  [doi:10.1016/j.combustflame.2017.02.005](https://doi.org/10.1016/j.combustflame.2017.02.005)
- [2017PP] W. Pazner, P.-O. Persson. Stage-parallel fully implicit Runge–Kutta solvers for
  discontinuous Galerkin fluid simulations. *J. Comput. Phys.* 335:700–717, 2017.
  [doi:10.1016/j.jcp.2017.01.050](https://doi.org/10.1016/j.jcp.2017.01.050)
- [2017RN] C. Rackauckas, Q. Nie. DifferentialEquations.jl – A performant and feature-rich
  ecosystem for solving differential equations in Julia. *J. Open Res. Softw.* 5(1):15, 2017.
  [doi:10.5334/jors.151](https://doi.org/10.5334/jors.151)
- [2018S] R. Speck. Parallelizing spectral deferred corrections across the method. *Comput. Vis.
  Sci.* 19(3–4):75–83, 2018.
  [doi:10.1007/s00791-018-0298-x](https://doi.org/10.1007/s00791-018-0298-x)
- [2018SAN] C. P. Stone, A. T. Alferman, K. E. Niemeyer. Accelerating finite-rate chemical kinetics
  with coprocessors: Comparing vectorization methods on GPUs, MICs, and CPUs. *Computer Physics
  Communications* 226:18–29, 2018.
  [doi:10.1016/j.cpc.2018.01.015](https://doi.org/10.1016/j.cpc.2018.01.015)
- [2019K] D. I. Ketcheson. Relaxation Runge–Kutta methods: Conservation and stability for
  inner-product norms. *SIAM J. Numer. Anal.* 57(6):2850–2870, 2019.
  [doi:10.1137/19M1263662](https://doi.org/10.1137/19M1263662)
- [2019KC] C. A. Kennedy, M. H. Carpenter. Diagonally implicit Runge–Kutta methods for stiff ODEs.
  *Appl. Numer. Math.* 146:221–244, 2019.
  [doi:10.1016/j.apnum.2019.07.008](https://doi.org/10.1016/j.apnum.2019.07.008)
- [2020KSSZ] D. I. Ketcheson, B. Seibold, D. Shirokoff, D. Zhou. DIRK schemes with high weak stage
  order. In *Spectral and High Order Methods for Partial Differential Equations ICOSAHOM 2018*, LNCSE
  134, pp. 453–463, Springer, 2020.
  [doi:10.1007/978-3-030-39647-3_36](https://doi.org/10.1007/978-3-030-39647-3_36)
- [2020M] The mathlib Community. The Lean mathematical library. In *CPP 2020*, pp. 367–381, ACM,
  2020. [doi:10.1145/3372885.3373824](https://doi.org/10.1145/3372885.3373824)
- [2020R+] H. Ranocha, M. Sayyari, L. Dalcin, M. Parsani, D. I. Ketcheson. Relaxation Runge–Kutta
  methods: Fully discrete explicit entropy-stable schemes for the compressible Euler and
  Navier–Stokes equations. *SIAM J. Sci. Comput.* 42(2):A612–A638, 2020.
  [doi:10.1137/19M1263480](https://doi.org/10.1137/19M1263480)
- [2020V+] P. Virtanen et al. SciPy 1.0: fundamental algorithms for scientific computing in
  Python. *Nat. Methods* 17:261–272, 2020.
  [doi:10.1038/s41592-019-0686-2](https://doi.org/10.1038/s41592-019-0686-2)
- [2021B+] C. J. Balos, D. J. Gardner, C. S. Woodward, D. R. Reynolds. Enabling GPU accelerated
  computing in the SUNDIALS time integration library. *Parallel Comput.* 108:102836, 2021.
  [doi:10.1016/j.parco.2021.102836](https://doi.org/10.1016/j.parco.2021.102836)
- [2021FKM] P. E. Farrell, R. C. Kirby, J. Marchena-Menéndez. Irksome: Automating Runge–Kutta
  time-stepping for finite element methods. *ACM Trans. Math. Softw.* 47(4), 2021.
  [doi:10.1145/3466168](https://doi.org/10.1145/3466168)
- [2021MU] L. de Moura, S. Ullrich. The Lean 4 theorem prover and programming language. In
  *CADE-28*, LNCS 12699, pp. 625–635, 2021.
  [doi:10.1007/978-3-030-79876-5_37](https://doi.org/10.1007/978-3-030-79876-5_37)
- [2021NRK] S. Nüßlein, H. Ranocha, D. I. Ketcheson. Positivity-preserving adaptive Runge–Kutta
  methods. *Communications in Applied Mathematics and Computational Science* 16(2):155–179, 2021.
  [doi:10.2140/camcos.2021.16.155](https://doi.org/10.2140/camcos.2021.16.155)
- [2021R+] M. M. Rana, V. E. Howle, K. Long, A. Meek, W. Milestone. A new block preconditioner for
  implicit Runge–Kutta methods for parabolic PDE problems. *SIAM J. Sci. Comput.* 43(5):S475–S495,
  2021. [doi:10.1137/20M1349680](https://doi.org/10.1137/20M1349680)
- [2021S+] K. Świrydowicz, J. Langou, S. Ananthan, U. Yang, S. Thomas. Low synchronization
  Gram–Schmidt and generalized minimal residual algorithms. *Numer. Linear Algebra Appl.*
  28(2):e2343, 2021. [doi:10.1002/nla.2343](https://doi.org/10.1002/nla.2343)
- [2022G+] D. J. Gardner, D. R. Reynolds, C. S. Woodward, C. J. Balos. Enabling new flexibility in
  the SUNDIALS suite of nonlinear and differential/algebraic equation solvers. *ACM Trans. Math.
  Softw.* 48(3), 2022. [doi:10.1145/3539801](https://doi.org/10.1145/3539801)
- [2022SKP] B. S. Southworth, O. Krzysik, W. Pazner. Fast solution of fully implicit Runge–Kutta
  and discontinuous Galerkin in time for numerical PDEs, Part II: nonlinearities and DAEs. *SIAM J.
  Sci. Comput.* 44(2):A636–A663, 2022.
  [doi:10.1137/21M1390438](https://doi.org/10.1137/21M1390438)
- [2022SKPD] B. S. Southworth, O. Krzysik, W. Pazner, H. De Sterck. Fast solution of fully implicit
  Runge–Kutta and discontinuous Galerkin in time for numerical PDEs, Part I: the linear setting.
  *SIAM J. Sci. Comput.* 44(1):A416–A443, 2022.
  [doi:10.1137/21M1389742](https://doi.org/10.1137/21M1389742)
- [2023AMF] R. Abu-Labdeh, S. MacLachlan, P. E. Farrell. Monolithic multigrid for implicit
  Runge–Kutta discretizations of incompressible fluid flow. *J. Comput. Phys.* 478:111961, 2023.
  [doi:10.1016/j.jcp.2023.111961](https://doi.org/10.1016/j.jcp.2023.111961)
- [2023BKSS] A. Biswas, D. I. Ketcheson, B. Seibold, D. Shirokoff. Design of DIRK schemes with high
  weak stage order. *Communications in Applied Mathematics and Computational Science* 18(1):1–28,
  2023. [doi:10.2140/camcos.2023.18.1](https://doi.org/10.2140/camcos.2023.18.1)
- [2023R+] D. R. Reynolds, D. J. Gardner, C. S. Woodward, R. Chinomona. ARKODE: A flexible IVP
  solver infrastructure for one-step methods. *ACM Trans. Math. Softw.* 49(2), 2023.
  [doi:10.1145/3594632](https://doi.org/10.1145/3594632)
- [2024ADN] O. Axelsson, I. Dravins, M. Neytcheva. Stage-parallel preconditioners for implicit
  Runge–Kutta methods of arbitrarily high order, linear problems. *Numer. Linear Algebra Appl.*
  31(1):e2532, 2024. [doi:10.1002/nla.2532](https://doi.org/10.1002/nla.2532)
- [2024GO] M. J. Gander, M. Outrata. Spectral analysis of implicit $s$-stage block Runge–Kutta
  preconditioners. *SIAM J. Sci. Comput.* 46(3):A2047–A2072, 2024.
  [doi:10.1137/23M1604266](https://doi.org/10.1137/23M1604266)
- [2024K] R. C. Kirby. On the convergence of monolithic multigrid for implicit Runge–Kutta time
  stepping of finite element problems. *SIAM J. Sci. Comput.* 46(5):S22–S45, 2024.
  [doi:10.1137/23M1569344](https://doi.org/10.1137/23M1569344)
- [2024M+] P. Munch, I. Dravins, M. Kronbichler, M. Neytcheva. Stage-parallel fully implicit
  Runge–Kutta implementations with optimal multilevel preconditioners at the scaling limit. *SIAM
  J. Sci. Comput.* 46(2):S71–S96, 2024.
  [doi:10.1137/22M1503270](https://doi.org/10.1137/22M1503270)
- [2025B+] C. J. Balos, M. Day, L. Esclapez, A. M. Felden, D. J. Gardner, M. Hassanaly, D. R.
  Reynolds, J. S. Rood, J. M. Sexton, N. T. Wimer, C. S. Woodward. SUNDIALS time integrators for
  exascale applications with many independent systems of ordinary differential equations. *Int. J.
  High Perform. Comput. Appl.* 39(1):123–146, 2025.
  [doi:10.1177/10943420241280060](https://doi.org/10.1177/10943420241280060)
- [2025CLGR] G. Čaklović, T. Lunet, S. Götschel, D. Ruprecht. Improving efficiency of parallel across
  the method spectral deferred corrections. *SIAM J. Sci. Comput.* 47(1):A430–A453, 2025.
  [doi:10.1137/24M1649800](https://doi.org/10.1137/24M1649800)
- [2025DKGK] I. Dravins, M. Koch, V. Griehl, K. Kormann. Performance evaluation of mixed-precision
  Runge–Kutta methods for the solution of partial differential equations. *Int. J. High Perform.
  Comput. Appl.* 40(4):487–507, 2025.
  [doi:10.1177/10943420251392963](https://doi.org/10.1177/10943420251392963)
- [2025ESR] S. Ekanathan, O. Smith, C. Rackauckas. A fully adaptive Radau method for the efficient
  solution of stiff ordinary differential equations at low tolerances. In *2025 IEEE High
  Performance Extreme Computing Conference (HPEC)*, pp. 1–9, 2025.
  [doi:10.1109/HPEC67600.2025.11196706](https://doi.org/10.1109/HPEC67600.2025.11196706)
- [2025KC] C. A. Kennedy, M. H. Carpenter. Singly-implicit Runge–Kutta methods for stiff, ordinary
  differential equations. NASA Technical Memorandum NASA/TM-20250008379, 2025 (not peer-reviewed).
  [ntrs.nasa.gov/citations/20250008379](https://ntrs.nasa.gov/citations/20250008379)
- [2025KM] R. C. Kirby, S. P. MacLachlan. Extending Irksome: improvements in automated Runge–Kutta
  time stepping for finite element methods. *ACM Trans. Math. Softw.* 51(3), 2025.
  [doi:10.1145/3759245](https://doi.org/10.1145/3759245)
- [2026DM] F. Durastante, M. Mazza. Stage-parallel implicit Runge–Kutta methods via low-rank matrix
  equation corrections. *J. Sci. Comput.* 107(1):4, 2026.
  [doi:10.1007/s10915-026-03185-z](https://doi.org/10.1007/s10915-026-03185-z)
- [2026LBMC] S. Leveque, L. Bergamaschi, Á. Martínez, E. Carson. Mixed precision solvers for the
  all-at-once Runge–Kutta discretization of the heat equation. arXiv:2609.16863, 2026 (preprint, not
  peer-reviewed). [arxiv.org/abs/2609.16863](https://arxiv.org/abs/2609.16863)

Software, cited for Table 6.1 (documentation and sources, not peer-reviewed):

- RADAU5 and RADAU: <https://www.unige.ch/~hairer/software.html>
- SciPy `Radau`: <https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.Radau.html>
- OrdinaryDiffEq.jl FIRK solvers: <https://docs.sciml.ai/OrdinaryDiffEq/stable/implicit/FIRK/>
- PETSc `TSIRK`: <https://petsc.org/release/manualpages/TS/TSIRK/>
- Irksome: <https://github.com/firedrakeproject/Irksome>
- PSIDE (CWI reports MAS-R9833 and MAS-R9834): <https://ir.cwi.nl/pub/4595>, <https://ir.cwi.nl/pub/4594>
- accelerInt: <https://github.com/SLACKHA/accelerInt>
- Bari IVP test set, release 2.4 (2008): <https://archimede.uniba.it/~testset/testsetivpsolvers/>
- Libraries without Radau IIA:
  [deal.II](https://www.dealii.org/current/doxygen/deal.II/classTimeStepping_1_1ImplicitRungeKutta.html),
  [MFEM](https://docs.mfem.org/html/classmfem_1_1ODESolver.html),
  [Trilinos Tempus](https://github.com/trilinos/Trilinos/tree/master/packages/tempus),
  [GSL](https://www.gnu.org/software/gsl/doc/html/ode-initval.html),
  [Boost.Odeint](https://www.boost.org/doc/libs/1_86_0/libs/numeric/odeint/doc/html/boost_numeric_odeint/odeint_in_detail/steppers.html),
  [diffrax](https://docs.kidger.site/diffrax/api/solvers/ode_solvers/),
  [Drake](https://drake.mit.edu/doxygen_cxx/classdrake_1_1systems_1_1_radau_integrator.html)

## Appendix A: stability plots

Figures 1–4 are drawn by `scripts/firkode_stability_plot.py` for the Radau IIA methods with
$1 \le s \le 9$. The script evaluates $R(z) = \frac{P(z)}{Q(z)}$, the $(s-1, s)$ Padé approximant of
$e^z$, from the exact rational coefficients

$$
P(z) = \sum_{j=0}^{s-1} \frac{(2s-1-j)!(s-1)!}{(2s-1)!j!(s-1-j)!} z^j,
\qquad
Q(z) = \sum_{j=0}^{s} \frac{(2s-1-j)!s!}{(2s-1)!j!(s-j)!} (-z)^j,
$$

and checks in 50-digit arithmetic that it agrees with $R(z) = 1 + z b^T (I - zA)^{-1} \mathbf{1}$
for the tables of the generator `scripts/firkode_radau_tables.py`, the reference for the C tables,
to $10^{-30}$ at four sample points; the script therefore requires the generator of Milestone M0. The
poles of $R$ are the eigenvalues of $A^{-1}$, and $\gamma_0$ follows the rule of Section 4.5; the
printed values are the literals of `firkode_tables.c` to all printed digits. For the Dahlquist
equation $\dot y = \lambda y$ with $z = h\lambda$, the stage increments are
$Z = \bigl((I - zA)^{-1} - I\bigr) \mathbf{1} y_n$, and the filtered estimate of Section 4.6 is

$$
\frac{\mathrm{err}}{y_n}
= \frac{\gamma_0 z + e^T \bigl((I - zA)^{-1} - I\bigr) \mathbf{1}}{1 - \gamma_0 z},
\qquad
e^T = (\hat b - b)^T A^{-1} = -\gamma_0 l(0)^T A^{-1},
$$

with $l_i(0)$ the Lagrange basis on the nodes $c$ evaluated at $0$. It is $O(z^{s+1})$ as $z \to 0$ and
tends to $-1$ as $|z| \to \infty$. Its numerator cancels $O(z)$ terms down to $O(z^{s+1})$, so it is
evaluated with `mpmath`. The order stars are evaluated in double precision except on the disk
$|z| < 2\left(\frac{u}{C_s}\right)^{\frac{1}{2s}}$, with $C_s = \frac{(s-1)! s!}{(2s-1)!(2s)!}$ the Padé error constant,
$e^z - R(z) = C_s z^{2s} + \dots$, where double precision cannot separate $|R(z)|$ from $|e^z|$; those
grid points are recomputed with `mpmath` (radius $3.6$ for $s = 9$).

Run from the repository root in an environment with `numpy`, `matplotlib` and `mpmath`:

```bash
python3 scripts/firkode_stability_plot.py
```

It writes the four PNG files to `doc/shared/figs/firkode/` and prints the versions it used. The
figures below were produced on October 1, 2026 with Python 3.14.7, numpy 2.5.3, matplotlib 3.11.2 and
mpmath 1.4.1. Every built-in check passed: Padé and table stability functions agree, odd $s$ have
exactly one real eigenvalue of $A^{-1}$, $\max_y |R(iy)| \le 1$, the filtered estimate
tends to $-1$, and the slope of the estimate for small $z$ on a log-log scale is $s + 1$ within
$0.01$ for every $s$.

![Stability regions of the Radau IIA methods](../../doc/shared/figs/firkode/RADAU_IIA_stability_region.png)

*Figure 1. The curves* $|R(z)| = 1$ *for* $s = 1, \dots, 9$. *Each method is stable outside its
curve, so every stability region contains the closed left half-plane. The bounded shaded sets*
$|R(z)| > 1$ *lie in* $\mathrm{Re} z > 0$ *and grow with* $s$; *for* $s = 9$ *the set reaches*
$\mathrm{Re} z = 75.1$ *and* $\lvert \mathrm{Im} z \rvert = 43.6$.

![Order stars of the Radau IIA methods](../../doc/shared/figs/firkode/RADAU_IIA_order_star.png)

*Figure 2. Order stars* $\{z : |R(z)| > |e^z|\}$ *[1978WHN]. Order* $2s - 1$ *gives* $2s$ *shaded
and* $2s$ *white sectors of angle* $\frac{\pi}{2s}$ *at the origin. Each bounded shaded finger
contains one pole of* $R$ *(circles) and each bounded white finger one zero (crosses). The shaded
set meets the imaginary axis only at the origin and* $R$ *has no pole in the left half-plane, which
is A-stability;* $\deg P < \deg Q$ *gives* $R(\infty) = 0$, *L-stability.*

![Stability function along the axes](../../doc/shared/figs/firkode/RADAU_IIA_axis_profile.png)

*Figure 3.* $|R(iy)|$ *(left) and* $|R(x)|$ *for* $x < 0$ *(right).* $|R(iy)| \le 1$ *with equality
only at* $y = 0$, *and* $|R(z)| \sim \frac{s}{|z|}$ *as* $|z| \to \infty$, *so the damping of stiff
components per step weakens as* $s$ *grows. The dips on the right are the real zeros of* $P$, *which
exist for even* $s$ *only.*

![Filtered error estimate along the axes](../../doc/shared/figs/firkode/RADAU_IIA_error_estimate.png)

*Figure 4. The filtered error estimate* $\left|\frac{\mathrm{err}}{y_n}\right|$ *for*
$\dot y = \lambda y$ *along the negative real axis (left) and the imaginary axis (right). The
slopes for small* $z$ *are* $s + 1$, *the estimate having order* $p = s$. *As* $|z| \to \infty$ *the
estimate tends to* $-y_n$ *instead of growing like* $\gamma_0 z y_n$, *which the unfiltered
estimate (dashed,* $s = 3$) *does.*

## Appendix B: computation of Table 4.1

Table 4.1 is `[derived]` by the following procedure, which the generator reproduces with `--table1`
and which needs only exact rational arithmetic and a root finder.

1. The eigenvalues $\lambda_j$ of $A^{-1}$ are the roots of $\det(I - zA)$, which equals the Padé
   denominator $Q(z)$ of Appendix A `[proved]`. They are found as the roots of $Q$ (Durand–Kerner
   iteration or any polynomial root finder); the eigenvalues of $A$ are $\mu_j = \frac{1}{\lambda_j}$.
2. $\gamma_0 = \frac{1}{\lambda_r}$ for odd $s$, with $\lambda_r$ the unique real root; for even $s$,
   $\gamma_0 = \frac{1}{\mathrm{Re} \lambda_\ast}$ with $\lambda_\ast$ the root of smallest modulus
   (Section 4.5).
3. The spectrum of $K P^{-1}$ on the model problem is the union of the disks $D_j$, the images of the
   closed left half-plane under $w_j(z) = \frac{1 - z\mu_j}{1 - \gamma_0 z}$. By the maximum principle
   it suffices to sample the boundary circles, the images of the imaginary axis
   $z = i \tan\Bigl(\frac{\pi (k + \frac{1}{2})}{n} - \frac{\pi}{2}\Bigr)$, $k = 0, \dots, n-1$, together with the stiff-limit points
   $\frac{\mu_j}{\gamma_0}$; $n = 4000$ samples suffice for three digits.
4. $\rho_\ast(\gamma_0) = \min_{c > 0} \max_j \max_{w \in \partial D_j} \Bigl\lvert 1 - \frac{w}{c} \Bigr\rvert$ is
   found by golden-section search on $c \in [0.05, 5]$, on which the function is unimodal in
   practice; a bracketing scan confirms the minimum.
5. $k = \Bigl\lceil \frac{\ln 0.05}{\ln \rho_\ast} \Bigr\rceil$ is the iteration count that guarantees a 20-fold
   residual reduction when $\kappa = 1$.

For $s = 1$ the single disk degenerates to the point $1$ and $\rho_\ast = 0$: the preconditioner is
exact. The $\gamma_0$ values of step 2 must agree with the generator's literals to every digit the
procedure prints; `--table1` performs that comparison as one of its checks.
