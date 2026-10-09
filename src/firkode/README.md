# FIRKODE
### Version 0.1.0 (Sep 2026)

**Yifan Hu, Department of Mathematics and Statistics, UMBC**

FIRKODE is a package for the solution of stiff ordinary differential equation
(ODE) systems (initial value problems) given in the form
```
M y' = f(t,y), y(t0) = y0,
```
using fully implicit Runge-Kutta methods. The integration methods of FIRKODE
are the Radau IIA collocation methods (`s` stages, order `2s-1`, stiffly
accurate, L-stable); version 0.1.0, experimental and numbered independently of
SUNDIALS, provides the three-stage method of order five, the method of RADAU5,
and the other stage counts are planned for later versions. The coupled stage
equations are solved with a simplified Newton iteration whose linear systems
are solved by a block-preconditioned flexible GMRES iteration that reuses a
single user-supplied SUNDIALS linear solver for the shifted systems
`M - h*gamma*J`.

FIRKODE is part of the SUNDIALS Suite of Nonlinear and Differential/Algebraic
equation Solvers which consists of ARKODE, CVODE, CVODES, FIRKODE, IDA, IDAS,
and KINSOL. It is written in ANSI standard C and runs with the serial,
parallel and many-vector `N_Vector` implementations and every SUNDIALS linear
solver, through the shared vector, matrix and linear solver APIs used across
SUNDIALS packages; it has not yet been exercised on accelerator-based (e.g.,
GPU) systems.

## Documentation

See the FIRKODE documentation at [Read the Docs](https://sundials.readthedocs.io/en/latest/firkode)
for more information about FIRKODE usage.

## Installation

For installation instructions see the
[SUNDIALS Installation Guide](https://sundials.readthedocs.io/en/latest/Install_link.html).

## Release History

Information on recent changes to FIRKODE can be found in the "Introduction"
chapter of the FIRKODE User Guide and a complete release history is available in
the "SUNDIALS Release History" appendix of the FIRKODE User Guide.

## References

* Y. Hu, "User Documentation for FIRKODE v0.1.0," Technical Report,
  UMBC, Sep 2026.

* E. Hairer and G. Wanner, "Solving Ordinary Differential Equations II: Stiff
  and Differential-Algebraic Problems," 2nd ed., Springer, 1996.
