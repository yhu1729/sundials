.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

Added a new, experimental package, FIRKODE, that integrates stiff ODE systems
:math:`M y' = f(t,y)` with the fully implicit Radau IIA Runge-Kutta methods. Its
version 0.1.0, numbered independently of SUNDIALS, provides the three-stage
method of order five, the method of RADAU5, with its simplified Newton
iteration, filtered error estimate, predictive step size controller, Newton
tolerance and Jacobian reuse policy; the other stage counts (orders 1 to 17)
are planned for later versions. The coupled stage equations are solved on the
stacked stage system by an internal flexible GMRES iteration preconditioned
block-wise with any SUNDIALS ``SUNLinearSolver`` (direct, matrix-based
iterative, or matrix-free), so that only systems of the size of the ODE are
ever factored. FIRKODE provides ``SUNAdaptController``-based step size
adaptivity, collocation dense output, rootfinding, and constant non-identity
mass matrices, and runs with the serial, parallel and many-vector ``N_Vector``
implementations. The package is enabled with the
:cmakeop:`SUNDIALS_ENABLE_FIRKODE` CMake option and is exported as the
``SUNDIALS::firkode`` target. Its API may change before version 1.0; Fortran
and Python interfaces are not yet available.

**New Features and Enhancements**

Added the utility function, :c:func:`SUNFileFlush` for flushing file
pointers. This is useful when using the Fortran 2003 interfaces.

**Bug Fixes**

Fixed a segfault that could occur when handling errors without a ``SUNContext``
provided.

Fixed a bug in ``FindMAGMA.cmake`` which didn't allow use of MAGMA versions with 
multiple digits in an identifier.

Fixed a bug in the LSRKStep module where the step post-processing function for
the RKL step was called with an incorrect time argument.

**Deprecation Notices**
