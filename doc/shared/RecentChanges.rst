.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

Added a new package, FIRKODE, that integrates stiff ODE systems
:math:`M y' = f(t,y)` with the fully implicit Radau IIA Runge-Kutta methods (orders
1 to 17, with the classical three-stage, fifth-order method of RADAU5 as the
default). The coupled stage equations are solved with a simplified Newton
iteration on the stacked stage system whose linear systems are solved by an
internal flexible GMRES iteration preconditioned block-wise with any SUNDIALS
``SUNLinearSolver`` (direct, matrix-based iterative, or matrix-free), so that
only systems of the size of the ODE are ever factored. FIRKODE provides the
RADAU5 filtered error estimate, ``SUNAdaptController``-based step size
adaptivity, collocation dense output, rootfinding, and constant non-identity
mass matrices. The package is enabled with the :cmakeop:`SUNDIALS_ENABLE_FIRKODE` CMake
option and is exported as the ``SUNDIALS::firkode`` target. Fortran and Python
interfaces are not yet available.

**New Features and Enhancements**

Added the utility function, :c:func:`SUNFileFlush` for flushing file
pointers. This is useful when using the Fortran 2003 interfaces.

**Bug Fixes**

Fixed a bug in ``FindMAGMA.cmake`` which didn't allow use of MAGMA versions with 
multiple digits in an identifier.

**Deprecation Notices**
