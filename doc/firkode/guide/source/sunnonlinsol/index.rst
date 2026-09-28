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

.. _SUNNonlinSol:

###########################
Nonlinear Algebraic Solvers
###########################

SUNDIALS time integration packages are written in terms of generic nonlinear
solver operations defined by the SUNNonlinSol API and implemented by a
particular SUNNonlinSol module of type ``SUNNonlinearSolver``.
Users can supply their own SUNNonlinSol module, or use one of the modules
provided with SUNDIALS. Depending on the package, nonlinear solver modules
can either target systems presented in a rootfinding (:math:`F(y) = 0`) or
fixed-point (:math:`G(y) = y`) formulation. FIRKODE solves the coupled
stage equations with the SUNNONLINSOL_NEWTON module internally and does not
currently accept a user-supplied nonlinear solver; this chapter is included
for reference on the nonlinear solver API and the module used by FIRKODE
(see :numref:`FIRKODE.Mathematics.Newton`).




For users interested in providing their own SUNNonlinSol module, the
following section presents the SUNNonlinSol API and its implementation
beginning with the definition of SUNNonlinSol functions in the
sections :numref:`SUNNonlinSol.API.CoreFn`, :numref:`SUNNonlinSol.API.SetFn` and
:numref:`SUNNonlinSol.API.GetFn`. This is followed by the definition of
functions supplied to a nonlinear solver implementation in the section
:numref:`SUNNonlinSol.API.SUNSuppliedFn`.  The nonlinear solver return
codes are given in the section :numref:`SUNNonlinSol.API.ReturnCodes`. The
``SUNNonlinearSolver`` type and the generic SUNNonlinSol module are defined
in the section :numref:`SUNNonlinSol.API.Generic`. Finally, the section
:numref:`SUNNonlinSol.API.Custom` lists the requirements for supplying a custom
SUNNonlinSol module. Users wishing to supply their own SUNNonlinSol module
are encouraged to use the SUNNonlinSol implementations provided with
SUNDIALS as templates for supplying custom nonlinear solver modules.




.. toctree::
   :maxdepth: 1

   SUNNonlinSol_API_link.rst
   SUNNonlinSol_links.rst
