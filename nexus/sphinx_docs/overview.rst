.. _overview:

Overview of Nexus
=================

What Nexus is
-------------

Nexus is a collection of tools, written in Python, to perform
complex electronic structure calculations and analyze the results.  The main
focus is currently on performing arbitrary Quantum Monte Carlo (QMC)
calculations with QMCPACK; however, VASP, Quantum Espresso, and GAMESS are
also supported.  A single QMC calculation typically requires several
previous calculations with other codes to produce a starting guess for the
many-body wavefunction and convert it into a form that QMCPACK understands.
Managing the resulting array of calculations, and the flow of information
between them, quickly becomes unwieldy to the researcher, demands a great
deal of human time, and increases the potential for human error.  Nexus
reduces both the human time required and potential for error by
automating the total simulation process.

What Nexus can do
-----------------

The capabilities of Nexus currently include crystal structure
generation, standalone Density Functional Theory (DFT) calculations with
PWSCF (Quantum Espresso) or VASP,  quantum chemical calculations with GAMESS,
Hartree-Fock (HF) calculations of atoms with the SQD code (packaged with
QMCPACK), complete QMC calculations with QMCPACK (including wavefunction
optimization, Variational Monte Carlo (VMC), and Diffusion Monte Carlo (DMC) in
periodic or open boundary conditions), automated job management on workstations
(by acting as a virtual queue) and clusters/supercomputers
including handling of dependencies
between calculations and job bundling,  and extraction of results from
completed calculations for analysis.  The integration of these capabilities
permits the user to focus on the high-level tasks of problem formulation and
interpretation of the results without (in principle) becoming too involved
in the time-consuming, lower level details.

How Nexus is used
-----------------

Use of Nexus currently involves writing a short Python script
describing the calculations to be performed.  This small script formed by the
user closely resembles an input file for electronic structure codes.  A key
difference is that this "input file" represents executable code, and so
variables are easily defined for use in expressions and more complicated
simulation workflows (*e.g.* an equation of state) can be constructed
with if/else logic and for loops.  Knowledge of the Python programming language
is helpful to perform complex calculations, but not essential for use of
Nexus.  Starting from working "input files" such as those covered
on the :ref:`examples` page is a good way to proceed.
