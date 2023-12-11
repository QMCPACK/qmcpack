.. _examples:

Examples
========

**WARNING: THESE EXAMPLES ARE NOT CONVERGED! YOU MUST CONVERGE
PARAMETERS (SIMULATION CELL SIZE, JASTROW PARAMETER NUMBER/CUTOFF, TWIST
NUMBER, DMC TIME STEP, DFT PLANE WAVE CUTOFF, DFT K-POINT MESH, ETC.)
FOR REAL CALCUATIONS!**

The following examples should run in serial on a modern workstation in a few hours.

Using QMCPACK directly
----------------------

In ``examples/molecules`` are the following examples.
Each directory also contains a ``README`` file with more details.

========= =====================================
Directory Description
``H20``   H2O molecule from GAMESS orbitals
``He``    Helium atom with simple wavefunctions
========= =====================================

Using Nexus
-----------

For more information about Nexus, see the User Guide in ``nexus/documentation``.

For Python to find the Nexus library, the PYTHONPATH environment variable should be set to ``<QMCPACK source>/nexus/library``.
For these examples to work properly, the executables for QE and QMCPACK either
need to be on the path, or the paths in the script should be adjusted.

These examples can be found under the ``nexus/examples/qmcpack`` directory.

================ ==========================================
Directory        Description
``diamond``      Bulk diamond with VMC
``graphene``     Graphene sheet with DMC
``c20``          C20 cage molecule
``oxygen_dimer`` Binding curve for O\ :math:`_2` molecule
``H2O``          H\ :math:`_2`\ O molecule with QE orbitals
``LiH``          LiH crystal with QE orbitals
================ ==========================================
