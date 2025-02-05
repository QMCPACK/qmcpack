.. _examples:

QMCPACK Examples
================

**WARNING: TO RUN QUICKLY ENOUGH, MOST OF THESE EXAMPLES ARE NOT CONVERGED TO RESEARCH STANDARDS! YOU MUST CONVERGE
PARAMETERS (STATISTICS, SIMULATION CELL SIZE, JASTROW PARAMETER NUMBER/CUTOFF, TWIST NUMBER, DMC TIME STEP,
DFT PLANE WAVE CUTOFF, DFT K-POINT MESH, ETC.) FOR REAL CALCULATIONS!**

The following examples should run on a modern workstation in a few hours.

Using QMCPACK directly
----------------------

In ``examples/molecules`` are the following examples.
Each directory also contains a ``README`` file with more details.

========= =====================================
Directory Description
``H20``   H2O molecule from GAMESS orbitals
``He``    Helium atom with simple wavefunctions
========= =====================================

In ``examples/solid`` is a simple bulk solid LiH example, together with
example Quantum ESPRESSO, PySCF, and RMG inputs to generate the orbitals.
For research calculations on solids it is strongly encouraged to use Nexus
for which many more examples are provided and e.g. twist averaging can be
performed automatically.

Using Nexus
-----------

Extensive Nexus workflow examples are provided including (1) demonstrations of the use of Quantum ESPRESSO, PySCF, and Quantum Package to generate trial wavefunctions and perform subsequent QMC calculations,
(2) application to isolated molecules and periodic systems with both plane-waves/splines and Gaussians, including periodic Gaussians,
and (3) use of different estimators and analysis of the samples results.

For information about using Nexus, see the User Guide in ``nexus/docs`` or at https://nexus-workflows.readthedocs.io/ .

For Python to find the Nexus library, the PYTHONPATH environment variable should be set to ``<QMCPACK source>/nexus/lib``.
The executables for both QMCPACK and Quantum ESPRESSO (or relevant density functional or quantum chemical code) should be
on the PATH, or the paths in the script should be adjusted.

These QMCPACK Nexus examples can be found under the ``nexus/examples/qmcpack`` directory.

======================================================= =================================================================================================================================================================
Directory                                               Description
``rsqmc_quantum_espresso/01_diamond_dft_vmc``           Carbon diamond primitive cell, run with pseudopotentials and spline orbitabls from Quantum ESPRESSO, then VMC optimization with 1 and 2 body Jastrow factors.
``rsqmc_quantum_espresso/02_diamond_dft_vmc_twistavg``  As above, but including automatic generation of twist averaging over multiple twists.
``rsqmc_quantum_espresso/03_diamond_dft_dmc_textrap``   Carbon diamond example including generation of DMC runs for time-step extrapolation.
``rsqmc_quantum_espresso/04_iron_dft_dmc_gcta``         Full grand-canonical twist averaging (GCTA) example for bulk Fe including VMC with 1,2 and 3 body Jastrow factors and spin-density accumulation.
``rsqmc_pyscf/01_h2o_hf_qmc``                           H\ :math:`_2`\ O molecule, run all electron using Hartree-Fock orbitals from PySCF, then VMC optimization with 1,2 body Jastrow factors, then final DMC.
``rsqmc_pyscf/02_diamond_hf_qmc``                       Periodic carbon diamond primitive cell, run with pseudopotentials and Quantum ESPRESSO, then VMC optimization with 1 and 2 body Jastrow factors, then final DMC.
``rsqmc_quantum_package/01_h2o_hf_qmc``                 H\ :math:`_2`\ O molecule, run all electron using Hartree-Fock orbitals from Quantum Package, then VMC optimization with 1,2 body Jastrow factors, then final DMC.
``rsqmc_quantum_package/02_o2_selci_qmc``               O\ :math:`_2` molecule using 5000 determinants via configuration interaction, followed by QMC as above.
``rsqmc_misc/H2O``                                      H\ :math:`_2`\ O molecule run with pseudopotentials and spline orbitabls from Quantum ESPRESSO, followed by VMC optimization and DMC.
``rsqmc_misc/H2O_pyscf``                                H\ :math:`_2`\ O molecule using PySCF. All electron then single VMC determinant only calculation.
``rsqmc_misc/diamond``                                  Carbon diamond supercell example including relaxation of atomic positions using Quantum ESPRESSO.
``rsqmc_misc/diamond_pyscf``                            Carbon diamond primitive cell examples run with pseudopotentials and spline orbitals from Quantum ESPRESSO.
``rsqmc_misc/diamond_radial_density``                   Spin density analysis example with computation of radial densities using qdens.
``rsqmc_misc/diamond_lowdin``                           As above, but Lowdin analysis.
``rsqmc_misc/oxygen_dimer``                             Binding curve for O\ :math:`_2` molecule using pseudopotentials and spline orbitals from Quantum ESPRESSO.
``rsqmc_misc/LiH``                                      Bulk LiH primitive cell VMC and DMC example using BFD and CASINO format pseudopotentials and orbitals from Quantum ESPRESSO.
``rsqmc_misc/excited``                                  Excited state example for 16 atom supercell carbon diamond.
``rsqmc_misc/O2_qp``                                    O\ :math:`_2` molecule using 5000 determinants via configuration interaction and Quantum Package, followed by QMC.
``rsqmc_misc/estimators``                               Complete set of estimator specification examples applied for bulk Fe using DFT+U based orbitals.
``rsqmc_misc/graphene``                                 Graphene sheet DMC example including use of Nexus analyzer to obtain total energy.
``rsqmc_misc/c20``                                      C\ :math:`_{20}` fullerene molecule using pseudopotentials and spline orbitals from Quantum ESPRESSO.
======================================================= =================================================================================================================================================================
